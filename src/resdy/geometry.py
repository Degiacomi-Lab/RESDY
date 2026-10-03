'''
Geometry check applied to a structure once curation has finished.

Curation rebuilds missing residues, splits and reassembles chains and renumbers
residues, and every one of those steps has at some point moved atoms that should not
have moved. None of it is visible in the measurements themselves, which is what makes a
geometric check worth having: a structure whose chains have drifted or whose rebuilt loop
runs into a neighbour still measures perfectly happily.

The module also holds :func:`distance_to_other_chains`, which records how close a residue
sits to a neighbouring chain, the metadata column Min_Dist_Other_Chain.
'''

import numpy as np
import pandas as pd
import biobox as bb
from .helper import require_biobox
require_biobox(bb)
from scipy.spatial import cKDTree
from .residues import PROTEIN_RESNAMES


#: Residue names treated as solvent.
WATER_RESNAMES = ('HOH', 'WAT', 'DOD', 'H2O')

#: Monoatomic ions. Their contacts with protein are coordination bonds, not clashes.
ION_RESNAMES = ('ZN', 'NI', 'CU', 'FE', 'MG', 'MN', 'NA', 'K', 'CA', 'CO', 'CL', 'MO',
                'CD', 'HG', 'PT', 'AU', 'AG', 'SR', 'BA', 'CS', 'RB', 'LI', 'BR', 'IOD', 'F')

#: The 20 standard residue names, used to tell protein from everything else.
STANDARD_RESIDUES = ('ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
                     'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL')

#: Metal elements. A short contact to one of these is a coordination bond, not a clash,
#: whether the metal is a free ion or an atom inside a cofactor such as the haem iron.
METAL_ELEMENTS = ('FE', 'ZN', 'CU', 'MG', 'MN', 'NI', 'CO', 'MO', 'CA', 'NA', 'K', 'LI',
                  'CD', 'HG', 'PT', 'AU', 'AG', 'W', 'V', 'CR', 'SR', 'BA', 'CS', 'RB')

def _read_altloc_occupancy(path, n_expected):
    '''
    Alternate-location indicator and occupancy of every atom record, in file order.

    Both are read from the text rather than from biobox, because biobox folds columns 12
    to 17 into its ``name`` field, so the altloc ends up glued to the atom name ('OE2A')
    and cannot be recovered from it unambiguously. Only the first model is read, which is
    what biobox loads.

    :returns: (altlocs, occupancies). If the record count does not match what biobox
        loaded, blanks and ones are returned so that no pair is excluded on this basis.
    '''
    alt, occ = [], []
    with open(path) as fh:
        for line in fh:
            if line.startswith('ENDMDL'):
                break
            if line.startswith(('ATOM', 'HETATM')) and len(line) > 16:
                alt.append(line[16])
                try:
                    occ.append(float(line[54:60]))
                except ValueError:
                    occ.append(1.0)
    if len(alt) != n_expected:
        print(f'>> could not align atom records for {path} '
              f'({len(alt)} in file, {n_expected} loaded); alternate conformers and '
              f'partial occupancies not excluded')
        return np.full(n_expected, ' '), np.ones(n_expected)

    return np.array(alt), np.array(occ)


def _load(path):
    '''
    Read a structure and return its heavy atoms.

    :returns: (atom records, (n, 3) coordinates, element symbols, altlocs, occupancies)
    '''
    M = bb.Molecule()
    M.import_pdb(path, include_hetatm=True)
    d = M.data.reset_index(drop=True)
    xyz = np.asarray(M.points, dtype=float)
    alt, occ = _read_altloc_occupancy(path, len(d))

    elem = d['atomtype'].astype(str).str.upper().str.strip()
    heavy = (~elem.isin(['H', 'D'])).to_numpy()
    return (d[heavy].reset_index(drop=True), xyz[heavy], elem[heavy].to_numpy(),
            alt[heavy], occ[heavy])


def _residue_ordinals(d):
    '''
    Position of each atom's residue within its chain, counted in file order.

    Sequence adjacency has to be decided this way rather than by residue number: after
    curation the numbering follows the Uniprot sequence and skips wherever the structure
    does, so two residues joined by a peptide bond can differ by more than one.

    :returns: (n,) array of ordinals.
    '''
    chain = d['chain'].astype(str).to_numpy()
    resid = d['resid'].astype(int).to_numpy()
    ordinal = np.empty(len(d), dtype=int)
    seen = {}
    for pos, (c, r) in enumerate(zip(chain, resid)):
        key = (c, r)
        if key not in seen:
            seen[key] = len([k for k in seen if k[0] == c])
        ordinal[pos] = seen[key]
    return ordinal


def _offending_pairs(d, xyz, elem, alt, occ, cutoff, min_occupancy=0.0):
    '''
    Atom pairs closer than ``cutoff`` that are not bonded, not in the same residue, and
    do not involve a monoatomic ion.

    Bonds are recognised from connectivity rather than from distance, because a peptide
    bond at 1.3 A and a severe clash at 1.3 A are not distinguishable geometrically. In
    a protein the only bonds that cross a residue boundary are the backbone C-N, the
    disulfide SG-SG, and any link into a cofactor; the first two are excluded by name
    and the third is reported with its own flag.

    :returns: (i, j, separations), all numpy arrays sorted by separation.
    '''
    pairs = cKDTree(xyz).query_pairs(cutoff, output_type='ndarray')
    if len(pairs) == 0:
        return np.array([], int), np.array([], int), np.array([], float)

    i, j = pairs[:, 0], pairs[:, 1]
    chain = d['chain'].astype(str).to_numpy()
    resid = d['resid'].astype(int).to_numpy()
    resname = d['resname'].astype(str).str.upper().str.strip().to_numpy()
    name = d['name'].astype(str).str.upper().str.strip().to_numpy()
    is_metal = np.isin(resname, ION_RESNAMES) | np.isin(elem, METAL_ELEMENTS)
    ordinal = _residue_ordinals(d)

    BACKBONE = ('N', 'CA', 'C', 'O')
    same_chain = chain[i] == chain[j]
    adjacent = same_chain & (np.abs(ordinal[i] - ordinal[j]) == 1)
    # Backbone atoms of sequence-adjacent residues are one or two bonds apart across the
    # peptide bond: C-N is the bond itself, and CA-N, C-CA and O-N are the 1-3 pairs it
    # creates. None of them is a steric contact, and the 1-3 pairs sit around 2.4 A, well
    # inside any clash cutoff.
    backbone_pair = np.isin(name[i], BACKBONE) & np.isin(name[j], BACKBONE)

    keep = ~(
        (same_chain & (resid[i] == resid[j]))        # atoms of one residue
        | (adjacent & backbone_pair)                 # the peptide bond and its 1-3 pairs
        | ((name[i] == 'SG') & (name[j] == 'SG'))    # disulfide
        | is_metal[i] | is_metal[j]                  # metal coordination
        | ((alt[i] != ' ') & (alt[j] != ' ') & (alt[i] != alt[j]))   # rival conformers
        | (occ[i] < min_occupancy) | (occ[j] < min_occupancy)        # partial occupancy
    )
    i, j = i[keep], j[keep]
    sep = np.linalg.norm(xyz[i] - xyz[j], axis=1)
    order = np.argsort(sep)
    return i[order], j[order], sep[order]


def distance_to_other_chains(chain, resname, element, xyz, anchor,
                             protein_resnames=PROTEIN_RESNAMES):
    '''
    Distance from each anchor atom to the closest protein heavy atom of any other chain.

    Only heavy atoms of protein residues count as a partner: an ion, a water or a ligand
    that carries a chain identifier of its own is not a chain. Which chain is the closest is
    not reported.

    :param chain: (n,) chain identifier of every atom.
    :param resname: (n,) residue name of every atom.
    :param element: (n,) element symbol of every atom, used to leave out hydrogens.
    :param xyz: (n, 3) coordinates.
    :param anchor: (n,) boolean mask of the atoms to measure from.
    :param protein_resnames: residue names counted as protein.
    :returns: (number of anchors,) distances in A, in the order of ``np.flatnonzero(anchor)``,
        NaN for an anchor whose structure holds no other protein chain.
    :rtype: numpy.ndarray
    '''
    chain = np.asarray(chain).astype(str)
    resname = np.char.upper(np.char.strip(np.asarray(resname).astype(str)))
    element = np.char.upper(np.char.strip(np.asarray(element).astype(str)))
    xyz = np.asarray(xyz, dtype=float)
    anchor = np.asarray(anchor, dtype=bool)

    partner = np.isin(resname, list(protein_resnames)) & ~np.isin(element, ['H', 'D'])
    anchor_idx = np.flatnonzero(anchor)
    out = np.full(len(anchor_idx), np.nan)
    for c in np.unique(chain[anchor_idx]):
        others = partner & (chain != c)
        if not others.any():
            continue
        mine = chain[anchor_idx] == c
        out[mine] = cKDTree(xyz[others]).query(xyz[anchor_idx[mine]])[0]
    return out


def check_geometry(path, clash_cutoff=2.0, modelled_residues=(), report_cutoff=4.0,
                   min_occupancy=0.0):
    '''
    Find heavy-atom contacts closer than ``clash_cutoff`` between atoms that are not
    bonded neighbours.

    Pairs excluded as legitimately close: atoms of the same residue; backbone atoms of
    sequence-adjacent residues, which covers the peptide bond and the 1-3 pairs it
    creates; disulfide SG-SG; and any pair involving a metal, since coordination is a
    bond. The metal test is on the element, so it covers a cofactor's metal centre as
    well as a free ion.

    Bonds are recognised from connectivity rather than from distance. A peptide bond at
    1.3 A and a severe clash at 1.3 A are not distinguishable geometrically, so any rule
    based on covalent radii either reports real bonds or excuses real overlaps.

    :param path: curated pdb file to check.
    :type path: str
    :param clash_cutoff: heavy-atom separation below which a contact is a clash, in A.
    :type clash_cutoff: float
    :param modelled_residues: (chain, resid) pairs built by Modeller, so that clashes
        involving rebuilt atoms can be reported separately. Atoms with an occupancy of
        zero, which is how curation marks the atoms Modeller built, count as modelled too.
    :type modelled_residues: iterable
    :param report_cutoff: when no clash is found, the radius searched to report the
        closest non-bonded contact anyway, in A.
    :type report_cutoff: float
    :param min_occupancy: contacts involving an atom with an occupancy below this value
        are not reported. Defaults to 0, reporting every contact: a curated file holds a
        single conformer, so its partially occupied atoms are atoms of that conformer, and
        the atoms Modeller built carry an occupancy of zero and are the ones most likely to
        clash.
    :type min_occupancy: float
    :returns: (summary dict, DataFrame of offending pairs, worst first)
    :rtype: tuple
    '''
    d, xyz, elem, alt, occ = _load(path)

    summary = {'file': path, 'n_atoms': int(len(d)), 'min_contact': np.nan,
               'n_clashes': 0, 'worst': '', 'involves_modelled': False,
               'involves_hetatm': False, 'involves_water': False, 'inter_chain': False}
    if len(d) < 2:
        return summary, pd.DataFrame()

    i, j, sep = _offending_pairs(d, xyz, elem, alt, occ, clash_cutoff, min_occupancy)
    if len(i) == 0:
        # nothing clashes: still report how close the structure comes
        _, _, wide = _offending_pairs(d, xyz, elem, alt, occ, report_cutoff, min_occupancy)
        summary['min_contact'] = float(wide[0]) if len(wide) else np.nan
        return summary, pd.DataFrame()

    chain = d['chain'].astype(str).to_numpy()
    resid = d['resid'].astype(int).to_numpy()
    resname = d['resname'].astype(str).str.upper().str.strip().to_numpy()
    name = d['name'].astype(str).str.upper().str.strip().to_numpy()
    is_protein = np.isin(resname, STANDARD_RESIDUES)
    is_water = np.isin(resname, WATER_RESNAMES)
    modelled = {(str(c), int(r)) for c, r in modelled_residues}
    is_modelled = np.array([(c, r) in modelled for c, r in zip(chain, resid)]) | (occ == 0.0)

    out = pd.DataFrame({
        'sep': np.round(sep, 3),
        'chain_a': chain[i], 'resname_a': resname[i], 'resid_a': resid[i], 'name_a': name[i],
        'chain_b': chain[j], 'resname_b': resname[j], 'resid_b': resid[j], 'name_b': name[j],
        'inter_chain': chain[i] != chain[j],
        'modelled': is_modelled[i] | is_modelled[j],
        'hetatm': ~(is_protein[i] & is_protein[j]),
        'water': is_water[i] | is_water[j],
    })

    summary.update({
        'min_contact': float(sep[0]),
        'n_clashes': int(len(out)),
        'worst': (f"{out.chain_a[0]}/{out.resname_a[0]}{out.resid_a[0]}:{out.name_a[0]} <-> "
                  f"{out.chain_b[0]}/{out.resname_b[0]}{out.resid_b[0]}:{out.name_b[0]}"),
        'involves_modelled': bool(out['modelled'].any()),
        'involves_hetatm': bool(out['hetatm'].any()),
        'involves_water': bool(out['water'].any()),
        'inter_chain': bool(out['inter_chain'].any()),
    })
    return summary, out
