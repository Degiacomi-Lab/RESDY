'''
Biological assembly of a curated structure.

A crystal structure deposits the asymmetric unit, which may hold part of the biological
assembly or several copies of it. REMARK 350 of the PDB file lists one or more
BIOMOLECULE entries, each giving the BIOMT operators that build an assembly from chains
of the asymmetric unit. Crystal packing operators (REMARK 290) are not used.

Curation works on the asymmetric unit, and the assembly is built from the curated
structure, so that every copy of a chain carries the same rebuilt residues. The copies
are environment: every feature is computed on the whole assembly, but only residues of
the original chains are reported, so that a residue is measured once however many copies
of it the assembly holds. The outcome is recorded in ``<structure>.assembly.json`` beside
the structure, and :func:`measured_chains` reads it back.
'''

import os
import json
import numpy as np
import biobox as bb
from .helper import require_biobox
require_biobox(bb)
from .geometry import check_geometry


#: Suffix of the record written beside each curated structure.
RECORD_SUFFIX = '.assembly.json'


def read_biomolecules(path):
    '''
    Read the BIOMOLECULE entries of REMARK 350.

    :param path: pdb file as downloaded from the PDB.
    :type path: str
    :returns: (biomatrix, determined). biomatrix is ``{biomolecule id: [(chains,
        matrices), ...]}`` as biobox stores it, with chains None when the operators apply
        to every chain and matrices an array of 3x4 [R|t] operators. determined maps each
        id to 'author' when the entry carries an AUTHOR DETERMINED BIOLOGICAL UNIT line,
        'software' when it carries only a SOFTWARE DETERMINED QUATERNARY STRUCTURE line,
        and 'unspecified' otherwise. Both are empty when the file has no BIOMT operator.
    :rtype: tuple
    '''
    determined = {}
    current = None
    has_biomt = False
    with open(path) as fin:
        for line in fin:
            if not line.startswith('REMARK 350'):
                continue
            text = line[10:].strip()
            if text.startswith('BIOMT'):
                has_biomt = True
            elif text.startswith('BIOMOLECULE:'):
                current = int(text.split(':')[1])
                determined[current] = 'unspecified'
            elif current is not None and text.startswith('AUTHOR DETERMINED'):
                determined[current] = 'author'
            elif (current is not None and text.startswith('SOFTWARE DETERMINED')
                  and determined[current] != 'author'):
                determined[current] = 'software'

    if not has_biomt:
        return {}, {}

    M = bb.Molecule()
    M.import_pdb(path, include_hetatm=True)
    biomatrix = M.properties.get('biomatrix', {})
    determined = {k: determined.get(k, 'unspecified') for k in biomatrix}
    return biomatrix, determined


def choose_biomolecule(biomatrix, determined):
    '''
    Choose the assembly to build: the first author-determined biomolecule, else the first
    listed.

    :param biomatrix: as returned by :func:`read_biomolecules`.
    :type biomatrix: dict
    :param determined: as returned by :func:`read_biomolecules`.
    :type determined: dict
    :returns: (biomolecule id, how it was determined), or (None, '') when there is none.
    :rtype: tuple
    '''
    ids = sorted(biomatrix)
    if not ids:
        return None, ''
    for i in ids:
        if determined.get(i) == 'author':
            return i, 'author'
    return ids[0], determined.get(ids[0], 'unspecified')


def record_path(path):
    '''
    Path of the assembly record of a structure. An energy-minimised copy
    (``<structure>_relaxed.pdb``) shares the record of the structure it was made from.

    :param path: pdb file of the structure.
    :type path: str
    :returns: path of ``<structure>.assembly.json``.
    :rtype: str
    '''
    stem = os.path.splitext(os.path.basename(path))[0]
    if stem.endswith('_relaxed'):
        stem = stem[:-len('_relaxed')]
    return os.path.join(os.path.dirname(path), f'{stem}{RECORD_SUFFIX}')


def measured_chains(path):
    '''
    Chains of a structure whose residues are reported: the chains of the asymmetric unit
    that are part of the assembly. Copies made by the assembly operators are left out.

    :param path: pdb file of the structure.
    :type path: str
    :returns: set of chain names, or None when the structure has no assembly record, in
        which case every chain is reported.
    :rtype: set
    '''
    rec = record_path(path)
    if not os.path.exists(rec):
        return None
    with open(rec) as fin:
        return set(json.load(fin)['measured_chains'])


def _is_identity(matrices):
    return all(np.allclose(m, np.eye(3, 4), atol=1e-4) for m in matrices)


def build_assembly(path, biomatrix, determined, max_atoms=None):
    '''
    Replace a curated structure with its biological assembly, and write its record.

    The assembly built is the one :func:`choose_biomolecule` picks. Chains the biomolecule
    does not list are left out, as in the assembly files of the PDB. The first copy of a
    chain keeps its name and later copies take names not used in the structure. The
    structure is left as it is when the operators are all identities covering every
    chain, when no listed chain is present, when the assembly would need chain names
    longer than one character (PROPKA and Biopython read only column 22 of a pdb file), or
    when it would hold more than ``max_atoms`` atoms. Every chain is then reported.

    The record holds the biomolecule id, how it was determined, the status ('built',
    'identity', 'none' when the file has no BIOMT operator, or 'kept' when the asymmetric
    unit was kept for one of the reasons above, given in 'reason'), the number of
    operators, the chains of the asymmetric unit and of the written structure, the
    measured chains, the number of atoms written and, for a built assembly, the number of
    non-bonded clashes ``resdy.geometry.check_geometry`` finds in it.

    :param path: curated pdb file, overwritten with the assembly.
    :type path: str
    :param biomatrix: as returned by :func:`read_biomolecules`.
    :type biomatrix: dict
    :param determined: as returned by :func:`read_biomolecules`.
    :type determined: dict
    :param max_atoms: largest assembly to build, in atoms; None for no limit.
    :type max_atoms: int
    :returns: the record written.
    :rtype: dict
    '''
    M = bb.Molecule()
    M.import_pdb(path, include_hetatm=True)
    original = list(dict.fromkeys(M.data['chain'].values))

    biomolecule, how = choose_biomolecule(biomatrix, determined)
    record = {'biomolecule': biomolecule, 'determined': how, 'status': '', 'reason': '',
              'n_operators': 0, 'original_chains': original, 'chains': original,
              'measured_chains': original, 'n_atoms': int(len(M.data)), 'n_clashes': None}

    if biomolecule is None:
        record['status'] = 'none'
        return _write_record(path, record)

    groups = biomatrix[biomolecule]
    record['n_operators'] = int(sum(len(mats) for _, mats in groups))
    if any(c is None for c, _ in groups):
        covered = set(original)
    else:
        covered = set().union(*(c for c, _ in groups)) & set(original)
    if covered == set(original) and all(_is_identity(mats) for _, mats in groups):
        record['status'] = 'identity'
        return _write_record(path, record)

    M.properties['biomatrix'] = {biomolecule: groups}
    try:
        A = M.apply_biomatrix(biomolecule=biomolecule)
    except ValueError as e:
        record['status'] = 'kept'
        record['reason'] = str(e)
        return _write_record(path, record)

    chains = list(dict.fromkeys(A.data['chain'].values))
    if any(len(c) > 1 for c in chains):
        record['status'] = 'kept'
        record['reason'] = f'the assembly has {len(chains)} chains, more than single-character chain names allow'
        return _write_record(path, record)
    if max_atoms is not None and len(A.data) > max_atoms:
        record['status'] = 'kept'
        record['reason'] = f'the assembly has {len(A.data)} atoms, more than max_atoms={max_atoms}'
        return _write_record(path, record)

    A.write_pdb(path)
    record.update({'status': 'built', 'chains': chains,
                   'measured_chains': [c for c in original if c in chains],
                   'n_atoms': int(len(A.data))})
    _write_record(path, record)
    record['n_clashes'] = int(check_geometry(path)[0]['n_clashes'])
    return _write_record(path, record)


def _write_record(path, record):
    with open(record_path(path), 'w') as fout:
        json.dump(record, fout, indent=2)
    return record
