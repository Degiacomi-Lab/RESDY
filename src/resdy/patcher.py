'''
Provide patching functions which are called from protein.py
These interact with Modeller to perform patching on protein structures.
'''
import fileinput
import glob
import os
import sys
import re
from textwrap import wrap
import shutil
from copy import deepcopy
import numpy as np
import biobox as bb

try:
    from modeller import *
    from modeller.automodel import *
    modeller_available = True
    _modeller_error = ''
except Exception as e:
    modeller_available = False
    _modeller_error = e

from .helper import ShutUp, require_biobox
from .geometry import check_geometry

require_biobox(bb)


def autopatch(tmp_folder, fbasename, gap_cutoff=8):
    '''
    Perform modelling on the given fbasename. Goes through the steps within individual functions:
    _pdb_to_seq(), _fasta_to_pir(), _full_align(), _trim_align(), _gap_check(), _patch_model(). Then
    removes any unnecessary files leftover from the patching.

    :param tmp_folder: Relative location path to tmp_folder, used to write the new alignment files
        to
    :type tmp_folder: str
    :param fbasename: basename to use for the file being used. In the overall code this includes the
        path to the temporary folder
    :type fbasename: str
    :param gap_cutoff: Number of amino acids to allow for as gaps in sequence for patching. Usually
        taken through from what was set in protein.py
    :type gap_cutoff: int
    :returns: patched chain name path
    :rtype: str
    '''
    if not modeller_available:
        raise ImportError(f'>> To patch a structure which has gaps in requires Modeller which is '
                          f'not currently installed. Install modeller with the instructions in '
                          f'the RESDY Github Repo for more details. Error: {str(_modeller_error)}')

    print('>> modelling missing residues')
    pdb_out = ''
    seq_name = ''
    built = []

    # Modeller is given fbasename as a PDB *code*, which it resolves against
    # env.io.atom_files_directory, and as the align_codes written into the PIR header.
    # Neither tolerates an absolute path: the lookup cannot resolve one, and on Windows
    # the drive letter adds a colon to a header whose fields are colon-separated. Work
    # from the parent of the temporary folder so that every name handed to Modeller is
    # relative and short.
    abs_tmp = os.path.abspath(tmp_folder)
    work_dir = os.path.dirname(abs_tmp)
    rel_tmp = os.path.basename(abs_tmp)
    rel_base = os.path.join(rel_tmp, os.path.basename(fbasename))
    cwd = os.getcwd()

    try:
        os.chdir(work_dir)

        _pdb_to_seq(rel_base)
        _fasta_to_pir(rel_base)
        seq_name = _full_align(rel_tmp, rel_base)
        _trim_align(rel_tmp, "alignment.seg.ali")
        patch_status = _gap_check(rel_tmp, "trimmed_align.ali", gap_cutoff)

        if patch_status:
            pdb_out, built = _patch_model(rel_tmp, rel_base, seq_name)
            pdb_out = os.path.join(abs_tmp, os.path.basename(pdb_out))
        else:
            pdb_out = ""

    except Exception as e:
        raise Exception(f">> ERROR on autopatching the structure: {e}") from e

    finally:
        # the intermediates were written relative to work_dir, which is still the
        # working directory here
        autopatch_files = [f'{rel_base}.seq', f'{rel_base}.pir',
                           os.path.join(rel_tmp, 'alignment.seg'),
                           os.path.join(rel_tmp, 'alignment.seg.ali'),
                           os.path.join(rel_tmp, 'trimmed_align.ali'), 'family.mat']
        for pattern in ('.ini', '.rsr', '.sch', '.V*', '.D*'):
            autopatch_files.extend(glob.glob(f'{rel_tmp}{pattern}'))

        for patch_file in autopatch_files:
            if not os.path.exists(patch_file):
                continue
            try:
                os.remove(patch_file)
            except Exception as e:
                print(f"cannot remove {patch_file}, error: {e}; continuing...")

        os.chdir(cwd)

    return pdb_out, built

#autopatch step 1a. pir format of AA from pdb
def _pdb_to_seq(fbasename):
    '''
    Autopatch step 1a: Use the modeller package to extract the sequence from the pdb file given and
    save this to a file with extension .seq. This file will be saved in the temporary directory
    created for the curation for this protein.

    :param fbasename: basename to use for the file being used. In the overall code this includes the
        path to the temporary folder
    :type fbasename: str
    '''
    env = Environ()
    env.io.two_char_chain = True
    mdl = Model(env, file=fbasename)
    aln = Alignment(env)
    aln.append_model(mdl, align_codes=fbasename)
    aln.write(file=f'{fbasename}.seq')

#autopatch step 1b. pir from complete AA fasta
def _fasta_to_pir(fbasename):
    '''
    Autopatch step 1b: Use the modeller package to transform the fasta file into a .pir file located
    in the temporary directory for the specific protein.

    :param fbasename: basename to use for the file being used. In the overall code this includes the
        path to the temporary folder
    :type fbasename: str
    '''
    env = Environ()
    env.io.two_char_chain = True
    a = Alignment(env, file=fbasename+".fasta", alignment_format='FASTA')
    a.write(file=f'{fbasename}.pir', alignment_format='PIR')

#autopatch step 2. add sequence name to 2nd line; copy the pir contents and structure info into alignment.seg; align sequences and generate model
def _full_align(tmp_folder, fbasename):
    '''
    Autopatch step 2: Take the curated .pir and .seq files and place them into an alignment.seg file
    located in the temporary directory

    :param tmp_folder: Relative location path to tmp_folder, used to write the new alignment files
        to
    :type tmp_folder: str
    :param fbasename: basename to use for the file being used. In the overall code this includes the
        path to the temporary folder
    :type fbasename: str
    '''
    pir_fname = f'{fbasename}.pir'
    seq_fname = f'{fbasename}.seq'

    f = open(pir_fname, "r")
    f1 = f.readlines()
    f.close()

    for i, line in enumerate(f1):
        if "P1;" in line:
            seq_name = line[4:8]
        if "sequence:" in line:
            seq_pos = i

    P1 = ">P1;"+seq_name+"\n"
    seq_line = "sequence:"+seq_name+":::::::-1.00:-1.00\n"
    AA_block = ''.join([str(elem) for elem in f1[seq_pos+1:]])
    seq_block = P1 + seq_line + AA_block
    f = open(pir_fname, 'w')
    f.writelines(seq_block)
    f.close()

    alignment_filename = f'{tmp_folder}{os.sep}alignment.seg'

    if sys.platform == "win32":
        my_cmd_a = f'type {pir_fname} {seq_fname} > {alignment_filename}'
    else:
        my_cmd_a = f'cat {pir_fname} {seq_fname} > {alignment_filename}'

    os.system(my_cmd_a)

    env = Environ()
    env.io.two_char_chain = True
    env.io.atom_files_directory = [tmp_folder, '.', f'..{os.sep}atom_files']
    a = AutoModel(env,
                  # file with template codes and target sequence
                  alnfile  = alignment_filename,
                  # PDB codes of the templates
                  knowns   = fbasename,
                  # code of the target
                  sequence = seq_name)
    a.auto_align() # get an automatic alignment (alignment.seg.ali)
    return seq_name

# autopatch step 3. trim the alignment by removing gaps for missing residues at the termini of the structure
def _trim_align(tmp_folder, align_file):
    '''
    Autopatch step 3: Take the alignment file, find the different sequence blocks, remove empty
    lines, split into seq and structure parts, reformat and then print structure parts and then seq
    parts after.

    :param tmp_folder: Relative location path to tmp_folder, used to write the new alignment files
        to
    :type tmp_folder: str
    :param fbasename: basename to use for the file being used. In the overall code this includes the
        path to the temporary folder
    :type fbasename: str
    '''
    align_file = f'{tmp_folder}{os.sep}alignment.seg.ali'
    f=open(align_file, "r")
    f1 = f.readlines()
    P1_pos = []

    #identify positions of different sequences (P1) blocks
    for i, line in enumerate(f1):
        if "P1;" in line:
            P1_pos.append(i)
    f.close()
    sec_1 = P1_pos[0]
    sec_2 = P1_pos[1]

    #remove empty lines and split into 2 blocks: seq and structure
    AA_struc = ''.join([str(elem.rstrip("\n").rstrip("*")) for elem in f1[sec_1+2:sec_2]])

    start_gaps = re.findall('^[-]+', AA_struc)
    end_gaps = re.findall('[-]+$', AA_struc)

    if len(start_gaps)>0:
        start_gaps_len = len(start_gaps[0])
    else:
        start_gaps_len = 0

    if len(end_gaps)>0:
        end_gaps_len = len(end_gaps[0])
    else:
        end_gaps_len = 0

    AA_struc_len = len(AA_struc)
    AA_struc_new = AA_struc[start_gaps_len : AA_struc_len - end_gaps_len]+"*"
    AA_struc_new = wrap(AA_struc_new, 75) #split after 75 characters

    AA_seq = ''.join([str(elem.rstrip("\n").rstrip("*")) for elem in f1[sec_2+2:]])
    AA_seq_len = len(AA_seq)
    AA_seq_new = AA_seq[start_gaps_len : AA_seq_len - end_gaps_len]+"*"
    AA_seq_new = wrap(AA_seq_new, 75) 

    AA_struc_new = '\n'.join([str(elem) for elem in AA_struc_new])
    AA_seq_new = '\n'.join([str(elem) for elem in AA_seq_new])

    struc_sec = f1[sec_1] + f1[sec_1+1] + AA_struc_new + "\n"
    seq_sec = f1[sec_2] + f1[sec_2+1] + AA_seq_new + "\n"
    trim_filename = f'{tmp_folder}{os.sep}trimmed_align.ali'
    f = open(trim_filename, 'w')
    f.writelines(struc_sec)
    f.writelines(seq_sec)
    f.close()

#autopatch step 4. Check if any gap is more than cutoff length in the trimmed_align.ali and if so set patch_status = "no"
def _gap_check(tmp_folder, align_file, gap_cutoff):
    '''
    Autopatch step 4: Take trimmed align file, construct amino acid sequence and find the gaps in
    the sequence, if any gap length bigger than cutoff return False, else return True for patching.

    :param tmp_folder: Relative location path to tmp_folder, used to write the new alignment files
        to
    :type tmp_folder: str
    :param align_file: path to the alignment file
    :type align_file: str
    :param gap_cutoff: Number of amino acids to allow for as gaps in sequence for patching. Usually
        taken through from what was set in protein.py
    :type gap_cutoff: int
    :returns: Returns true (default) if there are no gaps or gaps shorter than the gap cutoff, False
        if there is a gap bigger than cutoff which will not be patched
    :rtype: bool
    '''
    patch_status = True
    align_file = f'{tmp_folder}{os.sep}trimmed_align.ali'
    f=open(align_file, "r")
    f1 = f.readlines()
    P1_pos = []
    for i, line in enumerate(f1):
        if "P1;" in line: #identify positions of different sequences (P1) blocks
            P1_pos.append(i)
    f.close()
    sec_1 = P1_pos[0]
    sec_2 = P1_pos[1]

    AA_struc = ''.join([str(elem.rstrip("\n").rstrip("*")) for elem in f1[sec_1    +2:sec_2]])
    #check for gap_lengths in AA_struc
    struc_gaps = re.findall('[-]+', AA_struc)
    for i, struc_gap in enumerate(struc_gaps):
        gap_len = len(struc_gap)
        if gap_len > gap_cutoff:
            print(f">> structure not patched. Long sequence gap: {str(gap_len)}")
            patch_status = False

    return patch_status

#autopatch step 5. build missing residues
def _patch_model(tmp_folder, fbasename, seq_name):
    '''
    Autopatch step 5: Perform the patching on the structure, sets root name for temporary files
    created in normal directory (as all will be using same place on parallel), moves curated file
    for the chain into the temporary folder.

    :param tmp_folder: Relative location path to tmp_folder, used to write the new alignment files
        to
    :type tmp_folder: str
    :param fbasename: basename to use for the file being used. In the overall code this includes the
        path to the temporary folder
    :type fbasename: str
    :param seq_name: PDB name for target
    :type seq_name: str
    :returns: patched chain name path
    :rtype: str
    '''
    align_file = f'{tmp_folder}{os.sep}trimmed_align.ali'
    gaps = alignment_gaps(align_file)

    if not gaps:
        # the template covers the whole target, so there is nothing to build. Running
        # AutoModel anyway would rebuild and optimise coordinates that were measured
        # experimentally, moving every atom by around an Angstrom for no gain.
        print('>> nothing missing from this chain, leaving its coordinates untouched')
        return f'{fbasename}.pdb', []

    print(f">> patching model, building {sum(e - s + 1 for s, e in gaps)} residue(s) "
          f"in {len(gaps)} gap(s); the rest is held fixed")
    log.verbose()
    env = Environ()
    env.io.two_char_chain = True
    env.io.atom_files_directory = [tmp_folder, '.', f'..{os.sep}atom_files']

    class _GapModel(AutoModel):
        '''Optimise only the residues the template does not provide.'''

        def select_atoms(self):
            # select by position rather than by a 'number:chain' specification: the
            # model's own numbering and chain naming are Modeller's to choose, while
            # the positions come straight from the alignment
            sel = Selection()
            for first, last in gaps:
                for i in range(first - 1, last):
                    sel.add(self.residues[i])
            return sel

    a = _GapModel(env,
                  # file with template codes and target sequence
                  alnfile  = f'{tmp_folder}{os.sep}trimmed_align.ali',
                  # PDB codes of the templates
                  knowns   = fbasename,
                  # code of the target
                  sequence = seq_name,
                  assess_methods = (assess.DOPE, assess.GA341),
                  root_name=os.path.basename(tmp_folder))
    a.md_level = refine.fast #very_fast, fast, slow, very_slow, slow_large, refine
    #repeat whole cycle twice and do not stop unless obj. func > 1e6
    #a.repeat_optimization = 2
    a.max_molpdf = 1e6
    a.make()

    pdb_out = f'{fbasename}_PATCHED.pdb'

    os.rename(f'{os.path.basename(tmp_folder)}.B99990001.pdb', pdb_out)

    # Put the model back in the frame of the chain it was built from. This is done here
    # rather than in curate() because the correspondence between the two comes from the
    # alignment, which is only available at this point.
    pairs = alignment_pairs(align_file)
    n_fit, rmsd = superpose_onto(pdb_out, f'{fbasename}.pdb', pairs=pairs)
    if n_fit:
        print(f'>> model put back on the experimental coordinates over {n_fit} atoms, '
              f'rmsd {rmsd:.2f} A')

    # Modeller's own occupancy and B-factor columns hold neither, so take them from the
    # template and mark the atoms that were built
    n_template, n_built = restore_template_columns(pdb_out, f'{fbasename}.pdb', pairs=pairs)
    print(f'>> occupancy and B-factor taken from the template for {n_template} atoms, '
          f'{n_built} built atoms marked with occupancy 0')

    # the residues of the model that the template did not provide: these are the only
    # coordinates that are a prediction rather than a measurement
    built = [i for first, last in gaps for i in range(first - 1, last)]
    return pdb_out, built


def alignment_gaps(align_file):
    '''
    Residues of the model that the template does not provide coordinates for.

    These are the only residues that need building. Everything else was determined
    experimentally and should be left exactly as deposited.

    :param align_file: trimmed PIR alignment, template block first.
    :type align_file: str
    :returns: list of (first, last) residue numbers in the model's own numbering, which
        runs from 1 over the target sequence. Empty when nothing is missing.
    :rtype: list
    '''
    seqs, cur = [], []
    with open(align_file) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith('>P1;'):
                if cur:
                    seqs.append(''.join(cur))
                    cur = []
            elif line.startswith(('structureX:', 'structure:', 'sequence:')):
                continue
            elif line:
                cur.append(line.replace('*', ''))
    if cur:
        seqs.append(''.join(cur))
    if len(seqs) < 2:
        return []

    template, target = seqs[0], seqs[1]
    ranges, idx, start = [], 0, None
    for t_char, q_char in zip(template, target):
        if q_char == '-':
            continue                       # not part of the model at all
        idx += 1                           # model residues are numbered from 1
        if t_char == '-':
            if start is None:
                start = idx
        elif start is not None:
            ranges.append((start, idx - 1))
            start = None
    if start is not None:
        ranges.append((start, idx))
    return ranges


def merge_heteroatoms(source, target):
    '''
    Carry the heteroatom records of ``source`` into ``target``.

    Modeller returns only the polymer it modelled, and the complex is reassembled from
    those models, so anything the user asked to keep is lost on any chain that was
    patched. The coordinates are taken from ``source`` unchanged, which is correct only
    because a patched chain is superposed back onto the frame it was built from; without
    that, a ligand would end up several Angstrom from its site.

    :param source: file holding the heteroatoms, normally the structure curate() was given.
    :type source: str
    :param target: reassembled file to write them into. Any heteroatom already present is
        replaced, so calling this twice does not duplicate them.
    :type target: str
    :returns: number of heteroatom records carried over.
    :rtype: int
    '''
    het = [l for l in open(source) if l.startswith('HETATM')]
    if not het:
        return 0

    kept = [l for l in open(target) if not l.startswith('HETATM')]
    out, inserted = [], False
    for line in kept:
        if not inserted and line.startswith('END') and not line.startswith('ENDMDL'):
            out.extend(het)
            inserted = True
        out.append(line)
    if not inserted:
        out.extend(het)

    with open(target, 'w') as fh:
        fh.writelines(out)
    return len(het)


def analyze_protein(M):
    '''
    Take a biobox molecule instance, look for gaps in the sequence and return 3 elements list:
    [number of gaps, number of missing residues, largest sequence gap]

    Only polymer residues are considered, identified by carrying a CA atom. Waters,
    ions and ligands are numbered in their own range, often continuing past the end of
    the chain they sit in, so counting them would put the whole span between the last
    residue and the first heteroatom down as missing. On 3DBJ, whose chains run to
    residue 174 with waters numbered from 202, including them turns two real gaps of
    two residues into a single phantom gap of 28 and the structure is rejected.

    :param M: bb.Molecule instance for the chain of interest to get gaps over
    :type M: Biobox molecule instance
    :returns: [number of gaps, number of missing residues, largest sequence gap]. All
        zero if the chain holds no polymer residue at all.
    :rtype: list
    '''
    # biobox reads columns 12 to 17, so an alternate-location indicator arrives stuck to
    # the atom name ('CA A'); compare on the first whitespace-separated token
    names = M.data['name'].astype(str).str.split().str[0]
    res = np.unique(M.data.loc[names == 'CA', 'resid'].values)

    missing = []
    patch = []
    cnt = [0, 0, 0]
    if len(res) == 0:
        return cnt

    for r in range(int(np.min(res)), int(np.max(res)+1)):
        if r in res:
            if len(patch) > 0:
                missing.append(deepcopy(patch))
                cnt[0] += 1
                cnt[1] += len(patch)
                if len(patch) > cnt[2]:
                    cnt[2] = len(patch)

                patch = []

        else:
            patch.append(r)

    return cnt


def alignment_pairs(align_file):
    '''
    Which residue of the model corresponds to which residue of the template.

    Taken from the alignment itself rather than inferred. Neither residue numbers nor
    ordinal positions can be trusted: Modeller renumbers its model from 1, and once it
    inserts a residue every ordinal after the insertion shifts, so both proxies silently
    pair atoms of different residues.

    :param align_file: trimmed PIR alignment, template block first.
    :type align_file: str
    :returns: list of (model residue index, template residue index), both 0-based and
        counted in file order over the residues each actually contains.
    :rtype: list
    '''
    seqs, cur = [], []
    with open(align_file) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith('>P1;'):
                if cur:
                    seqs.append(''.join(cur))
                    cur = []
            elif line.startswith(('structureX:', 'structure:', 'sequence:')):
                continue
            elif line:
                cur.append(line.replace('*', ''))
    if cur:
        seqs.append(''.join(cur))
    if len(seqs) < 2:
        return []

    template, target = seqs[0], seqs[1]
    pairs, t_i, q_i = [], -1, -1
    for t_char, q_char in zip(template, target):
        if t_char != '-':
            t_i += 1
        if q_char != '-':
            q_i += 1
        if t_char != '-' and q_char != '-':
            pairs.append((q_i, t_i))
    return pairs


def superpose_onto(mobile, reference, out='', pairs=None):
    '''
    Rigidly move a patched chain back onto the coordinates it was built from.

    Modeller returns a model in its own frame, so a chain comes back from
    :func:`autopatch` translated and rotated with respect to the structure it was taken
    from. For a single chain that is harmless, but the chains of a complex are patched
    independently, so reassembling them without putting each one back destroys the
    quaternary structure. This performs a least-squares fit (Kabsch) on the atoms the two
    files share, keyed on residue number and atom name, and applies the resulting
    transformation to every atom of ``mobile``.

    Coordinates are rewritten in the text of the file rather than through biobox, so that
    the occupancy and B-factor columns are left exactly as they were.

    :param mobile: pdb file to move, normally the output of :func:`autopatch`.
    :type mobile: str
    :param reference: pdb file holding the original coordinates of the same chain.
    :type reference: str
    :param out: file to write. Defaults to overwriting ``mobile``.
    :type out: str
    :param pairs: (model residue index, reference residue index) correspondences, as
        :func:`alignment_pairs` returns them. Without it the two files are assumed to
        hold the same residues in the same order, which is only true when nothing was
        built.
    :type pairs: list
    :returns: (number of atoms fitted, RMSD over those atoms after the fit, in A). Both
        are zero when the fit could not be made, in which case ``mobile`` is left alone.
    :rtype: tuple
    '''
    out = out or mobile

    def _atoms(path):
        rows = []
        for line in open(path):
            if line.startswith(('ATOM', 'HETATM')):
                rows.append((line, line[21], int(line[22:26]), line[12:16].strip(),
                             (float(line[30:38]), float(line[38:46]), float(line[46:54]))))
        return rows

    mob, ref = _atoms(mobile), _atoms(reference)

    def _by_residue(rows):
        '''{residue index in file order: {atom name: coordinates}}'''
        out, ordinals = {}, {}
        for line, chain, resid, name, xyz in rows:
            res = (chain, resid)
            if res not in ordinals:
                ordinals[res] = len(ordinals)
            out.setdefault(ordinals[res], {})[name] = xyz
        return out

    mob_by, ref_by = _by_residue(mob), _by_residue(ref)
    if pairs is None:
        pairs = [(i, i) for i in sorted(set(mob_by) & set(ref_by))]
    P, Q = [], []
    for m_idx, r_idx in pairs:
        ref_res = ref_by.get(r_idx, {})
        for name, xyz in mob_by.get(m_idx, {}).items():
            if name in ref_res:
                P.append(xyz)
                Q.append(ref_res[name])
    if len(P) < 3:
        print(f'>> Only {len(P)} atoms shared between {os.path.basename(mobile)} and the '
              f'chain it was built from; cannot superpose it back, leaving it as Modeller '
              f'placed it')
        return 0, 0.0

    P, Q = np.array(P, dtype=float), np.array(Q, dtype=float)
    Pc, Qc = P.mean(axis=0), Q.mean(axis=0)
    U, _, Vt = np.linalg.svd((P - Pc).T @ (Q - Qc))
    # guard against the fit producing a reflection rather than a rotation
    D = np.diag([1.0, 1.0, np.sign(np.linalg.det(Vt.T @ U.T))])
    rot = Vt.T @ D @ U.T
    trans = Qc - rot @ Pc
    rmsd = float(np.sqrt((((rot @ P.T).T + trans - Q) ** 2).sum(axis=1).mean()))

    with open(out, 'w') as fh:
        for line, _, _, _, xyz in mob:
            v = rot @ np.array(xyz, dtype=float) + trans
            fh.write(f'{line[:30]}{v[0]:8.3f}{v[1]:8.3f}{v[2]:8.3f}{line[54:]}')
    return len(P), rmsd


#: Occupancy and B-factor columns (pdb columns 55-66) given to an atom that Modeller built.
#: An occupancy of zero marks it in the file itself, so that any later reader can tell a
#: predicted position from a measured one without a side file.
BUILT_ATOM_COLUMNS = f'{0.0:6.2f}{0.0:6.2f}'


def restore_template_columns(model, template, out='', pairs=None):
    '''
    Give every atom of a Modeller model the occupancy and B-factor of the template atom it
    was copied from, and mark the atoms Modeller built.

    Modeller does not keep these columns: its models carry the template B-factor in the
    occupancy column and a per-residue value of its own in the B-factor column. An atom
    the template provides takes the template's occupancy and B-factor, copied as text. An
    atom the template does not provide, either because its residue was missing or because
    the template residue lacked it, takes :data:`BUILT_ATOM_COLUMNS`.

    :param model: pdb file written by Modeller.
    :type model: str
    :param template: pdb file of the chain the model was built from.
    :type template: str
    :param out: file to write. Defaults to overwriting ``model``.
    :type out: str
    :param pairs: (model residue index, template residue index) correspondences, as
        :func:`alignment_pairs` returns them. Without it the two files are assumed to hold
        the same residues in the same order.
    :type pairs: list
    :returns: (number of atoms given template values, number of atoms marked as built).
    :rtype: tuple
    '''
    out = out or model

    def _residues(path):
        '''list of atom lines, and {residue index in file order: {atom name: columns 55-66}}'''
        lines, by_res, ordinals = [], {}, {}
        for line in open(path):
            if line.startswith(('ATOM', 'HETATM')):
                res = (line[21], line[22:27])
                if res not in ordinals:
                    ordinals[res] = len(ordinals)
                by_res.setdefault(ordinals[res], {})[line[12:16].strip()] = line[54:66]
                lines.append((line, ordinals[res]))
            else:
                lines.append((line, None))
        return lines, by_res

    model_lines, model_by = _residues(model)
    _, template_by = _residues(template)
    if pairs is None:
        pairs = [(i, i) for i in sorted(set(model_by) & set(template_by))]
    template_of = dict(pairs)

    n_template, n_built = 0, 0
    with open(out, 'w') as fh:
        for line, res_idx in model_lines:
            if res_idx is None:
                fh.write(line)
                continue
            columns = template_by.get(template_of.get(res_idx), {}).get(line[12:16].strip())
            if columns is None:
                columns = BUILT_ATOM_COLUMNS
                n_built += 1
            else:
                n_template += 1
            fh.write(f'{line[:54]}{columns}{line[66:]}')
    return n_template, n_built


def fragment(pdb, fasta, outfolder=".", include_hetatm=False):
    '''
    Take pdb file, check if more than 62 chains (can't patch this due to legacy pdb format problems)
    if so write to the given file, check if double letter chain names are present (as while we could
    patch, some measures features can't read them and would cause problems further on) if so write
    these pdb codes to a file. Split pdb up into the chains, run a gap calculation on each and add
    to the gap report list which it returns. Split up fasta file into chains with corrections for
    modified AAs and check pdb and fasta agree for chains. Some code remains for dealing with double
    letter chain names but this should never be reached with earlier block in place.

    :param pdb: Name of the pdb file of interest to fragment
    :type pdb: str
    :param fasta: Name of the fasta file corresponding to the pdb file of interest
    :type fasta: str
    :param outfolder: path of the temporary folder for the patching
    :type outfolder: str
    :param include_hetatm: Toggleable option to allow hetatms to pass through biobox
    :type include_hetatm: bool
    :returns: (gap report, chain names). The gap report has one row per chain, in the same
        order as the chain names, holding [number of gaps, number of missing residues,
        largest sequence gap].
    :rtype: tuple
    '''

    os.makedirs(outfolder, exist_ok=True)

    #split PDB file in chains using biobox
    M = bb.Molecule()
    M.import_pdb(pdb, include_hetatm=include_hetatm)
    chains = np.unique(M.data["chain"].values)

    # Catch cases which have more than 62 chains which are not able to be worked with using PDB files
    if len(chains) > 62:
        with open(os.path.join(os.path.split(outfolder)[0], 'excessive_chains_files.txt'), 'a') as f1:
            f1.write(f'PDB: {pdb}; Num Chains: {len(chains)}; Chains: {str(chains)}')
        raise Exception(f'>> More than 62 individual chains present in the structure, not possible to curate a file for this number of chain names. File: {pdb}')

    # Catch cases which have double letter chain names, used to work out if worth including more functionality for these, have left future double letter functionality in incase this is useful in future.
    if any(len(chain) > 1 for chain in chains):
        with open(os.path.join(os.path.split(outfolder)[0], 'double_letter_file_count.txt'), 'a') as f1:
            f1.write(f'PDB: {pdb}; Num Chains: {len(chains)}; Chains: {", ".join(chains)}\n')
        raise Exception(f'>> Double letter chain names are present within the structure, currently ignoring these, tally updated in \'double_letter_file_count.txt\'. File: {pdb}')

    gap_count = []
    for c in chains:
        idxs = M.atomselect([c], '*', '*', get_index=True)[1]
        M_two = M.get_subset(idxs)
        M_two.write_pdb(os.path.join(outfolder, f"chain{c}.pdb"))
        gap_count.append(analyze_protein(M_two))

    #split FASTA
    fin = open(fasta, 'r')
    fasta_headers = []
    fasta_chains = []
    fasta_chain_auth_mapping = {}
    sequences = []
    for line in fin:
        if ">" in line:
            fasta_headers.append(line)
            chain_raw_info = line.split("|")[1][6:].split(",")
            if len(chain_raw_info[0]) == 1:
                chain_info = chain_raw_info
            elif '[' in chain_raw_info[0]:
                chain_info = [str(a.split('[')[0].strip()) for a in chain_raw_info]
                for a in chain_raw_info:
                    fasta_chain_auth_mapping[str(a.split('[')[0].strip())] = str(a.split('[')[1].split(']')[0][4:].strip())
            else:
                chain_info = [a.strip()[0] for a in chain_raw_info]
            fasta_chains.append(chain_info)
            if "sequence" in locals():
                sequences.append(sequence)

            sequence = []

        else:
            #replace non-canonical aminoacids in FASTA sequence
            if "KCX" in line:
                line = line.replace('(KCX)', 'K')
            if "MSE" in line :
                line = line.replace('(MSE)', 'M')

            sequence.append(line)

    sequences.append(sequence)
    fin.close()

    #write FASTA files
    for i, header in enumerate(fasta_headers):
        for c in fasta_chains[i]:
            if c not in chains:
                raise Exception(f"chain mismatch between PDB and FASTA. {fasta_chains}, {chains}")
            fout = open(os.path.join(outfolder, f"chain{c}.fasta"), "w")
            fout.write(header)
            fout.writelines(sequences[i])
            fout.close()

    # This next section is a replacement for the conversion to iodata.two_char_chain Modeller format
    # Instead convert the chain names back to single character chain names so Modeller can read this properly and doesn't convert all double letter chain names to A
    # A new function is added at the end of patching to convert the single chain names back to the corresponding double chain names - call protein.py function which does this already?
    # find the double chain name files and store in list to iterate through when converting the single character names
    all_files = os.listdir(outfolder)
    doubleletter_pdb_files = []
    for file in all_files:
        temp_parts = file.split('.')
        if file[0:5] == 'chain' and temp_parts[1] == 'pdb' and len(temp_parts[0][5:]) == 2:
            doubleletter_pdb_files.append(file)

    # iterate through the files and for each one replace the chain name with a single letter lowercase version of the chain name: eg. CA would go to c
    used_sl_chain_names = []
    sl_to_dl_chain_map = {}
    for file in doubleletter_pdb_files:
        db_chain_name = file.split('.')[0][5:]

        candidate_names = [db_chain_name.lower()[0], db_chain_name.lower()[1]] + [str(a) for a in range(10)]
        sl_chain_name = next((c for c in candidate_names if c not in used_sl_chain_names), None)

        if sl_chain_name is None:
            raise Exception(f'>> No free single letter chain names available for renaming double '
                            f'letter chain names; chain: {db_chain_name}')

        used_sl_chain_names.append(sl_chain_name)
        sl_to_dl_chain_map[sl_chain_name] = db_chain_name

        try:
            temp_file_path = os.path.join(outfolder, file)
            with fileinput.FileInput(temp_file_path, inplace = True) as f:
                for line in f:
                    try:
                        #On lines with 'ATOM', 'TER' or 'HETATM' the current chain name is replaced with the
                        #new single letter chain name and adjusted to match the correct pdb file format
                        if (line[:4] == 'ATOM') or (line[:3] == 'TER') or (line[:6] == 'HETATM'):
                            if len(line) > 22:
                                line = line[:21] + sl_chain_name + ' ' + line[23:]  # space needed here as double letter previously used and so without it you will get an artefact leaving a double letter chain name anyway eg AA would go to aA otherwise
                                print(line, end ='')
                            else:
                                print(line, end='')
                        else:
                            print(line, end='')
                    except:
                        print(line, end='')

        #If the protein fails, print error message with the error
        except Exception as e:
            raise Exception(f'Failed replacing chains for file {file}. Could not convert double letter chain names while fragmenting. {e}')

    return np.array(gap_count), list(chains)


def reassemble(pdbs, labels, outname, outdir, include_hetatm=True):
    '''
    For each chain pdb file in the folder, check if the chain name in the filename is different to
    the one in the pdb data for double letter ones and if needs reverting from the temporary patches
    put in for patching. Take all the single chain pdb files and write these back to an overall pdb
    file for the curated structure.

    :param pdbs: List of the chain pdb files in the temporary folder
    :type pdbs: list
    :param labels: List of the chains labels corresponding to the list of pdbs
    :type labels: list
    :param outname: Path of the curated structure overall to write to
    :type outname: str
    :param outdir: Path of the output directory
    :type outdir: str
    :param include_hetatm: Toggleable option to allow hetatms to pass through biobox
    :type include_hetatm: bool
    '''
    for pdb_file in pdbs:
        # take pdb_file and extract the chain name - need to find the exact format it needs to go back into that the code is expecting to reconvert it from
        dbletter = False
        thousand_chain = False
        chainname = os.path.splitext(os.path.basename(pdb_file))[0].split('_')[0][5:]
        if len(chainname) > 1:
            dbletter = True

        M_tmp = bb.Molecule()
        M_tmp.import_pdb(pdb_file, include_hetatm=include_hetatm)
        max_res = int(M_tmp.data['resid'].max())
        if max_res > 999:
            thousand_chain = True

        try:
            with fileinput.FileInput(pdb_file, inplace = True) as f:
                for line in f:
                    try:
                        #On lines with 'ATOM', 'TER' or 'HETATM'
                        #the current chain name is replaced with the new single letter chain name and adjusted to match the correct pdb file format
                        if ((line[:4] == 'ATOM') or (line[:3] == 'TER') or (line[:6] == 'HETATM')) and len(line) > 26:
                            if dbletter:
                                raise Exception(f'>> Double letter chain names not currently supported in '
                                                f'carbamylation codebase. Protein will not be curated.')
                            else:
                                # thousand_chain doesn't actually matter here as double letter chain names aren't being curated, treat as normal pdb
                                resid = line[22:26].strip()
                                line = line[:21] + chainname + '%4s' % resid + line[26:]

                        print(line, end ='')

                    except:
                        print(line, end='')

        except Exception as e:
            raise Exception(f'Failed replacing chains for file {pdb_file}. Could not reassemble chains. Error: {e}')

    # Take all the files in pdbs and turn each into a biobox Molecule object and append this to a list of monomers
    # create a biobox multimer from the list of monomors and write a new pdb file combining all these
    # Slower solution but works
    with open(outname, 'wb') as f_outname:
        for i, (f, f_chain) in enumerate(zip(pdbs, labels)):
            T = bb.Molecule()
            T.import_pdb(f, include_hetatm=include_hetatm)
            T.data['chain'] = f_chain
            T.write_pdb(f)
            with open(f, 'rb') as f_new:
                for line in f_new:
                    if (line[:3] == 'TER'.encode() or line[:4] == 'ATOM'.encode() or line[:6] == 'HETATM'.encode()):
                        f_outname.write(line)
                    elif i == 0 and line.startswith('MODEL'.encode()):
                        f_outname.write(line)

    # The following commented section is preferred for writing files, however without resetting the chain
    # names will reset them all, and if you reset the chains the double letters dont print correctly
    '''
    monomers = []
    for f in pdbs:
        monomers.append(bb.Molecule(f))

    M_multi = bb.Multimer()
    #print(labels)
    M_multi.load_list(monomers, labels)
    print('m multi res', M_multi.chain_names)
    M_multi.chain_names = labels
    print('m multi res', M_multi.chain_names)
    print(list(set(list(M_multi.data['chain']))))
    print(M_multi.data)
    
    M_multi.write_pdb(outname)
    '''


def curate(pdb, fasta, outdir="result", gap=10,
           verbose=False, include_hetatm=False):
    '''
    Take pdb file, create specific temp folder for working in, call fragment to split into chains
    and get the gap counts, if gap counts are greater than max allowed, raise exception, otherwise
    launch autopatch on the structure. Check that these chains still have the same starting resid as
    before and then recombine them all back into the final curated pdb file.

    :param pdb: Name of the pdb file of interest to fragment
    :type pdb: str
    :param fasta: Name of the fasta file corresponding to the pdb file of interest
    :type fasta: str
    :param outdir: path of the output directory for the overall run
    :type outdir: str
    :param gap: max gap allowed in the protein sequences to be patched
    :type gap: int
    :param verbose: Toggleable option to allow for printing of full output from Modeller to the
        terminal (default False)
    :type verbose: bool
    :param include_hetatm: Toggleable option to allow hetatms to pass through biobox
    :type include_hetatm: bool
    :returns: (path of the curated file, largest gap patched, geometry report). The
        report is what :func:`resdy.geometry.check_geometry` returns for the assembled
        structure, with the number of rebuilt residues added.
    :rtype: tuple
    '''
    pdb_tmp_name = f'tmp_{os.path.splitext(os.path.basename(pdb))[0]}'
    tmp_folder = os.path.join(outdir, pdb_tmp_name)
    if os.path.exists(tmp_folder):
        shutil.rmtree(tmp_folder)
    os.makedirs(tmp_folder, exist_ok=True)

    # divide structure in individual chains
    # heteroatoms are deliberately kept out of the per-chain files: Modeller has no use
    # for them and drops them from its output anyway. They are merged back after the
    # complex is reassembled, from the structure this was given.
    gap_count, chain_order = fragment(pdb, fasta, tmp_folder, include_hetatm=False)
    gaps_by_chain = {c: int(g[0]) for c, g in zip(chain_order, gap_count)}

    # if there is a gap in the sequence greater than a specified amount, raise an exception and don't patch with Modeller
    largest = np.max(gap_count[:, 2])
    if largest>gap:
        raise Exception(f"large gap detected ({largest} residues)")

    files = glob.glob(os.path.join(tmp_folder, "chain*fasta"))
    chains = []
    fouts = []
    modelled = set()
    for f in files:

        # attempt modelling
        fbasename = os.path.splitext(f)[0]
        chain_name = os.path.basename(fbasename)[len('chain'):]

        built = []
        if gaps_by_chain.get(chain_name, 1) == 0:
            # nothing is missing from this chain, so there is nothing for Modeller to
            # build. Running it anyway would rebuild the chain and return it in its own
            # frame, which for a complex means the chains no longer sit correctly with
            # respect to one another.
            foutname = f'{fbasename}.pdb'
        else:
            if verbose:
                foutname, built = autopatch(tmp_folder, fbasename, gap)
            else:
                with ShutUp():
                    foutname, built = autopatch(tmp_folder, fbasename, gap)

            if foutname == "":
                raise Exception("Autopatching failed.")

        # ensure that sequences of AA starts from the correct resid
        M_raw = bb.Molecule()
        M_raw.import_pdb(f'{fbasename}.pdb', include_hetatm=include_hetatm)
        startval_raw = M_raw.data["resid"].values
        M_curated = bb.Molecule()
        M_curated.import_pdb(foutname, include_hetatm=include_hetatm)
        startval_clean = M_curated.data["resid"].values
        if startval_raw[0] != startval_clean[0]:
            startval_clean += startval_raw[0] - startval_clean[0]
            M_curated.data["resid"] = startval_clean
            M_curated.write_pdb(foutname)

        if built:
            # translate the model's residue positions into the numbering the file now
            # carries, so that they can be recognised in the reassembled complex
            order = []
            for line in open(foutname):
                if line.startswith('ATOM'):
                    key = (chain_name, int(line[22:26]))
                    if key not in order:
                        order.append(key)
            modelled.update(order[i] for i in built if i < len(order))

        chains.append(chain_name)
        fouts.append(foutname)

    # reassemble complex in final directory
    os.makedirs(outdir, exist_ok=True)

    fname = f"{os.path.basename(pdb).split('.')[0]}.pdb"
    outname = os.path.join(outdir, fname)

    sorting_pairs = sorted(zip(chains, fouts), key=lambda cf: cf[0])
    chains, fouts = zip(*sorting_pairs)

    reassemble(fouts, chains, outname, outdir, include_hetatm=False)

    if include_hetatm:
        n_het = merge_heteroatoms(pdb, outname)
        if n_het:
            print(f'>> {n_het} heteroatom records carried through to '
                  f'{os.path.basename(outname)}')

    shutil.rmtree(tmp_folder)

    # Check the assembled result. Each chain was modelled on its own, so a rebuilt loop
    # knows nothing about the chains packed against it, and a clash there is invisible to
    # everything downstream: the measurements are taken happily either way.
    geometry, offending = check_geometry(outname, modelled_residues=modelled)
    geometry['n_modelled_residues'] = len(modelled)
    if geometry['n_clashes']:
        against_neighbour = 0 if not len(offending) else int(
            (offending['modelled'] & offending['inter_chain']).sum())
        print(f'>> {geometry["n_clashes"]} non-bonded clash(es) in {os.path.basename(outname)}, '
              f'closest {geometry["min_contact"]:.2f} A at {geometry["worst"]}')
        if against_neighbour:
            print(f'>> {against_neighbour} of them put a rebuilt residue into a '
                  f'neighbouring chain; treat this structure with caution')

    return outname, largest, geometry



if __name__ == "__main__":

    pdb = f"curate_PDB{os.sep}conformations{os.sep}1U8F-alt1A.pdb"
    fasta = f"curate_PDB{os.sep}conformations{os.sep}1U8F.fasta"

    outfolder = f"curate_PDB{os.sep}curated"
    gap = 10

    fname = curate(pdb, fasta, outdir=outfolder, gap=gap)
    print(f"generated {fname}")
