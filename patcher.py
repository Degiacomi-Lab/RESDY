'''
Download files, patch them if necessary, and save result in folder "clean"
2 logfiles saved:
- gap_data.txt (reports on how many missing residues the protein had)
- patch_data.txt (reports on which files had to be patched with modeller, and whether the operation was successful)
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
import pandas as pd

from modeller import *
from modeller.automodel import *

from helper import get_download_tool, ShutUp


def autopatch(fbasename, gap_cutoff=8):

    print('>> modelling missing residues')
    #pdb_out = "%s_PATCHED.pdb"%fbasename; the output pdb file name (if successful, empty otherwise)
    pdb_out = ''
    seq_name = ''
    try:

        _pdb2seq(fbasename)
        _fasta2pir(fbasename)
        seq_name = _full_align(fbasename)
        _trim_align("alignment.seg.ali")
        patch_status = _gap_check("trimmed_align.ali", gap_cutoff)

        if patch_status == "yes":
            pdb_out = _patch_model(fbasename, seq_name)
        else:
            pdb_out = ""

    except Exception as e:
        raise Exception(f">> ERROR on autopatching the structure: {e}") from e

    finally:

        autopatch_files = [f'{fbasename}.seq', f'{fbasename}.pir', 'alignment.seg',
                'alignment.seg.ali', 'trimmed_align.ali', 'family.mat']
        autopatch_files.extend(glob.glob(f'{seq_name}.ini'))
        autopatch_files.extend(glob.glob(f'{seq_name}.rsr'))
        autopatch_files.extend(glob.glob(f'{seq_name}.sch'))
        autopatch_files.extend(glob.glob(f'{seq_name}.V*'))
        autopatch_files.extend(glob.glob(f'{seq_name}.D*'))

        for patch_file in autopatch_files:
            m = os.path.join(os.getcwd(), patch_file)
            try:
                os.remove(m)
            except Exception as e:
                print(f"cannot remove {m}, error: {e}; continuing...")
                continue

    return pdb_out

#autopatch step 1a. pir format of AA from pdb
def _pdb2seq(fbasename):
    env = Environ()
    env.io.two_char_chain = True
    mdl = Model(env, file=fbasename)
    aln = Alignment(env)
    aln.append_model(mdl, align_codes=fbasename)
    aln.write(file=fbasename+'.seq')

#autopatch step 1b. pir from complete AA fasta
def _fasta2pir(fbasename):
    env = Environ()
    env.io.two_char_chain = True  # TODO check locations of these to see if they do anything
    a = Alignment(env, file=fbasename+".fasta", alignment_format='FASTA')
    a.write(file=fbasename+'.pir', alignment_format='PIR')

#autopatch step 2. add sequence name to 2nd line; copy the pir contents and structure info into alignment.seg; align sequences and generate model
def _full_align(fbasename):
    pir_fname = fbasename+'.pir'
    seq_fname = fbasename+'.seq'
    f = open(pir_fname, "r")
    f1 = f.readlines()
    f.close()
    #f1 = [x.rstrip() for x in f1]
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

    if sys.platform == "win32":
        my_cmd_a = f'type {pir_fname} {seq_fname} > alignment.seg'
    else:
        my_cmd_a = f'cat {pir_fname} {seq_fname} > alignment.seg'

    os.system(my_cmd_a)

    env = Environ()
    env.io.two_char_chain = True  # TODO check locations of these to see if they do anything
    env.io.atom_files_directory = ['.', f'..{os.sep}atom_files']
    a = AutoModel(env,
                  # file with template codes and target sequence
                  alnfile  = 'alignment.seg',
                  # PDB codes of the templates
                  knowns   = fbasename,
                  # code of the target
                  sequence = seq_name)
    a.auto_align() # get an automatic alignment (alignment.seg.ali)
    return seq_name

# autopatch step 3. trim the alignment by removing gaps for missing residues at the termini of the structure
def _trim_align(align_file):
    align_file = "alignment.seg.ali"
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
    f = open("trimmed_align.ali", 'w')
    f.writelines(struc_sec)
    f.writelines(seq_sec)
    f.close()

#autopatch step 4. Check if any gap is more than cutoff length in the trimmed_align.ali and if so set patch_status = "no"
def _gap_check(align_file, gap_cutoff):
    patch_status = "yes"
    align_file = "trimmed_align.ali"
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
            print(f">> struture not patched. Long sequence gap: {str(gap_len)}")
            patch_status = "no"

    return patch_status

#autopatch step 5. build missing residues
def _patch_model(fbasename, seq_name):
    print(">> patching model...")
    log.verbose()
    env = Environ()
    env.io.two_char_chain = True  # TODO check locations of these to see if they do anything
    env.io.atom_files_directory = ['.', f'..{os.sep}atom_files']
    a = AutoModel(env,
                  # file with template codes and target sequence
                  alnfile  = 'trimmed_align.ali',
                  # PDB codes of the templates
                  knowns   = fbasename,
                  # code of the target
                  sequence = seq_name,
                  assess_methods = (assess.DOPE, assess.GA341))
    a.md_level = refine.fast #very_fast, fast, slow, very_slow, slow_large, refine
    #repeat whole cycle twice and do not stop unless obj. func > 1e6
    #a.repeat_optimization = 2
    a.max_molpdf = 1e6
    a.make()

    pdb_out = f"{fbasename}_PATCHED.pdb"

    os.rename(f'{seq_name}.B99990001.pdb', pdb_out)

    return pdb_out


def analyze_protein(M):
    '''
    look for gaps in the sequence and return 4 elements list:
    [number of gaps, number of missing residues, largest sequence gap]]
    '''

    res = np.unique(M.data["resid"].values)
    missing = []
    patch = []
    cnt = [0, 0, 0]
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


def fragment(pdb, fasta, outfolder="."):
    '''
    split a PDB file in individual chains
    split its associated FASTA file in FASTA of individual chains
    return information of gaps in protein structure (as per the function analyse_protein)
    '''

    if not os.path.exists(outfolder):
        os.mkdir(outfolder)

    #split PDB file in chains using biobox
    M = bb.Molecule(pdb)
    chains = np.unique(M.data["chain"].values)

    # Catch cases which have more than 62 chains which are not able to be worked with using PDB files
    if len(chains) > 62:
        with open(os.path.join(os.path.split(outfolder)[0], 'excessive_chains_files.txt'), 'a') as f1:
            f1.write(f'PDB: {pdb}; Chains: {str(chains)}')
        raise Exception(f'>> More than 62 individual chains present in the structure, not possible to curate a file for this number of chain names. File: {pdb}')

    # Catch cases which have double letter chain names, used to work out if worth including more functionality for these, have left future double letter functionality in incase this is useful in future.
    if any(len(chain) > 1 for chain in chains):
        with open(os.path.join(os.path.split(outfolder)[0], 'double_letter_file_count.txt'), 'a') as f1:
            f1.write(f'PDB: {pdb}; Chains: {", ".join(chains)}\n')
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
            #chain_info = [chain_raw_info[i].strip()[0] for i in range(len(chain_raw_info))]
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

    return np.array(gap_count)


def reassemble(pdbs, labels, outname, outdir):
    '''
    For each pdb file in pdbs go back through and take the chain name from the file
    name to replace the chain name back to a double letter. This reverses the
    temporary measure of changing to lower case single letter names to allow
    Modeller to patch before changing back.
    '''
    tmpfolder = os.path.join(outdir, "tmp")
    for pdb_file in pdbs:
        # take pdb_file and extract the chain name - need to find the exact format it needs to go back into that the code is expecting to reconvert it from
        dbletter = False
        thousand_chain = False
        chainname = os.path.splitext(os.path.basename(pdb_file))[0].split('_')[0][5:]
        if len(chainname) > 1:
            dbletter = True

        M_tmp = bb.Molecule()
        M_tmp.import_pdb(pdb_file, include_hetatm=True)
        num_res = len(set(list((M_tmp.data['resid']))))
        if num_res > 999:
            thousand_chain = True

        try:
            with fileinput.FileInput(pdb_file, inplace = True) as f:
                for line in f:
                    try:
                        #On lines with 'ATOM', 'TER' or 'HETATM'
                        #the current chain name is replaced with the new single letter chain name and adjusted to match the correct pdb file format
                        if ((line[:4] == 'ATOM') or (line[:3] == 'TER') or (line[:6] == 'HETATM')) and len(line) > 26:
                            if dbletter and thousand_chain:
                                raise Exception(f'>> Carbamylation code can\'t handle writing chains with more than '
                                                f'1000 residues and a double letter chain name at the same time. '
                                                f'Chain name: {dbletter}; Num resids: {num_res}')
                            elif dbletter and not thousand_chain:
                                resid = line[22:26].strip()
                                line = line[:21] + chainname + '%3s' % resid + line[26:]
                            elif not dbletter and thousand_chain:
                                resid = line[22:26].strip()
                                line = line[:21] + chainname + '%4s' % resid + line[26:]
                            else:
                                line = line[:21] + chainname + ' ' + line[23:]

                        print(line, end ='')

                    except:
                        print(line, end='')

        #If the protein fails, print error message with the error
        except Exception as e:
            raise Exception(f'Failed replacing chains for file {pdb_file}. Could not reassemble chains. Error: {e}')

    # Take all the files in pdbs and turn each into a biobox Molecule object and append this to a list of monomers
    # create a biobox multimer from the list of monomors and write a new pdb file combining all these
    # Slower solution but works
    with open(outname, 'ab') as f_outname:
        for i, (f, f_chain) in enumerate(zip(pdbs, labels)):
            T = bb.Molecule()
            T.import_pdb(f, include_hetatm=True)
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



def curate(pdb, fasta, outdir="result", gap=10, verbose=True):
    tmpfolder = os.path.join(outdir, "tmp")
    if os.path.exists(tmpfolder):
        shutil.rmtree(tmpfolder)
    os.makedirs(tmpfolder, exist_ok=True)

    # divide structure in individual chains
    gap_count = fragment(pdb, fasta, tmpfolder)

    # if there is a gap in the sequence greater than a specified amount, raise an exception and don't patch with Modeller
    largest = np.max(gap_count[:, 2])
    if largest>gap:
        raise Exception(f"large gap detected ({largest} residues)")

    #launch modeller on each individual chain  # TODO do we actually need to launch modeller if gaps are 0?
    files = glob.glob(os.path.join(tmpfolder, "chain*fasta"))
    chains = []
    fouts = []
    for f in files:

        # attempt modelling
        fbasename = os.path.splitext(f)[0]

        if verbose:
            foutname = autopatch(fbasename, gap)
        else:
            with ShutUp():
                foutname = autopatch(fbasename, gap)

        if foutname == "":
            raise Exception("Autopatching failed.")

        # ensure that sequences of AA starts from the correct resid
        M_raw = bb.Molecule(f"{fbasename}.pdb")
        startval_raw = M_raw.data["resid"].values
        M_curated = bb.Molecule(foutname)
        startval_clean = M_curated.data["resid"].values
        if startval_raw[0] != startval_clean[0]:
            startval_clean += startval_raw[0] - startval_clean[0]
            M_curated.data["resid"] = startval_clean
            M_curated.write_pdb(foutname)

        chains.append(fbasename.split('chain')[-1])
        fouts.append(foutname)

    # reassemble complex in final directory
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    fname = f"{os.path.basename(pdb).split('.')[0]}.pdb"
    outname = os.path.join(outdir, fname)

    sorting_pairs = sorted(zip(chains, fouts), key=lambda cf: cf[0])
    chains, fouts = zip(*sorting_pairs)

    reassemble(fouts, chains, outname, outdir)

    shutil.rmtree(tmpfolder)

    return outname

##############################################################################

if __name__ == "__main__":

    if True:

        pdb = f"curate_PDB{os.sep}conformations{os.sep}1U8F-alt1A.pdb"
        fasta = f"curate_PDB{os.sep}conformations{os.sep}1U8F.fasta"

        outfolder = f"curate_PDB{os.sep}curated"
        gap = 10

        #tmpfolder = "curate_PDB\\tmp"
        #fragment(pdb, fasta, tmpfolder)

        fname = curate(pdb, fasta, outdir=outfolder, gap=gap)
        print(f"generated {fname}")
        sys.exit()

    # fname should be a basename: expect to find both a .pdb and a .fasta file with that name
    fbasename = sys.argv[1]

    # gap should be an integer, if user fails to provide one, use a default of 10
    if len(sys.argv>2):
        try:
            gap = int(sys.argv[2])
        except:
            print("setting up default gap size=10")
            gap = 10

    foutname = autopatch(fbasename, gap)
    if foutname == "":
        print("autopatch failed")
    else:
        print(f"saved patched file {foutname}")
