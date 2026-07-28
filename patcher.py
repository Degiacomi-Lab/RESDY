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

from modeller import *
from modeller.automodel import *

from helper import get_download_tool, ShutUp


def autopatch(fbasename, gap_cutoff=8):

    print('>> modelling missing residues')
    #pdb_out = "%s_PATCHED.pdb"%fbasename; the output pdb file name (if successful, empty otherwise)
    pdb_out = ""
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
        raise Exception(f">> ERROR: {e}") from e

    myfiles = [f'{fbasename}.seq', f'{fbasename}.pir', 'alignment.seg',
               'alignment.seg.ali', 'trimmed_align.ali', 'family.mat']
    myfiles.extend(glob.glob('*.ini'))
    myfiles.extend(glob.glob('*.rsr'))
    myfiles.extend(glob.glob('*.sch'))
    myfiles.extend(glob.glob(f'{seq_name}.V*'))
    myfiles.extend(glob.glob(f'{seq_name}.D*'))

    for mfile in myfiles:
        m = os.path.join(os.getcwd(), mfile)
        try:
            if sys.platform == "win32":
                os.remove(m)
            else:
                os.system(f"rm {m} &> /dev/null")
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

    if sys.platform == "win32":
        my_cmd_mv = f'MOVE /Y {seq_name}.B99990001.pdb {pdb_out}'
    else:
        my_cmd_mv = f'mv {seq_name}.B99990001.pdb {pdb_out}'

    os.system(my_cmd_mv)

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
    for r in range(1, np.max(res)+1):
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
    gap_count = []
    for c in chains:
        _, idxs = M.atomselect([c], '*', '*', get_index=True)
        M_two = M.get_subset(idxs)
        M_two.write_pdb(os.path.join(outfolder, f"chain{c}.pdb"))
        gap_count.append(analyze_protein(M_two))

    #split FASTA
    fin = open(fasta, 'r')
    fasta_headers = []
    fasta_chains = []
    sequences = []
    for line in fin:
        if ">" in line:
            fasta_headers.append(line)
            chain_raw_info = line.split("|")[1][6:].split(",")
            #chain_info = [chain_raw_info[i].strip()[0] for i in range(len(chain_raw_info))]
            if len(chain_raw_info[0]) == 1:
                chain_info = chain_raw_info
            elif '[' in chain_raw_info[0]:
                chain_info = [str(chain_raw_info[0].split('[')[0].strip())]
            else:
                chain_info = [chain_raw_info[i].strip()[0] for i in range(len(chain_raw_info))]
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


    # The following section commented out turns the fragmented files into the format specified by Modeller to work with double letter chain names
    # This should work with the parameter iodata.two_char_chain, however this doesn't appear to be working in this implementation
    # Therefore this code is left incase the parameter will start to work, new method is added below - GW_04.03.24
    '''
    # change the position of the chain names from single character format at 22 to double character 21-22 
    all_files = os.listdir(outfolder) 
    pdb_files = []
    for file in all_files:
        temp_parts = file.split('.')
        if file[0:5] == 'chain':
            if temp_parts[1] == 'pdb' and len(temp_parts[0][5:]) == 2:
                pdb_files.append(file)
        
    for file in pdb_files:
        # This will move the chain name 1 position to the left such that it will start at position 21 rather than 22 as it currently is
        # This allows it to work with the two_chain_char modification to allow double chain names to be read by Modeller
        # This style may be different for other programmes than Modeller
        # TODO: check if Modeller output of double chain names works for other structures.
        

        with fileinput.FileInput(os.path.join(outfolder, file), inplace = True) as f:
            for line in f:
                try:
                    #On lines with 'ATOM', 'TER' or 'HETATM' if the chain is in the auth_list
                    #the auth chain name is replaced with the RCSB chain name.
                    if (line[:4] == 'ATOM') or (line[:6] == 'HETATM'):
                            line = line[:20] + line[21:23] + ' ' + line[23:]
                            print(line, end ='')

                    else:
                        print(line, end='')
                except:
                    print(line, end='')
    '''

    # This next section is a replacement for the conversion to iodata.two_char_chain Modeller format
    # Instead convert the chain names back to single character chain names so Modeller can read this properly and doesnt convert all double letter chain names to A
    # A new function is added at the end of patching to convert the single chain names back to the corresponding double chain names - call protein.py function which does this already?
    # find the double chain name files and store in list to iterate through when converting the single character names
    all_files = os.listdir(outfolder)
    doubleletter_pdb_files = []
    for file in all_files:
        temp_parts = file.split('.')
        if file[0:5] == 'chain':
            # next take only files which are pdb and have 2 character chain names
            if temp_parts[1] == 'pdb' and len(temp_parts[0][5:]) == 2:
                doubleletter_pdb_files.append(file)

    # iterate through the files and for each one replace the chain name with a single letter lowercase version of the chain name: eg. CA would go to c
    used_sl_chain_names = []
    for file in doubleletter_pdb_files:
        db_chain_name = file.split('.')[0][5:]
        if db_chain_name.lower()[0] not in used_sl_chain_names:
            sl_chain_name = db_chain_name.lower()[0]
        elif db_chain_name.lower()[1] not in used_sl_chain_names:
            sl_chain_name = db_chain_name.lower()[1]
        else:
            numb_chains = [int(a) for a in used_sl_chain_names if a.isdigit()]
            if numb_chains: str(sl_chain_name = max(numb_chains) + 1)
            else: sl_chain_name = '0'

        try:
            temp_file_path = os.path.join(outfolder, file)
            with fileinput.FileInput(temp_file_path, inplace = True) as f:
                for line in f:
                    try:
                        #On lines with 'ATOM', 'TER' or 'HETATM'
                        #the current chain name is replaced with the new single letter chain name and adjusted to match the correct pdb file format
                        if (line[:4] == 'ATOM') or (line[:3] == 'TER') or (line[:6] == 'HETATM'):
                            new_chain_name = ' ' + sl_chain_name
                            line = line[:20] + new_chain_name + ' ' + line[23:]
                            print(line, end ='')
                        else:
                            print(line, end='')
                    except:
                        print(line, end='')

        #If the protein fails, print error message with the error
        except Exception as e:
            raise Exception(f'Failed replacing chains for file {file}. Could not convert double letter chain names while fragmenting. {e}')

    return np.array(gap_count)


def reassemble(pdbs, labels, outname, outdir):

    # For each pdb file in pdbs go back through and take the chain name from the file name to replace the chain name back to a double letter
    # This reverses the temporary measure of changing to lower case single letter names to allow Modeller to patch before changing back
    tmpfolder = os.path.join(outdir, "tmp")
    for pdb_file in pdbs:
        # take pdb_file and extract the chain name - need to find the exact format it needs to go back into that the code is expecting to reconvert it from
        dbletter = False
        chainname = pdb_file.split('.')[0].split('/')[-1].split('_')[0][5:]
        if len(chainname) > 1:
            dbletter = True

        try:
            if dbletter:
                with fileinput.FileInput(pdb_file, inplace = True) as f:
                    for line in f:
                        try:
                            #On lines with 'ATOM', 'TER' or 'HETATM'
                            #the current chain name is replaced with the new single letter chain name and adjusted to match the correct pdb file format
                            if (line[:4] == 'ATOM') or (line[:3] == 'TER') or (line[:6] == 'HETATM'):
                                line = line[:20] + chainname + '' + line[22:]
                                print(line, end ='')
                            else:
                                print(line, end='')
                        except:
                            print(line, end='')

        #If the protein fails, print error message with the error
        except Exception as e:
            raise Exception(f'Failed replacing chains for file {pdb_file}. Could not reassemble chains. {e}')

    # Take all the files in pdbs and turn each into a biobox Molecule object and append this to a list of monomers
    # create a biobox multimer from the list of monomors and write a new pdb file combining all these
    monomers = []
    for f in pdbs:
        monomers.append(bb.Molecule(f))

    M = bb.Multimer()
    M.load_list(monomers, labels)
    M.write_pdb(outname)


    '''
    # After the pdb file has been reassembled, with the current method of replacing double letter chain names with lower case single ones temporarily
    # need to convert back to the format it was given to the patcher.py script in
    # Iterate through the file line by line and check for lower case letter chain names, replace with upper case version + 'A'

    try:
        outfile_path = os.path.join(outfolder, pdb)
        with fileinput.FileInput(outfile_path, inplace = True) as f:
            for line in f:
                try:
                    #On lines with 'ATOM', 'TER' or 'HETATM'
                    #the current chain name is replaced with the new single letter chain name and adjusted to match the correct pdb file format
                    if (line[:4] == 'ATOM') or (line[:3] == 'TER') or (line[:6] == 'HETATM'):
                            current_chainname = line[21]
                            if len(current_chainname)
                            new_db_chainname = current_chainname.upper() + 'A'
                            line = line[:20] + new_db_chainname + ' ' + line[22:]
                            
                            print(line, end ='')
                            
                    else:
                        print(line, end='')
                except:
                    print(line, end='')

    #If the protein fails, print error message with the error 
    except Exception as e:
        raise Exception(f'Failed replacing chains for file {file}. %s'%e)
    '''



def curate(pdb, fasta, outdir="result", gap=10, verbose=True):
    tmpfolder = os.path.join(outdir, "tmp")
    if os.path.exists(tmpfolder):
        shutil.rmtree(tmpfolder)
        os.mkdir(tmpfolder)

    # divide structure in individual chains
    gap_count = fragment(pdb, fasta, tmpfolder)

    # if there is a gap in the sequence greater than a specified amount, raise an exception and don't patch with Modeller
    largest = np.max(gap_count[:, 2])
    if largest>gap:
        raise Exception(f"large gap detected ({largest} residues)")

    #launch modeller on each individual chain
    files = glob.glob(os.path.join(tmpfolder, "chain*fasta"))
    chains = []
    fouts = []
    for f in files:

        # attempt modelling
        fbasename = f.split(".")[0]

        if verbose:
            foutname = autopatch(fbasename, 10)
        else:
            with ShutUp:
                foutname = autopatch(fbasename, 10)

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

        # GW 18.03.24 - change chains appending to account for double letter chain names
        chains.append(fbasename.split('chain')[-1])
        fouts.append(foutname)

    # reassemble complex in final directory
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    fname = f"{os.path.basename(pdb).split('.')[0]}.pdb"
    outname = os.path.join(outdir, fname)

    #possible bug here (fixed by sorting the lists alphabetically)
    chains = sorted(chains)
    fouts = sorted(fouts, key=lambda x: x.split('_')[-2][-1])

    #possible bug here (fixed by sorting the lists alphabetically)

    reassemble(fouts, chains, outname, outdir)
    #TODO: check whether patching process caused clashing with lysine
    shutil.rmtree(tmpfolder)

    return outname

##############################################################################

if __name__ == "__main__":

    if True:

        pdb = "curate_PDB\\conformations\\1U8F-alt1A.pdb"
        fasta = "curate_PDB\\conformations\\1U8F.fasta"

        outfolder = "curate_PDB\\curated"
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
