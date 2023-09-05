# Download files, patch them if necessary, and save result in folder "clean" (ready to be used as training set)
# 2 logfiles saved:
# - gap_data.txt (reports on how many missing residues the protein had)
# - patch_data.txt (reports on which files had to be patched with modeller, and whether the operation was successful)

import glob
import os
import sys
import re
from textwrap import wrap
import shutil

import numpy as np
import biobox as bb
from copy import deepcopy

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
        raise Exception(">> ERROR: %s"%e)

    myfiles = ['%s.seq'%fbasename, '%s.pir'%fbasename, 'alignment.seg',
               'alignment.seg.ali', 'trimmed_align.ali', 'family.mat']
    myfiles.extend(glob.glob('*.ini'))
    myfiles.extend(glob.glob('*.rsr'))
    myfiles.extend(glob.glob('*.sch'))
    myfiles.extend(glob.glob('%s.V*'%seq_name))
    myfiles.extend(glob.glob('%s.D*'%seq_name))
    
    for mfile in myfiles:
        m = os.path.join(os.getcwd(), mfile)
        try:
            if sys.platform == "win32":
                os.remove(m)
            else:
                os.system("rm %s &> /dev/null"%m)
        except:
            print("cannot remove %s, continuing..."%m)
            continue    

    return pdb_out

#autopatch step 1a. pir format of AA from pdb
def _pdb2seq(fbasename):
    env = Environ()
    mdl = Model(env, file=fbasename)
    aln = Alignment(env)
    aln.append_model(mdl, align_codes=fbasename)
    aln.write(file=fbasename+'.seq')

#autopatch step 1b. pir from complete AA fasta 
def _fasta2pir(fbasename):
    env = Environ()
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
    for i in range(len(f1)):
        if ("P1;" in f1[i]):
            seq_name = f1[i][4:8]
        if ("sequence:" in f1[i]):
            seq_pos = i 
    P1 = ">P1;"+seq_name+"\n"
    seq_line = "sequence:"+seq_name+":::::::-1.00:-1.00\n"
    AA_block = ''.join([str(elem) for elem in f1[seq_pos+1:]])
    seq_block = P1 + seq_line + AA_block
    f = open(pir_fname, 'w')
    f.writelines(seq_block)
    f.close()

    if sys.platform == "win32":
        myCmd_A = 'type %s %s > alignment.seg'%(pir_fname, seq_fname)
    else:
        myCmd_A = 'cat %s %s > alignment.seg'%(pir_fname, seq_fname)
    
    os.system(myCmd_A)

    env = Environ()
    env.io.atom_files_directory = ['.', '..%satom_files'%(os.sep)]
    a = AutoModel(env,
                  # file with template codes and target sequence
                  alnfile  = 'alignment.seg',
                  # PDB codes of the templates
                  knowns   = fbasename,
                  # code of the target
                  sequence = seq_name) 
    a.auto_align() # get an automatic alignment (alignment.seg.ali)
    return seq_name

#autopatch step 3. trim the alignment by removing gaps for missing residues at the termini of the structure
def _trim_align(align_file):
    align_file = "alignment.seg.ali"
    f=open(align_file, "r")
    f1 = f.readlines()
    P1_pos = []

    #identify positions of different sequances (P1) blocks
    for i in range(len(f1)):
        if ("P1;" in f1[i]):
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
    for i in range(len(f1)):
        if ("P1;" in f1[i]): #identify positions of different sequances (P1) blocks
            P1_pos.append(i)
    f.close()
    sec_1 = P1_pos[0]
    sec_2 = P1_pos[1]

    AA_struc = ''.join([str(elem.rstrip("\n").rstrip("*")) for elem in f1[sec_1    +2:sec_2]])
    #check for gap_lengths in AA_struc
    struc_gaps = re.findall('[-]+', AA_struc)
    for gaps in range(len(struc_gaps)):
        gap_len = len(struc_gaps[gaps])
        if gap_len > gap_cutoff:
            print(">> struture not patched. Long sequence gap: "+ str(gap_len))
            patch_status = "no"

    return patch_status

#autopatch step 5. build missing residues
def _patch_model(fbasename, seq_name):
    print(">> patching model...")
    log.verbose()
    env = Environ()
    env.io.atom_files_directory = ['.', '..%satom_files'%(os.sep)]
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

    pdb_out = "%s_PATCHED.pdb"%fbasename
    
    if sys.platform == "win32":
        myCmd_mv = 'MOVE /Y %s.B99990001.pdb %s'%(seq_name, pdb_out)
    else:
        myCmd_mv = 'mv %s.B99990001.pdb %s'%(seq_name, pdb_out)
    
    os.system(myCmd_mv)
    
    return pdb_out

############################################

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
   
    #split PDB file in chains
    M = bb.Molecule(pdb)
    chains = np.unique(M.data["chain"].values)
    gap_count = []
    for c in chains:
        _, idxs = M.atomselect([c], "*", "*", get_index=True)
        M2 = M.get_subset(idxs)
        M2.write_pdb(os.path.join(outfolder, "chain%s.pdb"%c))
        gap_count.append(analyze_protein(M2))
        
    #split FASTA
    fin = open(fasta, "r")
    headers = [] # fasta headers
    fasta_chains = [] #chains associated with header
    sequences = [] #collection of sequences
    for line in fin:
        if ">" in line:
            headers.append(line)
            chain_rawinfo = line.split("|")[1][6:].split(",")
            chain_info = [chain_rawinfo[i].strip()[0] for i in range(len(chain_rawinfo))]
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
    for i in range(len(headers)):
        for c in fasta_chains[i]:
            if c not in chains:
                
                raise Exception("chain mismatch between PDB and FASTA. %s, %s"%(fasta_chains, chains))
            fout = open(os.path.join(outfolder, "chain%s.fasta"%c), "w")
            fout.write(headers[i])
            fout.writelines(sequences[i])
            fout.close()
            
    return np.array(gap_count)
            

def reassemble(pdbs, labels, outname):

    monomers = []
    for f in pdbs:
        monomers.append(bb.Molecule(f))

    M = bb.Multimer()
    M.load_list(monomers, labels)
    M.write_pdb(outname)


def curate(pdb, fasta, outdir="result", gap=10, verbose=True):
        
    tmpfolder = os.path.join(outdir, "tmp")
    if os.path.exists(tmpfolder):
        shutil.rmtree(tmpfolder)
        os.mkdir(tmpfolder)

    # divide structure in individual chains
    gap_count = fragment(pdb, fasta, tmpfolder)

    largest = np.max(gap_count[:, 2])
    if largest>gap:
        raise Exception("large gap detected (%s residues)"%largest)
                
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
        M_raw = bb.Molecule("%s.pdb"%fbasename)
        startval_raw = M_raw.data["resid"].values
        M_curated = bb.Molecule(foutname)
        startval_clean = M_curated.data["resid"].values
        if startval_raw[0] != startval_clean[0]:
            startval_clean += startval_raw[0] - startval_clean[0]
            M_curated.data["resid"] = startval_clean
            M_curated.write_pdb(foutname)
        
        chains.append(fbasename[-1])
        fouts.append(foutname)

    # reassemble complex in final directory
    if not os.path.exists(outdir):
        os.makedirs(outdir)
        
    fname = "%s.pdb"%os.path.basename(pdb).split(".")[0]
    outname = os.path.join(outdir, fname)
    
    """possible bug here (fixed by sorting the lists alphabetically)"""
    chains = sorted(chains)
    
    fouts = sorted(fouts, key=lambda x: x.split('_')[-2][-1])
    
    """possible bug here (fixed by sorting the lists alphabetically)"""
    
    reassemble(fouts, chains, outname)
    #TODO: check whether patching process caused clashing with lysine  
    shutil.rmtree(tmpfolder)
    
    return outname


'''
def correct_resid(pdb, chain):
    #
    #Fix any resid numbering issues that arise as a result of the autopatcher.
    #
    
    #opens each patched file in clean.
    cleanfiles = np.array(glob.glob(os.path.join("curate_PDB", "clean", "*.pdb")))

    for cleanfile in cleanfiles:
        name = pdb + '_' + chain + '_patched'
        
        if (cleanfile[17:-4]) == name:
            print('> fixing patched file resid')
            try:
                #Next it gets the first resid.
                M = bb.Molecule()
                M.import_pdb(cleanfile)
                cleandf = M.data
                cleanresid = cleandf.at[0, 'resid']

                #Next it opens the corresponding raw file.
                rawfiles = np.array(glob.glob(os.path.join("curate_PDB", "raw", "*.pdb")))
                for rawfile in rawfiles:
                    name2 = pdb + '_' + chain

                    #It opens it in biobox and gets the first resid
                    if rawfile[15:-4] == name2:
                        S = bb.Molecule()
                        S.import_pdb(rawfile)
                        rawdf = S.data
                        rawresid = rawdf.at[0, 'resid']
                        
                        #If they aren't the same, it shifts every resid in the patched file so they match and writes a new pdb file.
                        if cleanresid != rawresid:

                            cleandf['resid'] = cleandf['resid'] + rawresid - 1
                            M.write_pdb(cleanfile)
                        else:
                            return
                        
            except Exception as e:
                raise Exception('Failed correcting resid. %s'%e)
                
    return
'''        

'''
def patch_structure(pdb, chain, download=True, download_fasta=True, gap_cutoff=8):
        
        #load PDB file of choice, and return a biobox structure.
        #if needed (outfile != ""), save the cleaned file in a new PDB. 
        
    
        # if folders containing raw (downloaded) and clean (ready for training) PDBs, create them
        if not os.path.exists("curate_PDB"):
                os.mkdir("curate_PDB")

        if not os.path.exists(os.path.join("curate_PDB", "raw")):
                os.mkdir(os.path.join("curate_PDB", "raw"))

        if not os.path.exists(os.path.join("curate_PDB", "clean")):
                os.mkdir(os.path.join("curate_PDB", "clean"))

        # load PDBs and save only chain of interest (pdbcode_chainname.pdb)
        # replace False with True to launch download from PDB databank
        if download:
       
                # data columns stored in data are:
                #NAME FAMILY GROUPS PDB CHAIN ALTERNATE_MODEL SPECIES LIGAND PDB_IDENTIFIER ALLOSTERIC_NAME ALLOSTERIC_PDB DFG AC_HELIX
        
                # get PDBs (chain reported in database)
                fin = "%s.pdb"%pdb
                fout = "%s_%s.pdb"%(pdb, chain)
                oldpwd = os.getcwd()
                
                conf_file = open(os.path.join("curate_PDB", "conformations", fin))                
                raw_file = open(os.path.join("curate_PDB", "raw", fin), 'w')
                for line in conf_file:
                        raw_file.write(line)
                        
                raw_file.close()

                # if file does not exist, download file, and get the subset out
                # (only the chain indicated in KLIFS database)
                if not os.path.exists(fout):

                    try:
                            M = bb.Molecule()
                            pdb_location = os.path.join("curate_PDB", "raw", fin)
                            M.import_pdb(pdb_location, include_hetatm=True)
                            _, idxs = M.atomselect(chain, "*", "*", get_index=True)
                            path = os.path.join("curate_PDB", "raw", fout)
          
                            M.write_pdb(path, index=idxs, split_struc=False)
          
                    except Exception as e:
                            raise Exception("File %s: %s"%(fin, e))
    
                # remove the downloaded PDB file (we already saved the chain subset we need)
                os.remove(os.path.join("curate_PDB", "raw", fin))

        # download FASTA sequences of proteins of interest in "raw" folder (not *alt files)
        if download_fasta:
                for f in glob.glob(os.path.join("curate_PDB", "raw", "*pdb")):

                    try:

                        file_name = f[15:19]
                        chain = f[-5]

                        if "PATCHED" in f:
                            continue

                        file_name_url = file_name + '.' + chain
                        file_name_fasta = file_name + '_' + chain + '.fasta'
                        web_url = "https://www.rcsb.org/fasta/chain/" + file_name_url + '/download'
                        
                    
                        os.chdir(os.path.join("curate_PDB", "raw"))
                        if sys.platform == "win32":
                            line = "curl -s -o " + file_name_fasta + " " + web_url
                        else:
                            line = "wget -O " + file_name_fasta + " " + web_url
                        
                        if os.path.exists(file_name_fasta):
                            os.remove(file_name_fasta)
                        
                        subprocess.check_call(line, shell=True)

                        new_name = (fin[:-4]) + '_' + chain + '.fasta'
                        
                        os.rename(file_name_fasta, new_name)
                    
                        os.chdir(oldpwd)
                        clean_fasta(new_name)
                        
                    except:
                        os.chdir(oldpwd)
                        pass
            
        # clean and analyze downloaded structures: report on gaps and missing residues
        # note: cleaned files are not saved (pass an additional parameter to the load_and_clean function to write them out
        if True:#not os.path.exists("gap_data.txt"):
                print("> Sequence gap analysis")
                # result will contain output, 4 numbers per protein:
                #sequence gap cnt., sequence missing residues cnt., sequence max gap size, geometric gap count (using M.guess_chain_split())
                result = []
                files = np.array(glob.glob(os.path.join("curate_PDB", "raw", "*pdb")))
                patchstat = [] # 0 = not needed, 1 = successful, 2 = failed

                for k, f in enumerate(files):

                        if "PATCHED" in f:
                                continue

                        success = True
                        fout = os.path.basename(f)
                        print("> working on %s"%fout.split(".")[0])

                        # load protein and assess its structure
                        cnt, mol = analyze_protein(f)
                        print(">> geom.gaps: %s. seq.gaps: %s. seq.missing resid: %s. seq.largest gap: %s"%(cnt[3], cnt[0], cnt[1], cnt[2]))

                        # if a small amount of geometric gaps are present (or C-N atomcount mismatch), send the structure to patching, and re-analyze result
                        if cnt[3] > 0 or cnt[3] == -1 or cnt[0] > 0:
                                # call patching code in autopatch module

                                f_patched = autopatch(f.split(".")[0], gap_cutoff)
                                # test if patching has been successful (empty string returned = failure)
                                if f_patched != "":
                                        cnt2, mol2 = analyze_protein(f_patched)
                                        print(">> geom.gaps: %s. seq.gaps: %s. seq.missing resid: %s. seq.largest gap: %s"%(cnt2[3], cnt2[0], cnt2[1], cnt2[2]))
                                        # use data of patched molecule, save patched molecule, edit output name to indicate molecule is patched
                                        cnt = cnt2[:]
                                        mol = deepcopy(mol2)
                                        fout = "%s_patched.pdb"%fout.split(".")[0]
                                        patchstat.append(1)
                                        #pdbchain = pdb + chain

                                else:
                                        print(">> Patching failed, keeping original results in file %s"%f)
                                        patchstat.append(2)
                                        success = False

                        elif cnt[3] == 0:
                                print(">> Patching not needed")
                                patchstat.append(0) # if no gap is present

                        else:
                                print(">> Patching not applicable (molecule loading failed)")
                                patchstat.append(2) # if protein loading failed
                                success = False

                        result.append(cnt)

                        # write clean PDB in "clean" folder (unless protein loading failed)
                        try:
                                if success:
                                        foutname = os.path.join("curate_PDB","clean", fout)
                                        print(">> Saving patched protein in %s"%foutname)
                                        mol.write_pdb(foutname, split_struc=False)
                                        
                                else:
                                        print(">> protein not saved (patching failed)")
                        except:
                                print(">> protein not saved (writing error)")
                                continue

                list_of_files = files.tolist()

                for f in list_of_files:
                        if re.search('PATCHED', f):
                                list_of_files.remove(f)

                files = np.array(list_of_files)

                result = np.array(result)
                patchstat = np.array(patchstat)

                # report on whether patching was needed and, if so, successful
                outdata = np.concatenate((np.array([files]).T, result), axis=1)

                np.savetxt("gap_data.txt", outdata, fmt="%s")
                np.savetxt("patch_data.txt", patchstat)

        else:
                # load precalculated gap and patch data for statistics
                indata = np.loadtxt("gap_data.txt", dtype=str)
                files = indata[:, 0]
                result = indata[:, 1:].astype(int)
                patchstat = np.loadtxt("patch_data.txt")

        return np.max(result[:, 2])
    
'''
'''
def assemble_multimer(gap_dict, pdb_code, list_chains):
    
    #get_data break the protein up into chains.
    #This code reassembles the protein into a multimer from the chains in the clean folder given the pdb code and chains the protein consists of.
    
    
    #makes sure there is an assembled folder.
    #rc.rename_chains(pdb_code)
    if not os.path.exists("assembled"):
        os.mkdir("assembled")

    list_to_remove = []
    print('> Assembling')
    
    try:
        name_of_assembly = pdb_code + '_assembled.pdb'
        Multi = bb.Multimer()

        #opens the pdb file for each chain and appends to Multi
        for chain in list_chains:
            name_of_pdb_file = pdb_code + '_' + chain + '.pdb'
            patched_pdb_file = pdb_code + '_' + chain + '_' + 'patched.pdb'
            
            M = bb.Molecule()
            
            try:
                path = os.path.join("curate_PDB", "clean", patched_pdb_file)
                M.import_pdb(path, include_hetatm=True)
                Multi.append(M, chain)

            except:
                try:
                    path = os.path.join("curate_PDB", "clean", name_of_pdb_file)
                    M.import_pdb(path, include_hetatm=True)
                    Multi.append(M, chain)

                except Exception as e:
                    print("%s"%e)
                    continue

        #write Multi as a single .pdb file
        path = os.path.join("assembled", name_of_assembly)
        Multi.write_pdb(path)
        list_to_remove = check_missing_chains(gap_dict, pdb_code, list_chains)

        with fileinput.FileInput(path, inplace = True) as f:
            for line in f:
                line = line.replace("TER","")
                print(line, end ='') 

    except Exception as e:
        raise Exception("failed to assemble. %s"%e)

    return list_to_remove
'''
'''
def check_missing_chains(gap_dict, pdb_code, list_chains):
    
    #These next two modules (as well as one in rename_chains) sort out issues that arise when a chain fails the autopatcher.
    #In particular they do two things:
    #- Make sure the chain naming is correct.
    #- Mark any points that are close to the missing chains to be removed.
    

    clean_chains = []
    failed_chains =[]

    #Firstly the chains that passed the autopatch (i.e. those in the clean folder) are appended to a list.
    for f in glob.glob(os.path.join("curate_PDB","clean","*.pdb")):

        try:
            if f[-12:] == '_patched.pdb':
                clean_chain = f[-13]
            else:
                clean_chain = f[-5]
            clean_chains.append(clean_chain)

        except Exception as e:
            print('>> Issue identifying missing chains. %s'%e)
            continue

    #If any chain is not in clean but is in the original protein it is appended to another list (failed_chains)
    for chain in list_chains:
        if chain not in clean_chains:
            failed_chains.append(chain)

    #If any chains have failed, rename_chains_assembled and note_residues_near_missing _chain sort out chain naming issues and remove residues close to the missing chains respectively.
    if len(failed_chains) != 0:
        print('>> Fixing chain issue....')
        try:
            rename_chains_assembled(failed_chains, list_chains)
            list_to_remove = note_residues_near_missing_chain(gap_dict, pdb_code, failed_chains)
        except Exception as e:
            print('Error dealing with chain failure- chain labelling may be incorrect in data. %s'%e)
            return
        #print('Chain issue fixed.')
    
    else:
        print('>> no missing chains!')
        list_to_remove = list()

    return list_to_remove
'''
'''
def rename_chains_assembled(failed_chains, list_chains):
    
    #Rename the chains in an assembled structure.
    #This is needed as if any of the chains fail the autopatch, the naming of any chains further down the alphabet shifts.
    #For example if the original chains were ABC and B failed then in the assembled the naming would be AB whereas this program renames B to C.
    
    
    #create a copy of the list of chains (a list of all the chains including those which have failed).
    dict_replacements = dict()
    list_chains_updated = list_chains.copy()

    for chain in list_chains:

       #If any chain has failed, it is removed from the updated copy.
        if chain in failed_chains:

            list_chains_updated.remove(chain)

    #Next, it checks if the position of any chains has changed in the updated list.
    for chain in list_chains_updated:
        
        new_index = list_chains_updated.index(chain)
        replacing = list_chains[new_index]

        #If it has, the dictionary is updated to include the chain name to be replaced and the chain name which is going to replace it.
        if replacing != chain:
            dict_replacements[replacing] = chain


    replace_list = list(dict_replacements.keys())

    #It then goes and opens the assembled file.
    files = list((glob.glob(os.path.join("assembled", "*.pdb"))))
    for file in files:
        
        #It then rewrites the file with the replacement chain name.
        try:
            with fileinput.FileInput(file, inplace = True) as f:
                for line in f:
                    try:
                        if (line[:4] == 'ATOM') or (line[:3] == 'TER'):
                            chain_name = line[21]
                            if chain_name in replace_list:
                                replacement_chain_name = dict_replacements.get(chain_name)
                                line = line[:21] + replacement_chain_name + line[22:]
                                print(line, end ='')
                            else:
                                print(line, end = '')
                        else:
                            print(line, end='')
                    except:
                        print(line, end='')

        except Exception as e:
            print('>> %s'%e)
                  
    return
'''

'''
def note_residues_near_missing_chain(gap_dict, pdb_code, failed_chains):
    
    #Works as follows:
    #- Opens pdb file in curate_PDB/conformations/ with Biobox
    #- Check each lysine to see if they are near failed chain.
    #- If they are add chain_resid to list (list_to_remove)
    #- Later in main those in list_to_remove are removed from the results.
    
    
    path = os.path.join("curate_PDB", "conformations", pdb_code)
    list_of_chains = list()
    list_of_resid = list()
    list_to_remove = []
    
    print('> Breaking up molecule...')
    try:

        #open file in biobox
        M = bb.Molecule()
        M.import_pdb(path, include_hetatm=True)
        df = M.data

        #note of all the lysine index/coordinates
        lys_coords, lys_idx = M.atomselect('*','LYS', 'CA', use_resname=True, get_index=True)

        #note of each lysines chain and resid.
        for entry in lys_idx:
            chain = df.at[entry, 'chain']
            list_of_chains.append(chain)

        for entry in lys_idx:
            resid = df.at[entry, 'resid']
            list_of_resid.append(resid)

        #get coordinates and index for all atoms
        all_coords, idx = M.atomselect('*','*','*', get_index=True)

    except Exception as e:
        raise Exception('Failed identifying residues near missing chain %s'%e)

    #goes through all the lysines
    for failedchain in failed_chains:
        gapsize = gap_dict.get(failedchain)
        gapsize = int(gapsize)
        cutoff = 0.5*gapsize + 3.5

        for j in range(len(lys_coords)):

            #Does not bother if they are on a chain which failed the autopatch
            #as they won't be in the final structure anyway.
            chain_lys = list_of_chains[j]
            if chain_lys == failed_chains:
                continue
            
            else:

                try:
                    #Work out the distance between all the atoms and each lysine.
                    for i in range(len(all_coords)):
                        x_dist = (((lys_coords[j])[0] - (all_coords[i])[0])**2)
                        y_dist = (((lys_coords[j])[1] - (all_coords[i])[1])**2)
                        z_dist = (((lys_coords[j])[2] - (all_coords[i])[2])**2)
                        distance = np.sqrt(x_dist + y_dist + z_dist)

                        #If the distance is less than 10 A AND the atoms chain is one of the failed chains,
                        #the lysine's Chain_Resid is appended to the list_to_remove.
                        if distance < cutoff:
                            index_of_aa = idx[i]
                            aa_chain = df.at[index_of_aa, 'chain']

                            if aa_chain in failed_chains:
                                chain_resid = str(list_of_chains[j]) + str(list_of_resid[j])
                                list_to_remove.append(chain_resid)

                except Exception as e:
                    print('> Error %s'%e)
        
    return list_to_remove
'''
'''
def clean_fasta(fname):
    #clean the fasta file, getting rid of any non-canonical AAs that may cause an issue.
    
    #open the fasta file
    #path = os.path.join("curate_PDB","raw", new_name)
    #rewrite it replacing the amino acids that may cause an issue.
    with fileinput.FileInput(fname, inplace = True) as f:
        for line in f:
            if("KCX" in line):
                line = line.replace('(KCX)', 'K')
                print(line, end = '')
            if("MSE" in line):
                line = line.replace('(MSE)', 'M')
                print(line, end = '')
            else:
                print(line, end = '')

    return
'''
    
##############################################################################

if __name__ == "__main__":

    if True:

        pdb = "curate_PDB\\conformations\\1U8F-alt1A.pdb"
        fasta = "curate_PDB\\conformations\\1U8F.fasta"
        
        outfolder = "curate_PDB\\curated"
        gap = 10
        
        #tmpfolder = "curate_PDB\\tmp"
        #fragment(pdb, fasta, tmpfolder) 
        
        fname = curate(pdb, fasta, outfolder=outfolder, gap=gap)
        print("generated %s"%fname)
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
        print("saved patched file %s"%foutname)
        
