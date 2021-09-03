#! /usr/bin/python
# Download all files from KLIFS database, saving only the chain indicated in the database in folder "raw"
# clean downloaded files, patch them if necessary, and save result in folder "clean" (ready to be used as training set)
# 2 logfiles saved:
# - gap_data.txt (reports on how many missing residues the protein had)
# - patch_data.txt (reports on which files had to be patched with modeller, and whether the operation was successful)

import glob
import os
import subprocess
import numpy as np
import biobox as bb
import pandas as pd
import fileinput
from copy import deepcopy
# offers automatic automatic structure patching function
from curate_PDB.autopatch import autopatch
import re
import cleanfasta as cf
########################################
# control flags

download = True # download PDBs from databank (even if files are downloaded already)
download_fasta = True # download FASTA sequences associated to non-alternative structures
reprocess = True # force PDB processing and analysis even if results have already been logged
cutoff = 10 # autopatch cutoff (attempt adding residues to a protein if its gaps are all smaller than this amount of residues)

########################################

# load PDB file of choice, and return a biobox structure.
# if needed (outfile != ""), save the cleaned file in a new PDB. 
def get_data(pdb, chain):
        

# report on gaps on a given PDB file
# returns:
# - 4 elements list, [sequence gap cnt., sequence missing residues cnt., sequence max gap size, geometric gap count]]
# - biobox.Molecule, of loading was successful, nothing otherwise
        def analyze_protein(f):

                #attempt loading the protein (error: -2 if unloadable)
 
 

        # check backbone geometric split (error:-1 if N and C atoms count mismatch)
                M = bb.Molecule()

                M.import_pdb(f, include_hetatm=True)
                print('after')
                try:
                        c_cnt, _, _ = M.guess_chain_split(distance=3.5)
                        c_cnt -= 1

                except:
                        c_cnt = 1


        # check sequence split (note: we avoid residues with negative numbers)
                res = np.unique(M.data["resid"].values)
                res = res[res>0]
                missing = []
                patch = []
                cnt = [0, 0, 0, c_cnt]
                for r in range(np.min(res), np.max(res)+1):
                        if r in res:
                                gap = False
                                if len(patch) > 0:
                                        missing.append(deepcopy(patch))
                                        cnt[0] += 1
                                        cnt[1] += len(patch)
                                        if len(patch) > cnt[2]:
                                                cnt[2] = len(patch)

                                        patch = []
                   
                        else:
                                gap = True
                                patch.append(r)
                return cnt, M


##############################################################################
##############################################################################

# if folders containing raw (downloaded) and clean (ready for training) PDBs, create them

        if not os.path.exists("curate_PDB/raw"):
                os.mkdir("curate_PDB/raw")

        if not os.path.exists("curate_PDB/clean"):
                os.mkdir("curate_PDB/clean")


# load PDBs and save only chain of interest (pdbcode_chainname.pdb)
# replace False with True to launch download from PDB databank
        if download:
        # data columns stored in data are:
        #NAME FAMILY GROUPS PDB CHAIN ALTERNATE_MODEL SPECIES LIGAND PDB_IDENTIFIER ALLOSTERIC_NAME ALLOSTERIC_PDB DFG AC_HELIX
        
        # download and save PDBs (chain reported in database)
                fin = "%s.pdb"%pdb
                fout = "%s_%s.pdb"%(pdb, chain)

                conf_file = open('curate_PDB/conformations/' + fin)
                raw_file = open('curate_PDB/raw/' + fin, 'w')
                for line in conf_file:
                        raw_file.write(line)
                raw_file.close()

                # if file exists, skip
                if os.path.exists(fout):
                        print("skipping %s"%fout)
                        pass

                # else, download file, and get the subset out (only the chain indicated in KLIFS database)

                try:
                        M = bb.Molecule()
                        pdb_location = 'curate_PDB/raw/' + pdb + '.pdb'
                        M.import_pdb(pdb_location, include_hetatm=True)
                        _, idxs = M.atomselect(chain, "*", "*", get_index=True)
                        path = "curate_PDB/raw/" + fout
                        print(path)
                        M.write_pdb(path, index=idxs, split_struc=False)


                except Exception as e:
                        print("ERROR: %s"%e)
                        print("> issue with file %s"%fin)

                # remove the downloaded PDB file (we already saved what we need)
                os.remove('curate_PDB/raw/' + fin)


# download FASTA sequences of proteins of interest in "raw" folder (not *alt files)
        if download_fasta:
                for f in glob.glob("curate_PDB/raw/*pdb"):
                        
                        file_name = f[15:19]
                        print(file_name)
                        chain = f[-5]
                        print(f)
                        if "PATCHED" in f:
                                continue

                        file_name_url = file_name + '.' + chain
                        file_name_fasta = file_name + '_' + chain + '.fasta'
                        oldpwd=os.getcwd()
                        os.chdir('curate_PDB/raw')
                        web_url = "https://www.rcsb.org/fasta/chain/" + file_name_url + '/download'
                        subprocess.check_call("wget -O " + file_name_fasta + " " + web_url, shell=True)
                        new_name = (fin[:-4]) + '_' + chain + '.fasta'
                        print(file_name_fasta)
                        print(new_name)
                        os.rename(file_name_fasta, new_name)

                        os.chdir(oldpwd)
                        try:
                                cf.clean_fasta(new_name)
                        except:
                                pass

                

        

# clean and analyze downloaded structures: report on gaps and missing residues
# note: cleaned files are not saved (pass an additional parameter to the load_and_clean function to write them out
        if not os.path.exists("gap_data.txt") or reprocess:
                # result will contain output, 4 numbers per protein:
                #sequence gap cnt., sequence missing residues cnt., sequence max gap size, geometric gap count (using M.guess_chain_split())
                result = []
                files = np.array(glob.glob("curate_PDB/raw/*pdb"))
                patchstat = [] # 0 = not needed, 1 = successful, 2 = failed

                for k, f in enumerate(files):

                        if "PATCHED" in f:
                                continue

                        success = True
                        fout = os.path.basename(f)
                        print("\n%s: %s"%(k, fout.split(".")[0]))

                        # load protein and assess its structure
                        cnt, mol = analyze_protein(f)
                        print("> geom.gaps: %s. seq.gaps: %s. seq.missing resid: %s. seq.largest gap: %s"%(cnt[3], cnt[0], cnt[1], cnt[2]))

                        # if a small amount of geometric gaps are present (or C-N atomcount mismatch), send the structure to patching, and re-analyze result
                        if cnt[3] > 0 or cnt[3] == -1 or cnt[0] > 0:
                                # call patching code in autopatch module

                                f_patched = autopatch(f.split(".")[0], cutoff)
                                # test if patching has been successful (empty string returned = failure)
                                if f_patched != "":
                                        cnt2, mol2 = analyze_protein(f_patched)
                                        print(">> PATCH: geom.gaps: %s. seq.gaps: %s. seq.missing resid: %s. seq.largest gap: %s"%(cnt2[3], cnt2[0], cnt2[1], cnt2[2]))
                                        # use data of patched molecule, save patched molecule, edit output name to indicate molecule is patched
                                        cnt = cnt2[:]
                                        mol = deepcopy(mol2)
                                        fout = "%s_patched.pdb"%fout.split(".")[0]
                                        patchstat.append(1)
                                        pdbchain = pdb + chain



                                else:
                                        print(">> PATCH: failed, keeping original results in file %s"%f)
                                        patchstat.append(2)
                                        success = False


                        elif cnt[3] == 0:
                                print(">> PATCH: not needed")
                                patchstat.append(0) # if no gap is present

                        else:
                                print(">> PATCH: not applicable (molecule loading failed)")
                                patchstat.append(2) # if protein loading failed
                                success = False

                        result.append(cnt)

                # write clean PDB in "clean" folder (unless protein loading failed)
                        try:
                                if success:
                                        print(">> SAVING PROTEIN in curate_PDB/clean/%s"%fout)
                                        mol.write_pdb("curate_PDB/clean/%s"%fout, split_struc=False)
                                        
                                else:
                                        print(">> protein not saved (patching failed)")
                        except:
                                print(">> protein not saved (writing error)")
                                continue

                print(files)


                list_of_files = files.tolist()

                for f in list_of_files:
                        if re.search('PATCHED', f):
                                list_of_files.remove(f)

                files = np.array(list_of_files)
                print(files)

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


##############################################################################
##############################################################################

### report on pdb gaps stats, and success of patching ###

        print("\n")
        print("%s protein chains processed"%(len(result)))
        print("%s proteins not loadable"%(np.sum(result[:, 0] == -2)))

        # note: we want to use as dataset all proteins having patchstat equal to 0 or 1
        print("\nPatching:")
        print("   %s not needed"%np.sum(patchstat == 0))
        print("   %s successful"%np.sum(patchstat == 1))
        print("   %s failed"%np.sum(patchstat == 2))

        print("\nSequence gaps count:")
        for n in np.unique(result[:, 0]):
                if n == -2:
                        continue
                else:
                        print("   %s proteins with %s gaps"%(np.sum(result[:, 0] == n), n))


        print("\nGeometry gaps count:")
        for n in np.unique(result[:, 3]):

                if n == -2:
                        continue
                elif n == -1:
                        print("   %s proteins C-N atoms mismatch"%(np.sum(result[:, 3] == n)))
                        continue
                else:
                        print("   %s proteins with %s gaps"%(np.sum(result[:, 3] == n), n))


        print("\nSequence maximal gap size:")
        for n in np.unique(result[:, 2]):

                if n == -2:
                        continue
                else:
                        print("   %s proteins with %s gap size"%(np.sum(result[:, 2] == n), n))

        print(n)
        return n

