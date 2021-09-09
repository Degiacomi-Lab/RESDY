import os
import subprocess
import re
import sys
import numpy as np
import glob
import fileinput

import pandas as pd
import biobox as bb

#OVERALL STRUCTURE...
#First downloads file into conformations folder
#Next cleans (i.e. removes heteroatoms that aren't metal ions)
#Next writes file for alt NMR confs
#Next writes file for alt AA confs


#This script downloads the file and removes all the hetereoatoms that aren't metal ions.
#The results is saved into a file with the structure *PDB code*-clean.pdb.
def clean_and_split_alt_conformations(pdb, done_pdbs):

    if not os.path.exists("curate_PDB"):
        os.mkdir("curate_PDB")

    if not os.path.exists("curate_PDB%sconformations"%os.sep):
        os.mkdir("curate_PDB%sconformations"%os.sep)

    def checks(pdb, done_pdbs):
        print('Checking')
        keep = True

#Firstly it checks whether they are in the log_file (the file saying what has already been done). If they are then the structure is not included.

        for entry in done_pdbs:

            if re.search(pdb, entry):
                keep = False
                print('\n')
                print('**********************************************************************')
                print(pdb + ' is in the log_file therefore will not be further investigated.')
                print('**********************************************************************')
                print('\n')
                if not os.path.exists('log_file.csv'):
                    f = open('log_file.csv', 'w')
                    f.write('PDB Code,Result')
                    f.write('\n')
                    f.write(pdb + ',Failed as pdb file is already in log_file- remove to continue.')
                    f.close()

                elif os.path.exists('log_file.csv'):
                    f = open('log_file.csv', 'a')
                    f.write('\n')
                    f.write(pdb + ',Failed as pdb file is already in log_file- remove to continue.')
                    f.close()
                break


#Next it downloads the relevant pdb file.

        if sys.platform == "win32":
            line = "curl -o %s.pdb https://files.rcsb.org/download/%s.pdb"%(pdb, pdb)
        else:
            line = "wget https://files.rcsb.org/download/" + pdb + '.pdb'
        subprocess.check_call(line, shell=True)
        
        
        f = open(pdb + '.pdb', 'r')
        list_of_metals = ['ZN', 'NI', 'CU', 'FE', 'MG', 'MN', 'NA', 'K', 'CA', 'CO', 'CL', 'MO']

#Here it cleans the file by only including lines if they include protein atoms or metal ions.
#Note- if you want to discard any files with non-metal ligands make keepp untrue and uncomment out the next ~18 lines

        for line in f:
                if line[:6] == 'HETATM':
                    words = line.split()
                    if words[3] in list_of_metals:
                        continue
                    elif words[2] == 'HOH':
                        continue
                    else:
                        keep = True
                        
                        #print('\n')
                        #print('****************************************************')
                        #print('NON-METAL HETATM FOUND- will not continue with this pdb file')
                        #print('****************************************************')
                        #print('\n')

                        #if not os.path.exists('log_file.csv'):
                            #f = open('log_file.csv', 'w')
                            #f.write('PDB Code,Result')
                            #f.write('\n')
                            #f.write(pdb + ',Failed due to presence of Hetatm-replace Hetatm with ATOM in .pdb file to include.')
                            #f.close()

                        #elif os.path.exists('log_file.csv'):
                            #f = open('log_file.csv', 'a')
                            #f.write('\n')
                            #f.write(pdb + ',Failed due to presence of Hetatm-replace Hetatm with ATOM in .pdb file to include.')
                            #f.close()
                            #break
                            
        f.close()
        os.remove(pdb + '.pdb')
                    
        return keep


    def clean(pdb):

        print('CLEANING')
        try:
            #Firstly it goes into curate_PDB/conformations and downloads the .pdb file.
            cwd = os.getcwd()
            #print(cwd)

            os.chdir("curate_PDB%sconformations"%os.sep)
            try:
                
                if sys.platform == "win32":
                    line = "curl -o %s.pdb https://files.rcsb.org/download/%s.pdb"%(pdb, pdb)
                else:
                    line = "wget https://files.rcsb.org/download/" + pdb + '.pdb'
                
                subprocess.check_call(line, shell=True)
                          
                
            except Exception as e:
                print('Error %s'%e)
                os.chdir(cwd)
            #Next it returns to the carbamylation folder.

            os.chdir(cwd)
            try:
                rename_chains(pdb)
            except Exception as e:
                print('Error %s'%e)
                print('Failure renaming chains')
                return('FAIL')

            
        except Exception as e:
            print('Error %s'%e)
            print('Print failure opening file.')
            return

        try:
            replace_mse(pdb)
        except:
            pass
    #Next it opens and starts reading the .pdb file and starts writing a new file with the ending '-clean.pdb'.

        list_of_metals = ['ZN', 'NI', 'CU', 'FE', 'MG', 'MN', 'NA', 'K', 'CA', 'CO', 'CL', 'MO']

        read_file_path = "curate_PDB%sconformations%s%s.pdb"%(os.sep, os.sep, pdb)

        read_file = open(read_file_path)
        write_file_path = "curate_PDB%sconformations%s%s-clean.pdb"%(os.sep, os.sep, pdb)
        write_file = open(write_file_path, 'w')

    #Next it writes the clean file, including HETATMs (if they are metal ions), all atoms and lines starting with TER and END.
        unknown_hetatm = False
        for line in read_file:
            if line[:6] == 'HETATM':
                words = line.split()
                if words[3] in list_of_metals:
                    write_file.write(line)
                    continue

            if line[:4] == 'ATOM':
                words = line.split()
                if words[2] != 'H':
                    write_file.write(line)
                    continue

            if line[:3] == 'TER':
                write_file.write(line)
                continue
            if line[:3] == 'END':
                write_file.write(line)
                continue
        write_file.close()

    #Lastly it removes the original pdb file as it is not needed anymore.

        read_file.close()
        os.remove(read_file_path)
        print('Clean ' + pdb)

        return pdb

    def replace_mse(pdb):
        print('Replacing MSE')
        files = np.array(glob.glob("curate_PDB%sconformations%s*pdb"%(os.sep, os.sep)))
        list_of_files = []
        for file in files:
            file_name = file[25:]
            if file_name[:4] == pdb:
                list_of_files.append(file_name)
        
#For each it then replaces any atoms beloning to KCX with LYS and HETATM with ATOM.

        for f in list_of_files:

            try:
                path = "curate_PDB%sconformations%s%s"%(os.sep, os.sep, f)

                with fileinput.FileInput(path, inplace = True) as f:
                    for line in f:
                        if("MSE" in line):
                            line = line.replace("HETATM","ATOM  ")
                            line = line.replace('MSE', 'MET')
                            line = line.replace('SE', ' S')
                            print(line, end ='') 
                        else:
                            print(line, end ='') 

            except:
                return
        return


    #This writes a new file for each alternate NMR structure.

    def split_struc_NMR(pdb):

    #Firstly it defines the name format which will be followed for each file.

        number = 1
        name = pdb + '-alt-' + str(number) + '.pdb'
        path = "curate_PDB%sconformations%s%s-clean.pdb"%(os.sep, os.sep, pdb)

    #Next it opens the file produced from the cleaning script and opens a new file to write in.

        f = open(path)

        endmdls = list()

    #Next it searches to see if there are any 'ENDMDL' statements in the folder (i.e. if there are multiple models).
        for line in f:
            if re.search('ENDMDL', line):
                endmdls.append(line)

        f.close()
    #If there are no 'ENDMDL' statements the clean file is renamed to suit the new format.

        if len(endmdls) == 0:
            path_rename = "curate_PDB%sconformations%s%s"%(os.sep, os.sep, name)
            os.rename(path, path_rename)
            print('NO ALTERNATE WHOLE STRUCTURES FOUND FOR ' + pdb)


    #If there are 'ENDMDL' statements a new file is written for each model.

        elif len(endmdls) != 0:
            
            print('ALTERNATE WHOLE STRUCTURE FOUND FOR ' + pdb)

    #First it opens the clean file in the conformations folder and opens a new folder to write in

            f = open(path)
            path_rename = "curate_PDB%sconformations%s%s"%(os.sep, os.sep, name)
            f_write = open(path_rename, 'w')

    #Next it writes the new file.
    #It includes every line until it gets to 'ENDMDL', where it opens a new file to write in.
    #The process stops when it gets to 'MASTER'
            for line in f:
                try:
                    line = str(line)

                    if re.search('END ', line):
                        f.close()
                        f_write.close()
                        os.remove(path)
                        os.remove(path_rename)
                        break

                    else:

                        if re.search('ENDMDL', line):
                            f_write.write(line)
                            f_write.close()
                            number = number + 1
                            name = pdb + '-alt-' + str(number) + '.pdb'
                            path_rename = "curate_PDB%sconformations%s%s"%(os.sep, os.sep, name)
                            f_write = open(path_rename, 'w')

                        else:
                            f_write.write(line)
                except:
                    continue
        
            f_write.close()
            f.close()
        return




    #Writes a new file for each alternative amino acid conforation present.

    def split_struc_alt_aa(pdb):

#Firstly it gets a list of all the .pdb files present in curate_PDB/conformations and selects those that belong to the pdb we are interested in.

        files = np.array(glob.glob("curate_PDB%sconformations%s*pdb"%s(os.sep, os.sep)))
        list_of_files = []
        for file in files:
            file_name = file[25:]
            if file_name[:4] == pdb:
                list_of_files.append(file_name)
        
        for f in list_of_files:

                ABC_list = ['A', 'B', 'C', 'D']
                ABC_dict = {"A": 0, "B": 0, "C": 0}

#Next it checks if there are any alternate amino acid conformations present (i.e. if line[16 == A, B or C]).

                for i in range(len(ABC_list)):
                    path = "curate_PDB%sconformations%s%s"%(os.sep, os.sep, f)
                    read_file = open(path)

                    for line in read_file:
                        if (line[:4] == 'ATOM') and (line[16] == ABC_list[i]):
                            ABC_dict[ABC_list[i]] = 1
                

                list_of_values = ABC_dict.values()

#If there aren't any for this file it moves on to the next.

                if 1 not in list_of_values:
                    continue
                else:
                    print('FOUND ALTERNATE AMINO ACID CONFORMATIONS FOR ' + pdb)

#If there are it rewrites a new file for all the As, Bs and Cs.

                for i in range(len(ABC_list)):
                    try:
                        if ABC_dict.get(ABC_list[i]) == 1:
                            name = f[:-6] + f[-5] + ABC_list[i] + '.pdb'

                            write_path = "curate_PDB%sconformations%s%s"%(os.sep, os.sep, name)
                            f_write = open(write_path, 'w')

                            target_letter = ABC_list[i]
                            non_target_letters = []
                            for letter in ABC_list:
                                if letter != target_letter:
                                    non_target_letters.append(letter)

                            read = open(path)
                            for line in read:
                                
                                if (line[:4] == 'ATOM') and (line[16] == target_letter):
                                    newline = line[:16] + ' ' + line[17:]
                                    f_write.write(newline)
                                    continue
                                if (line[:4] == 'ATOM') and (line[16] in non_target_letters):
                                    continue
                                if (line[:4] == 'ATOM'):
                                    f_write.write(line)
                                if (line[:6] == 'HETATM'):
                                    f_write.write(line)
                                if (line[:3] == 'TER'):
                                    f_write.write(line)
                            f_write.close()

                    except Exception as e:
                        print("Error %s"%e)
                        
#The original is then removed if it has been replaced.
                read_file.close()
                os.remove(path)
                
        return



#This script removes the carbamate from any structure

    def remove_kcx(pdb):

#Firstly it gets a list of files in curate_PDB/conformations that belong to the pdb of interest.

        files = np.array(glob.glob("curate_PDB%sconformations%s*pdb"%(os.sep, os.sep)))
        list_of_files = []
        for file in files:
            file_name = file[25:]
            if file_name[:4] == pdb:
                list_of_files.append(file_name)
        
        for file in list_of_files:
            path = "curate_PDB%sconformations%s%s"%(os.sep, os.sep, file)
            M = bb.Molecule()
            M.import_pdb(path, include_hetatm=True)

#Next it removes the atoms belonging to KCX from the structure.

            idxs = []
            pos ,idx = M.atomignore('*', 'KCX', 'CX', get_index=True, use_resname=True)
            A = M.get_subset(idxs=idx)
            pos, idx = A.atomignore('*', 'KCX', 'OQ1', get_index=True, use_resname=True)
            B = A.get_subset(idxs=idx)
            pos, idx = B.atomignore('*', 'KCX', 'OQ2', get_index=True, use_resname=True)
            C = B.get_subset(idxs=idx)
            try:
                C.write_pdb(path, split_struc=True)
            except:
                C.write_pdb(path, split_struc=False)


#For each it then replaces any atoms beloning to KCX with LYS and HETATM with ATOM.

        for f in list_of_files:

            try:
                path = "curate_PDB%sconformations%s%s"%(os.sep, os.sep, f)

                with fileinput.FileInput(path, inplace = True) as f:
                    for line in f:
                        if("KCX" in line):
                            line = line.replace("HETATM","ATOM  ")
                            line = line.replace('KCX', 'LYS')
                            print(line, end ='') 
                        else:
                            print(line, end ='') 

            except:
                return


#This script removes any hydrogens from the files in curate_PDB/conformations
    def remove_hydrogens(pdb):
        print('REMOVING HYDROGENS......')
#Firstly it puts each file name that belongs to the pdb of interest into a list
        files = np.array(glob.glob("curate_PDB%sconformations%s*pdb"%(os.sep, os.sep)))
        list_of_files = []
        for file in files:
            file_name = file[25:]
            if file_name[:4] == pdb:
                list_of_files.append(file_name)
#Next it opens them in biobx and gets the index of each non-hydrogen atom
        for file in list_of_files:
            
            path = "curate_PDB%sconformations%s%s"%(os.sep, os.sep, file)
            M = bb.Molecule()
            M.import_pdb(path)
            df = M.data
            list_of_names = df['name'].to_list()
            clean_names = list()

            for name in list_of_names:
                if name[0] != 'H':
                    clean_names.append(name)
            pts, idx = M.atomselect('*', '*', clean_names, get_index=True)
#Next it writes a new pdb including all the atoms except the hydrogens.
            try:
                M.write_pdb(path, index=idx, split_struc=True)
            except:
                M.write_pdb(path, index=idx, split_struc=False)

        return

#This script cleans every file (i.e. makes sure it can be read by biobox and makes sure the N and C count are the same).

    #Fistly it checks if the pdb is in the log file and does ligand check (if enabled).
    keep = checks(pdb, done_pdbs)
    try:

        if keep == True:
            try:
        #Next it cleans the structure and splits into all alternative conformations.
                clean(pdb)
                if pdb == 'FAIL':
                    print('Error cleaning ' + pdb)
                    f = open('log_file.csv', 'a')
                    f.write('\n')
                    f.write(pdb + ', Failed cleaning process')
                    pass
                else:
                    split_struc_NMR(pdb)
                    split_struc_alt_aa(pdb)
                    remove_kcx(pdb)
                    remove_hydrogens(pdb)
            except Exception as e:
                print('Error cleaning ' + pdb)
                f = open('log_file.csv', 'a')
                f.write('\n')
                f.write(pdb + ', Failed cleaning process')


    except Exception as e:
        raise Exception('Error %s'%e)

    return keep


#This module renames the protein's chains durin the cleaning process so that they match the chain names given in the fasta file.
#This is required as pdb files name their chains using the 'auth' name and fasta with the normal chain name.
#Therefore to avoid confusion we rename them all to what is used in the fasta file.

def rename_chains(pdb_code):

    def get_chain_replacement(pdb_code):
        #Fistly the fasta file is downloaded
        try:
            web_url = "https://www.rcsb.org/fasta/entry/" + pdb_code + '/download'

            name = pdb_code + '.fasta'

            if sys.platform == "win32":
                line = "curl -o " + name + " " + web_url
            else:
                line = "wget -O " + name + " " + web_url
                
            subprocess.check_call(line, shell=True)

        except Exception as e:
            print('Error %s'%e)
            print('Error downloading FASTA sequence for chain name comparison.')
            return {}

        try:
            #Next it is opened and all the chain information is appended into the chains_raw list

            f = open(name, 'r')
            list_of_chains = []
            chains_raw =[]
            replacement = []
            need_replacing = []
            replacement_dict = dict()
            for line in f:
                m = re.findall('Chain[ a-z , A-Z \[\]]*|', line)
                for entry in m:
                    if (len(entry) > 0):
                        chains_raw.append(entry)

            #The fasta file is then removed
            f.close()
            os.remove(name)
            
        except Exception as e:
            print('Error %s'%e)
            print('Error parsing through fasta file.')
            f.close()
            os.remove(name)
            return {}

        
#If a chain has two different names (i.e. an auth name and the name given by the RCSB) the chain name will have the following format in the fasta file:
#Chains K[auth M], L[auth N]
#Therefore the code here puts the auth name and name given by RCSB into a dictionary which is later used to replace the auth names.
#e.g. the two examples here would be appended into the dictionary as {M:K, N:L}
#If the auth name and RCSB name are the same nothing is appended to the dictionary.

        for entry in chains_raw:
            try:
                chains_raw_2 = entry.split(',')
                for entry in chains_raw_2:
                    if re.search('auth', entry):
                        chains_wrong = entry.split('auth')
                        for entry in chains_wrong:
                            select = entry.split(' ')
                            for entry in select:
                                if entry[-1:] == '[':
                                    replacement = entry[:-1]
                                    
                                if entry[-1:] == ']':
                                    need_replacing = entry[:-1]
                        replacement_dict[need_replacing] = replacement
                        
                    else:
                        chain = entry.split(' ')
                        for entry in chain:
                            if 'Chain' not in entry and (len(entry) > 0):
                                list_of_chains.append(entry)
            except Exception as e:
                print('Error %s'%e)
                print('Failure renaming chains.')
                continue
            
#This dictionary is then fed onto the next function (replace_chains) where it is used to replace the auth chain names with the RCSB chain names for files in curate_PDB/conformations.

        return replacement_dict




    def replace_chains(pdb_code, replacement_dict):

#Firstly, the auth chain names are put in a list.

        auth_list = list(replacement_dict.keys())

        path = "curate_PDB%sconformations%s%s.pdb"%(os.sep, os.sep, pdb_code)
        

#The relevant file in conformations is then opened and rewritten.

        try:
            with fileinput.FileInput(path, inplace = True) as f:
                for line in f:
                    try:
#On lines with 'ATOM', 'TER' or 'HETATM' if the chain is in the auth_list the auth chain name is replaced with the RCSB chain name.

                        if (line[:4] == 'ATOM') or (line[:3] == 'TER') or (line[:6] == 'HETATM'):
                                chain_name = line[21]
                                if chain_name in auth_list:
                                    replacement_chain_name = replacement_dict.get(chain_name)
                                    line = line[:21] + replacement_chain_name + line[22:]
                                    print(line, end ='')
                                else:
                                    print(line, end = '')

                        else:
                            print(line, end='')
                    except:
                        print(line, end='')

#If the protein fails the replacement process it is removed from the conformations folder.

        except Exception as e:
            print('Error %s'%e)
            print('Error renaming chains')
            os.remove(path)
            pass
        
        return


#The both functions are called within a larger function (rename_chains)
#Here, if there are no chains to be replaced the replace_chains function won't be called.

    replacement_dict = get_chain_replacement(pdb_code)
    if len(replacement_dict) == 0:
        return
    else:
        replace_chains(pdb_code, replacement_dict)

    return


if __name__ == "__main__":

    #try:
    pdb = '1ci4'
    done_pdbs = []
    clean_and_split_alt_conformations(pdb, done_pdbs)
        
    #except Exception as e:
    #    print("ERROR: %s"%e)