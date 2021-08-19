import os
import subprocess
import re
import numpy as np
import glob
import fileinput
import biobox as bb
import pandas as pd
#OVERALL STRUCTURE...
#First downloads file into conformations folder
#Next cleans (i.e. removes heteroatoms that aren't metal ions)
#Next writes file for alt NMR confs
#Next writes file for alt AA confs




#This script downloads the file and removes all the hetereoatoms that aren't metal ions.
#The results is saved into a file with the structure *PDB code*-clean.pdb.
def clean_and_split_alt_conformations(pdb):
    def clean(pdb):

        print('CLEANING')
        try:
            #Firstly it goes into curate_PDB/conformations and downloads the .pdb file.
            cwd = os.getcwd()
            print(cwd)

            os.chdir('curate_PDB/conformations')
            subprocess.check_call("wget https://files.rcsb.org/download/" + pdb + '.pdb', shell=True)

            #Next it returns to the carbamylation folder.

            os.chdir(cwd)
        except Exception as e:
            print('Error %s'%e)
            print('Print failure opening file.')
            return

    #Next it opens and starts reading the .pdb file and starts writing a new file with the ending '-clean.pdb'.

        list_of_metals = ['ZN', 'NI', 'CU', 'FE', 'MG', 'MN', 'NA', 'K', 'CA', 'CO', 'CL', 'MO']

        read_file_path = 'curate_PDB/conformations/' + pdb + '.pdb'

        read_file = open(read_file_path)

        write_file_path = 'curate_PDB/conformations/' + pdb + '-clean.pdb'
        write_file = open(write_file_path, 'w')

    #Next it writes the clean file, including HETATMs (if they are metal ions), all atoms and lines starting with TER and END.

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

        os.remove(read_file_path)
        print('Clean ' + pdb)

        return






    #This writes a new file for each alternate NMR structure.

    def split_struc_NMR(pdb):

    #Firstly it defines the name format which will be followed for each file.

        number = 1
        name = pdb + '-alt-' + str(number) + '.pdb'
        path = 'curate_PDB/conformations/' + pdb + '-clean.pdb'

    #Next it opens the file produced from teh cleaning script and opens a new file to write in.

        f = open(path)

        endmdls = list()

    #Next it searches to see if there are any 'ENDMDL' statements in the folder (i.e. if there are multiple models).


        for line in f:
            if re.search('ENDMDL', line):
                endmdls.append(line)

    #If there are no 'ENDMDL' statements the clean file is renamed to suit the new format.

        if len(endmdls) == 0:
            path_rename = 'curate_PDB/conformations/' + name
            os.rename(path, path_rename)
            print('NO ALTERNATE WHOLE STRUCTURES FOUND FOR ' + pdb)


    #If there are 'ENDMDL' statements a new file is written for each model.

        elif len(endmdls) != 0:
            
            print('ALTERNATE WHOLE STRUCTURE FOUND FOR ' + pdb)

    #First it opens the clean file in the conformations folder and opens a new folder to write in

            f = open(path)
            path_rename = 'curate_PDB/conformations/' + name
            f_write = open(path_rename, 'w')

    #Next it writes the new file.
    #It includes every line until it gets to 'ENDMDL', where it opens a new file to write in.
    #The process stops when it gets to 'MASTER'

            for line in f:
                try:
                    line = str(line)

                    if re.search('END ', line):
                        print(path)
                        os.remove(path)
                        os.remove(path_rename)
                        break

                    else:

                        if re.search('ENDMDL', line):
                            f_write.write(line)
                            f_write.close()
                            number = number + 1
                            name = pdb + '-alt-' + str(number) + '.pdb'
                            path_rename = 'curate_PDB/conformations/' + name
                            f_write = open(path_rename, 'w')

                        else:
                            f_write.write(line)
                except:
                    continue
        
        
        return()




    #Writes a new file for each alternative amino acid conforation present.

    def split_struc_alt_aa(pdb):

#Firstly it gets a list of all the .pdb files present in curate_PDB/conformations and selects those that belong to the pdb we are interested in.

        files = np.array(glob.glob("curate_PDB/conformations/*pdb"))
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
                    path = 'curate_PDB/conformations/' + f
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

                            write_path = 'curate_PDB/conformations/' + name
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
                os.remove(path)
        return()





#This script removes the carbamate from any structure

    def remove_kcx(pdb):

#Firstly it gets a list of files in curate_PDB/conformations that belong to the pdb of interest.

        files = np.array(glob.glob("curate_PDB/conformations/*pdb"))
        list_of_files = []
        for file in files:
            file_name = file[25:]
            if file_name[:4] == pdb:
                list_of_files.append(file_name)
        
        for file in list_of_files:
            path = 'curate_PDB/conformations/' + file
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
                path = 'curate_PDB/conformations/' + f

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
        files = np.array(glob.glob("curate_PDB/conformations/*pdb"))
        list_of_files = []
        for file in files:
            file_name = file[25:]
            if file_name[:4] == pdb:
                list_of_files.append(file_name)
#Next it opens them in biobx and gets the index of each non-hydrogen atom
        for file in list_of_files:
            path = 'curate_PDB/conformations/' + file
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
                M.write_pdb(path, index=idx, split_struc=True),
            except:
                M.write_pdb(path, index=idx, split_struc=False)

        return


        return




#This script cleans every file (i.e. makes sure it can be read by biobox and makes sure the N and C count are the same).



    try:
        clean(pdb)
        split_struc_NMR(pdb)
        split_struc_alt_aa(pdb)
        remove_kcx(pdb)
        remove_hydrogens(pdb)

    except Exception as e:
        print('Error %s'%e)


    return()

        






if __name__ == "__main__":

    try:
        pdb = '4XBJ'
        clean_and_split_alt_conformations(pdb)
        
    except Exception as e:
        print("ERROR: %s"%e)