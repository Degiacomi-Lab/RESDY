from os import remove
from re import X
import biobox as bb
import pandas as pd
import fileinput
import re
import os
import subprocess

def clean_metalloprotein(pdb, chain):
    outcome = 'FAILURE'
    try:
        path = 'curate_PDB/raw/' + pdb + '_' + chain + '.pdb'
        f = open(path, 'r')
        file = f.read()
        list_of_chains_no_dup = []
        list_of_chains = re.findall('\w ([A-Z]\S*) \d', file)
        for chain_name in list_of_chains:
            if (len(chain_name) == 1):
                if chain_name not in list_of_chains_no_dup:
                    list_of_chains_no_dup.append(chain_name)
        list_of_metals = ['ZN', 'K', 'NA', 'CA', 'CO', 'MG', 'MN', 'FE', 'CU', 'NI']
    except Exception as e:
        print('Error %s'%e)
        return

    for metal in list_of_metals:
        try:
            for chain_name in list_of_chains_no_dup:
                print(chain_name)
                word_to_be_replaced = metal + '  ' + chain_name + ' ' + metal
                replacement = metal + '   ' + metal
                word_to_be_replaced_2 = '   ' + metal + '  ' + metal
                replacement_2 = '  ' + metal + '  ' + metal
                with fileinput.FileInput(path, inplace = True) as f:
                    for line in f:

                        if(word_to_be_replaced in line):
                            line = line.replace(word_to_be_replaced, replacement)
                            print(line, end='')

                        else:
                            print(line, end ='')

            outcome = 'SUCCESS'
        except Exception as e:
            print('Error %s'%e)
            continue
    return(outcome)







def clean(f):
    try:
        print('CLEANING ')
        #Need to change:
        #remove_kcx(pdb_code, chain)

        M = bb.Molecule()
        M.import_pdb(f, include_hetatm=True)
        print('success')

        pos, idx = M.query('name != ["H"]', get_index=True)


        M.write_pdb("tmp", index=idx)

    


        # call a shell cleaning script (removes hydrogens and alternate side chain conformations)
        # saves a cleaned temporary file called "tmp"


        Mtmp = bb.Molecule()
        Mtmp.import_pdb("tmp", include_hetatm=True)

        if len(Mtmp.coordinates) > 1:
                Mtmp.coordinates = Mtmp.coordinates[0:1]
                Mtmp.set_current(0)

        # remove amino acids with resid < 1
        _, idx = Mtmp.query("resid > 0", get_index=True)
        Mtmp = Mtmp.get_subset(idx)


        M.write_pdb("tmp", index=idx)


        #extract only protein atoms (no water, ligands, DNA, ...)

        list_of_metals = ['ZN', 'NI', 'CU', 'FE', 'MG', 'MN', 'NA', 'K', 'CA', 'CO']
        file_write = open(f, 'w')

        file = open('tmp')
        for line in file:
            if line[:4] == 'ATOM':
                file_write.write(line)
                continue
            if line[:6] == 'HETATM':
                words = line.split()
                print(words)
                if words[3] in list_of_metals:
                    file_write.write(line)
                else:
                    continue
            else:
                file_write.write(line)


        os.remove("tmp")
        print(clean)
    except Exception as e:
        print('Error %s'%e)

    return 


#This code removes the carbamate from any lysine that is already carbamylated in the .pdb structure.

def remove_kcx(pdb_code, chain):
    try:
    #Firstly it imports the pdb file from raw.
        prot_code = pdb_code + '_' + chain + '.pdb'
        path = 'curate_PDB/raw/' + prot_code
        M = bb.Molecule()
        f = open(path, 'r')
        with fileinput.FileInput(path, inplace = True) as f:
            for line in f:
                if("A ZN" in line):
                    line = line.replace(' A ZN', ' ZN')
                    print(line, end ='')
                else:
                    print(line, end ='')

        M.import_pdb(path, include_hetatm=True)

        #Next it removes the carbon and two oxygens from the structure

        idxs = []
        pos ,idx = M.atomignore('*', 'KCX', 'CX', get_index=True, use_resname=True)
        A = M.get_subset(idxs=idx)
        pos, idx = A.atomignore('*', 'KCX', 'OQ1', get_index=True, use_resname=True)
        B = A.get_subset(idxs=idx)
        pos, idx = B.atomignore('*', 'KCX', 'OQ2', get_index=True, use_resname=True)
        C = B.get_subset(idxs=idx)

        #It then rewrites the structure without the carbamate

        C.write_pdb(path, split_struc=False)
        filename = path

        #Lastly it renames the residue as a lysine and relabels it as an atom instead of a hetereoatom.

        with fileinput.FileInput(filename, inplace = True) as f:
            for line in f:
                if("KCX" in line):
                    line = line.replace("HETATM","ATOM  ")
                    line = line.replace('KCX', 'LYS')
                    print(line, end ='') 
                else:
                    print(line, end ='') 
    except:
        return

    return









def remove_kcx_from_fasta(pdb_code, chain):
    path = 'curate_PDB/raw/' + pdb_code + '_' + chain + '.fasta'
    with fileinput.FileInput(path, inplace = True) as f:
        for line in f:
            if("KCX" in line):
                line = line.replace('(KCX)', 'K')
                print(line, end ='')
            else:
                print(line, end ='')
    return()

if __name__ == "__main__":
    try:
        print(clean('4E3T', 'A'))
    except Exception as e:
        print("ERROR: %s"%e)


#def clean_hetatm()