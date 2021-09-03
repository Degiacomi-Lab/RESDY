import biobox as bb
import pandas as pd
import os
import fileinput
import numpy as np
import glob
import re
import subprocess
import math
#get_data.py and the autopatcher break the protein up into chains.
#This code 'reassembles' the protein into a multimer from the chains in the clean folder given the pdb code and chains the protein consists of.


def assemble_multimer(gap_dict, pdb_code, list_chains):
    #First makes sure there is an assembled folder.
    #rc.rename_chains(pdb_code)
    if not os.path.exists("assembled"):
        os.mkdir("assembled")

    #Next open defines the name of the assembly.
    try:
        name_of_assembly = pdb_code + '_assembled.pdb'
        Multi = bb.Multimer()

    #Next opens the pdb file for each chain in turn and appends to Multi
        for chain in list_chains:
            name_of_pdb_file = pdb_code + '_' + chain + '.pdb'
            patched_pdb_file = pdb_code + '_' + chain + '_' + 'patched.pdb'
            M = bb.Molecule()
            try:
                path = 'curate_PDB/clean/' + patched_pdb_file
                M.import_pdb(path, include_hetatm=True)
                Multi.append(M, chain)

            except:
                try:
                    path = 'curate_PDB/clean/' + name_of_pdb_file
                    M.import_pdb(path, include_hetatm=True)
                    Multi.append(M, chain)

                except Exception as e:
                    print("Error: %s"%e)
                    continue

    #Lastly writes out Multi as a .pdb file.

        path = 'assembled/' + name_of_assembly
        Multi.write_pdb(path)
        print('Success assembling ' + pdb_code)
        filename = path
        list_to_remove = check_missing_chains(gap_dict, pdb_code, list_chains)

        with fileinput.FileInput(filename, inplace = True) as f:
            for line in f:
                line = line.replace("TER","")
                print(line, end ='') 

    except Exception as e:
        print("Error: %s"%e)
        print('Failed to assemble ' + pdb_code)

    return list_to_remove


#These next two modules (aswell as one in rename_chains) sort out issues that arise when a chain fails the autopatcher.
#In particular they do two things:

    #Make sure the chain naming is correct.

    #Mark any points that are close to the missing chains to be removed.

def check_missing_chains(gap_dict, pdb_code, list_chains):

    clean_chains = []
    failed_chains =[]

    #Firstly the chains that passed the autopatch (i.e. those in the clean folder) are appended to a list.

    for f in glob.glob("curate_PDB/clean/*.pdb"):

        try:

            if f[-12:] == '_patched.pdb':
                clean_chain = f[-13]
            else:
                clean_chain = f[-5]
            clean_chains.append(clean_chain)

        except Exception as e:
            print('Error %s'%e)
            print('Issue identifying missing chains.')
            continue

    #If any chain is not in clean but is in the original protein it is appended to another list (failed_chains)

    for chain in list_chains:
        if chain not in clean_chains:
            failed_chains.append(chain)


    #If any chains have failed, rename_chains_assembled and note_residues_near_missing _chain sort out chain naming issues and remove residues close to the missing chains respectively.

    if len(failed_chains) != 0:
        print('Chain issue- fixing....')
        try:
            rename_chains_assembled(failed_chains, list_chains)
            list_to_remove = note_residues_near_missing_chain(gap_dict, pdb_code, failed_chains)
        except Exception as e:
            print('Error %s'%e)
            print('Error dealing with chain failure- chain labelling may be incorrect in data.')
            return()
        print('Chain issue fixed.')
    
    else:
        print('No missing chains found')
        list_to_remove = list()

    return list_to_remove



#This program works to rename the chains in an assembled structure.
#This is needed as if any of the chains fail the autopatch, the naming of any chains further down the alphabet shifts.
#For example if the original chains were ABC and B failed then in the assembled the naming would be AB whereas this program renames B to C.

def rename_chains_assembled(failed_chains, list_chains):

#Firstly it creates a copy of the list of chains (a list of all the chains including those which have failed).


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

    files = list((glob.glob("assembled/*pdb")))

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

            print('Error %s'%e)
        
            
    return

















def note_residues_near_missing_chain(gap_dict, pdb_code, failed_chains):
    #Works as follows:
        #Opens pdb file in curate_PDB/conformations/ with Biobox
        #Check each lysine to see if they are near failed chain.
        #If they are add chain_resid to list (list_to_remove)
        #Later in main those in list_to_remove are removed from the results.

  
    path = 'curate_PDB/conformations/' + pdb_code + '.pdb'
    list_of_chains = list()
    list_of_resid = list()
    list_to_remove = []
    
    print('BREAKING UP MOLECULE')
    print(pdb_code)
    try:

    #First opens file in biobox

        M = bb.Molecule()
        M.import_pdb(path, include_hetatm=True)
        df = M.data

    #Next makes note of all the lysine index/coordinates

        lys_coords, lys_idx = M.atomselect('*','LYS', 'CA', use_resname=True, get_index=True)

#Next makes a note of each lysines chain and resid.

        for entry in lys_idx:
            chain = df.at[entry, 'chain']
            list_of_chains.append(chain)

        for entry in lys_idx:
            resid = df.at[entry, 'resid']
            list_of_resid.append(resid)

#Next gets coordinates and index for all atoms

        all_coords, idx = M.atomselect('*','*','*', get_index=True)

    except Exception as e:
        print('Error %s'%e)
        print('Failure identifying residues near missing chain')
        return()

#Next goes through all the lysines
    for failedchain in failed_chains:
        gapsize = gap_dict.get(failedchain)
        gapsize = int(gapsize)
        cutoff = 0.5*gapsize + 3.5

        for j in range(len(lys_coords)):

#Doesn't bother if they are on a chain which failed the autopatch as they won't be in the final structure anyway.

            chain_lys = list_of_chains[j]
            if chain_lys == failed_chains:
                continue
            else:

                try:

#Works out the distance between all the atoms and each lysine.
                    
                    for i in range(len(all_coords)):
                        x_dist = (((lys_coords[j])[0] - (all_coords[i])[0])**2)
                        y_dist = (((lys_coords[j])[1] - (all_coords[i])[1])**2)
                        z_dist = (((lys_coords[j])[2] - (all_coords[i])[2])**2)
                        distance = math.sqrt(x_dist + y_dist + z_dist)

#If the distance is less than 10 A AND the atoms chain is one of the failed chains, the lysine's Chain_Resid is appended to the list_to_remove.
                        if distance < cutoff:
                            index_of_aa = idx[i]
                            aa_chain = df.at[index_of_aa, 'chain']

                            if aa_chain in failed_chains:
                                chain_resid = str(list_of_chains[j]) + str(list_of_resid[j])
                                list_to_remove.append(chain_resid)


                except Exception as e:
                    print('Error %s'%e)
        
    return list_to_remove
        




if __name__ == "__main__":
    try:
        pdb_code = '7K5X'
        list_chains = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M']
        check_missing_chains(pdb_code, list_chains)
    except Exception as e:
        print("ERROR: %s"%e)