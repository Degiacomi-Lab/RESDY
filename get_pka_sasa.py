import inquirer
import urllib.request, urllib.parse, urllib.error
import re
from bs4 import BeautifulSoup
import csv
import pandas as pd
import os
import subprocess
import numpy as np
import biobox as bb
from biobox.measures.calculators import sasa
import math


#This function calls propka to calculate the pKa of all groups in the assembled/*.pdb file and then parses the .pka file for the pka of lysine side chains in the protein.

def calculate_pKa(code):

    try:
        code_for_df = code[:4]
        path = 'assembled/' + code  
        print('Obtaining pKa data for ' + code)
        process = subprocess.Popen(['python', '-m', 'propka', path],
                            stdout=subprocess.PIPE, 
                            stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()
    except Exception as e:
        print('Failed to obtain pKa data for ' + code)
        return()

    
    try:
        pkafile = code_for_df + '_assembled.pka'
        propres = open(pkafile)
        AF_struc = False
    except:
        try:
            code_pka = (code[:-4]) + '.pka'
            pkafile = code_pka
            propres = open(pkafile)
            AF_struc = True
        except Exception as e:
            print("Error: %s" %e)
            print('Failed to find ' + code + '.pka')
            return()

    lys_number = list()
    pkas = list()
    chain = list()
    chain_resid_list = list()
    buried_percentage = list()
    
    try:
        for line in propres:
            if re.search('^   LYS' , line):
                try:
                    print(line)
                    line = line[6:]
                    line = line.split()
                    
                    chain_resid_list.append(line[1] + line[0])
                    lys_number.append(line[0])
                    chain.append(line[1])
                    pkas.append(line[2])
                except Exception as e:
                    print("Error %s"%e)
                    continue
    except Exception as e:
        print("Error %s"%e)
        print('Failure parsing ' + code + '.pka')
        return()

    if AF_struc == True:
        code_for_df = code[:-4]
    elif AF_struc == False:
        code_for_df = code_for_df


    try:
        df = pd.DataFrame({'resid':lys_number, 'chain':chain, 'pKa':pkas, 'PDB Code':code_for_df, 'Chain_Resid': chain_resid_list})
        df.sort_values(by=['pKa'], inplace=True)
        convert_dict = {'pKa': float}
        df = df.astype(convert_dict)
        convert_dict = {'chain': float}
        df = df.where(df['pKa'] < 20)
        df = df.dropna()
        df = df.drop_duplicates(subset=None, keep='first', inplace=False, ignore_index=True)
    except Exception as e:
        print("Error %s"%e)
        print('Failed to construct pKa dataframe')
        return()
    return(df)






#Forms small structures which include just the atoms surrounding the lysine of interest and computes SASA from that.

def break_up_and_calculate_sasa(pdb_code):

    try:
        print('BREAKING UP MOLECULE')
        print(pdb_code)
        list_of_index = list()
        list_of_sasa = list()
        list_of_resid = list()
        list_of_chains = list()
        chain_resid_list = list()

#Accesses .pdb file from assembled/
        M = bb.Molecule()
        path = 'assembled/' + pdb_code
        print(path)
        M.import_pdb(path, include_hetatm=True)
        df = M.data

#Finds the coordinates and index of all lysine residues in the protein.

        lys_coords, lys_idx = M.atomselect('*',['LYS', 'BLYS'], 'NZ', use_resname=True, get_index=True)
        df = M.data
#Use get_subset...
#Finds the chain and resid of each lysine.

        for entry in lys_idx:
            chain = df.at[entry, 'chain']
            list_of_chains.append(chain)

        for entry in lys_idx:
            resid = df.at[entry, 'resid']
            list_of_resid.append(resid)
        print(list_of_chains)

#Finds the coordinates and index of every atom in the molecule.

        all_coords, idx = M.atomselect('*','*','*', get_index=True)

    except Exception as e:
        print("Error %s"%e)
        return()
#For each lysine it works out the distance between the lys NZ and the each atom in the protein.

    for j in range(len(lys_coords)):
        list_close_points = list()
        
        for i in range(len(all_coords)):
            try:
                x_dist = (((lys_coords[j])[0] - (all_coords[i])[0])**2)
                y_dist = (((lys_coords[j])[1] - (all_coords[i])[1])**2)
                z_dist = (((lys_coords[j])[2] - (all_coords[i])[2])**2)
                distance = math.sqrt(x_dist + y_dist + z_dist)
                if distance < 15:
                    list_close_points.append(idx[i])
            except:
                continue

#If the atoms are close to the lys NZ they are included in a small .pdb structure.
        try:
            M.write_pdb('temp_struc.pdb', index=list_close_points, split_struc=False)

            S = bb.Molecule()
            S.import_pdb('temp_struc.pdb', include_hetatm=True)
            chain = list_of_chains[j]
            resid = list_of_resid[j]

#sasa is calculated for that lysine in the small molecule.

            pts_2, indx_2 = S.atomselect(chain, [resid], ["CB", "CG", "CD", "CE", "NZ"],  use_resname=False, get_index=True)
            x = sasa(S, targets=indx_2, probe=1.4, n_sphere_point=960, threshold=0.05)
            print(x[0])
            chain_resid_list.append(str(chain + str(resid)))
            list_of_sasa.append(x[0])

        except:
            print('Error obtaining SASA for ' + pdb_code + ' index value ' + str(j))
            list_of_sasa.append(None)
            continue

#The results are appended to a df which is given as outpit
    try:
        df = pd.DataFrame({'Assembled Index': lys_idx, 'chain': list_of_chains, 'resid': list_of_resid, 'sasa': list_of_sasa, 'Chain_Resid':chain_resid_list})
        print(df)

    except Exception as e:
        print("Error %s"%e)
        print('Error obtaining sasa data for ' + pdb_code)
        return()

    try:
        os.remove('temp_struc.pdb')
    except Exception as e:
        print("Error %s"%e)
        print('Failed to remove temporary pdb structure for ' + pdb_code)
        pass
        

    return(df)







def identity_of_local_amino_acids(pdb_code):
    list_of_chains = list()
    list_of_resid = list()
    columns = ['PDB Code', 'Lysine Index', 'Local Amino Acid', 'Distance']
    local_aa_df = pd.DataFrame(columns=columns)
    print('BREAKING UP MOLECULE')
    print(pdb_code)

    M = bb.Molecule()
    path = 'assembled/' + pdb_code
    print(path)
    M.import_pdb(path, include_hetatm=True)
    df = M.data
    df.to_csv('1dpm_info.csv')
    print(df)
    lys_coords, lys_idx = M.atomselect('*','LYS', 'NZ', use_resname=True, get_index=True)


    for entry in lys_idx:
        chain = df.at[entry, 'chain']
        list_of_chains.append(chain)

    for entry in lys_idx:
        resid = df.at[entry, 'resid']
        list_of_resid.append(resid)

    all_coords, idx = M.atomselect('*','*','*', get_index=True)


    for j in range(len(lys_coords)):
        try:
            list_close_points = list()
            
            for i in range(len(all_coords)):
                x_dist = (((lys_coords[j])[0] - (all_coords[i])[0])**2)
                y_dist = (((lys_coords[j])[1] - (all_coords[i])[1])**2)
                z_dist = (((lys_coords[j])[2] - (all_coords[i])[2])**2)
                distance = math.sqrt(x_dist + y_dist + z_dist)
                if distance < 10:
                    index_of_aa = idx[i]
                    aa_resid = df.at[index_of_aa, 'resid']
                    lys_resid = list_of_resid[j]

                    if aa_resid == lys_resid:
                        continue
                    else:
                        identity = df.at[index_of_aa, 'resname']
                        d = ({'PDB Code': pdb_code, 'Lysine Index': lys_idx[j], 'Local Amino Acid': identity, 'Distance':distance})
                        local_aa_df = local_aa_df.append(d, ignore_index=True)
                    #print(local_aa_df)
        except Exception as e:
            print('Error %s'%e)
    local_aa_df.to_csv('local_aa_near_lysine.csv')
    return(local_aa_df)





if __name__ == "__main__":
    try:    
        print(break_up_and_calculate_sasa('6LVN_assembled.pdb'))
        print(calculate_pKa('4XBJ_assembled.pdb'))
    except Exception as e:
        print("ERROR: %s"%e)