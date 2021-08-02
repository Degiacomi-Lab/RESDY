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

def calculate_pKa_and_SASA(code , uniprot_code, AF_lysines_df):
    try:
        pdb_code = code + '.pdb'
        print('Obtaining pKa data for ' + pdb_code)

        process = subprocess.Popen(['python', '-m', 'propka', pdb_code],
                            stdout=subprocess.PIPE, 
                            stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()

        lys_number = list()
        pkas = list()
        pkafile = code + '.pka'
        propres = open(pkafile)
        chain = list()

        for line in propres:
            if re.search('^   LYS' , line):

                line = line[6:]
                line = line.split()
                lys_number.append(line[0])
                chain.append(line[1])
                pkas.append(line[2])

        df = pd.DataFrame({'Uniprot Code': uniprot_code, 'resid':lys_number, 'chain':chain, 'pKa':pkas, 'PDB Code':code})
        df.sort_values(by=['pKa'], inplace=True)
        convert_dict = {'pKa': float}
        df = df.astype(convert_dict)
        convert_dict = {'chain': float}
        df = df.where(df['pKa'] < 20)
        df = df.dropna()
        df = df.drop_duplicates(subset=None, keep='first', inplace=False, ignore_index=True)

        if (len(AF_lysines_df) != 0):
            df_merge = pd.merge(AF_lysines_df, df, how='outer', on='resid')
            df_merge.drop('chain_y', inplace=True, axis=1)
            df_merge.rename(columns={'chain_x': 'chain'}, inplace=True)
            df = df_merge
            df['plddt'] = pd.to_numeric(df['plddt'],errors='coerce')
            df = df.where(df['plddt'] > 70)
            df = df.dropna()

        if df.empty:
            print('No lysine residues with an epsilon amino group pKa < 9 were found in ' + pdb_code + '\n')

        else:
            print('The following lysine residues have epsilon amino group pKa value(s) < 9' + '\n')
            print(df)
            M = bb.Molecule()
            M.import_pdb(pdb_code)
            df_2 = M.data

            resid_list = df['resid'].tolist()
            chain_list = df['chain'].tolist()

            df.index = pd.RangeIndex(len(df.index))
            df.index = range(len(df.index))
            print('\n Obtaining sasa data for ' + pdb_code + '\n')
            list_of_sasa = list()

            for i in range(len(resid_list)):
                percentage_through = (i/(len(resid_list)))*100
                percentage_through = (str(percentage_through))[:4]
                print(percentage_through + '%')


                chain = chain_list[i]
                ResidueID = resid_list[i]
                ResidueID = [int(ResidueID)]
                pts, indices = M.atomselect(chain,  ResidueID, ["CB", "CG", "CD", "CE", "NZ"],  use_resname=False, get_index=True)
                x = sasa(M, targets=indices, probe=1.4, n_sphere_point=960, threshold=0.05)
                acc_surf_area = str(x[0])
                list_of_sasa.append(acc_surf_area)


            df['sasa'] = list_of_sasa
            print(df)
                #s = bb.Structure(p = x[0])
                #s.write_pdb('teststructure.pdb')

    except Exception as e:
        print("ERROR: %s"%e)
        print(pdb_code, 'Failed testing')

    return(df)

if __name__ == "__main__":
    try:    
        print(calculate_pKa_and_SASA('3bg3'))
    except:
        print('Error')