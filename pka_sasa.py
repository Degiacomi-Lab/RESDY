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

def calculate_pKa_and_SASA(code):
    #Need to feed list of chains into here...
    try:
        code_for_df = code[:4]
        path = 'assembled/' + code  
        #pdb_code = code + '.pdb'
        print('Obtaining pKa data for ' + code)
        process = subprocess.Popen(['python', '-m', 'propka', path],
                            stdout=subprocess.PIPE, 
                            stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()

        lys_number = list()
        pkas = list()
        pkafile = code_for_df + '_assembled.pka'
        propres = open(pkafile)
        chain = list()

        for line in propres:
            if re.search('^   LYS' , line):

                line = line[6:]
                line = line.split()
                lys_number.append(line[0])
                chain.append(line[1])
                pkas.append(line[2])

        df = pd.DataFrame({'resid':lys_number, 'chain':chain, 'pKa':pkas, 'PDB Code':code_for_df})
        df.sort_values(by=['pKa'], inplace=True)
        convert_dict = {'pKa': float}
        df = df.astype(convert_dict)
        convert_dict = {'chain': float}
        df = df.where(df['pKa'] < 20)
        df = df.dropna()
        df = df.drop_duplicates(subset=None, keep='first', inplace=False, ignore_index=True)


        if df.empty:
            print('No lysine residues with an epsilon amino group pKa < 9 were found in ' + code + '\n')

        else:
            
            print('The following lysine residues have epsilon amino group pKa value(s) < 9' + '\n')
            print(df)
            M = bb.Molecule()
            M.import_pdb(path)
            df_2 = M.data

            resid_list = df['resid'].tolist()
            chain_list = df['chain'].tolist()

            df.index = pd.RangeIndex(len(df.index))
            df.index = range(len(df.index))
            print('\n Obtaining sasa data for ' + code + '\n')
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
        print(code, 'Failed testing')

    return(df)

if __name__ == "__main__":
    try:    
        print(calculate_pKa_and_SASA('3bg3', 'A','P11498'))
    except Exception as e:
        print("ERROR: %s"%e)