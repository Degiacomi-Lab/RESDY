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
    try:
        pdb_code = code + '.pdb'
        print('Obtaining pKa data for ' + pdb_code)

        process = subprocess.Popen(['python', '-m', 'propka', pdb_code],
                            stdout=subprocess.PIPE, 
                            stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()

        lysines = list()
        pkas = list()
        pkafile = code + '.pka'
        propres = open(pkafile)
        chain = list()

        for line in propres:
            if re.search('^   LYS' , line):
                line = line.strip()

                res = re.findall('LYS(.* [A-Z])', line)
                chain.append((res[0])[-1:])
                res = (res[0])[:-1]

                
                if res[0] == ' ':
                    res = res[1:]

                lysines.append(res)

                pka = line[13:18]
                pkas.append(pka)

        df = pd.DataFrame({'resid':lysines, 'chain':chain, 'pKa':pkas})
        df.sort_values(by=['pKa'], inplace=True)
        convert_dict = {'pKa': float}
        df = df.astype(convert_dict)
        convert_dict = {'chain': float}
        df = df.where(df['pKa'] < 20)
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

            columns = ['resid', 'chain', 'pKa', 'sasa', 'PDB Code']
            pka_sasa_df = pd.DataFrame(columns=columns)
            df.index = pd.RangeIndex(len(df.index))
            df.index = range(len(df.index))
            print('\n Obtaining sasa data for ' + pdb_code + '\n')

            for i in range(len(resid_list)):
                chain = chain_list[i]
                ResidueID = resid_list[i]
                ResidueID = [int(ResidueID)]
                pts, indices = M.atomselect(chain,  ResidueID, ["NZ"],  use_resname=False, get_index=True)
                x = sasa(M, targets=indices, probe=1.4, n_sphere_point=960, threshold=0.05)
                pka_value = df.at[i, 'pKa']
                data = ({'PDB Code': code, 'resid': ResidueID[0], 'chain': chain, 'pKa': pka_value, 'sasa': x[0]})
                pka_sasa_df = pka_sasa_df.append(data, ignore_index=True)
                #s = bb.Structure(p = x[0])
                #s.write_pdb('teststructure.pdb')

    except Exception as e:
        print("ERROR: %s"%e)
        print(pdb_code, 'Failed testing')

    return(pka_sasa_df)

if __name__ == "__main__":
    try:    
        print(calculate_pKa_and_SASA('3bg3'))
    except:
        print('Error')