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

def calculate_pKa_and_SASA(pdb):
    try:
        pdb_code = pdb + '.pdb'

        process = subprocess.Popen(['python', '-m', 'propka', pdb_code],
                            stdout=subprocess.PIPE, 
                            stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()

        lysines = list()
        pkas = list()
        pkafile = pdb + '.pka'
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
        df = df.where(df['pKa'] < 9)
        df = df.dropna()
        if df.empty:
            print('No lysine residues with an epsilon amino group pKa < 9 were found in ' + pdb + '\n')

        else:
            print('The following lysine residues have epsilon amino group pKa value(s) < 9' + '\n')
            print(df)
            M = bb.Molecule()
            M.import_pdb(pdb_code)
            df_2 = M.data


            resid_list = df['resid'].tolist()


            for ResidueID in resid_list:
                ResidueID = str(ResidueID)
                pts, indices = M.atomselect("*",  ResidueID, ["CA"],  use_resname=False, get_index=True)
                x = sasa(M, targets=indices, probe=1.4, n_sphere_point=960, threshold=0.05)
                print(x)

    except:
        print(pdb_code, 'Failed testing')
    return()

if __name__ == "__main__":
    try:    
        calculate_pKa_and_SASA('1htx')
    except:
        print('Error')