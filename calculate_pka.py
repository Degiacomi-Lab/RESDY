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

def calculate_pKa(code):
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
    except Exception as e:
        print('Failed to obtain pKa data for ' + code)

    
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
    lys_number = list()
    pkas = list()
    chain = list()
    chain_resid_list = list()
    buried_percentage = list()
    for line in propres:
        try:
            if re.search('^   LYS' , line):
                print(line)
                line = line[6:]
                line = line.split()
                
                chain_resid_list.append(line[1] + line[0])
                lys_number.append(line[0])
                chain.append(line[1])
                pkas.append(line[2])
                buried_percentage.append(line[3])
        except Exception as e:
            print("Error %s"%e)
            print('Failure parsing ' + code + '.pka')

    if AF_struc == True:
        code_for_df = code[:-4]
    elif AF_struc == False:
        code_for_df = code_for_df
    df = pd.DataFrame({'resid':lys_number, 'chain':chain, 'pKa':pkas, 'PDB Code':code_for_df, 'Chain_Resid': chain_resid_list})
    df.sort_values(by=['pKa'], inplace=True)
    convert_dict = {'pKa': float}
    df = df.astype(convert_dict)
    convert_dict = {'chain': float}
    df = df.where(df['pKa'] < 20)
    df = df.dropna()
    df = df.drop_duplicates(subset=None, keep='first', inplace=False, ignore_index=True)

    return(df)







if __name__ == "__main__":
    try:    
        print(calculate_pKa('4XBJ_assembled.pdb'))
    except Exception as e:
        print("ERROR: %s"%e)