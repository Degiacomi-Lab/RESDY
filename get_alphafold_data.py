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
import subprocess

#This code downloads alphafold structures into the assembled folder.

def download_AF_struc(pdb):
    oldcwd = os.getcwd()
    if not os.path.exists("assembled"):
        os.mkdir("assembled")
    os.chdir('assembled')
    subprocess.check_call("wget https://alphafold.ebi.ac.uk/files/" + pdb + ".pdb", shell=True)
    os.chdir(oldcwd)
    print('PDB Structure for ' + pdb + ' downloaded')
    return()



#Code obtains the plddt (a measure of certainty where 100 is high and 70 low) value for each lysine in an alphafold structure.



def find_AF_plddt(AF_code_full):
    try:
        #Opens .pdb file in assembled folder
        columns = ['resid', 'chain', 'plddt']
        AF_lysines_df = pd.DataFrame(columns=columns)

    except Exception as e:
        print("Error %s"%e)
        return(AF_lysines_df)

    try:
        f = open('assembled/' + AF_code_full, "r")

        #Parses though file to find plddt value.

        for line in f:
            try:
                if re.search('CA  LYS', line):
                    strline = str(line)
                    data = strline.split()
                    if len(data[4]) > 1:
                        resid = (data[4])[1:]
                        plddt = data[9]
                        chain = (data[4])[0]
                    elif len(data[4]) == 1:
                        data = strline.split()
                        resid = data[5]
                        plddt = data[10]
                        chain = data[4]

                    data = ({'resid': resid, 'chain': chain, 'plddt': plddt})

            except Exception as e:
                print("Error %s"%e)
                data = ({'resid': resid, 'chain': chain, 'plddt': '0'})
                continue

        #Appends to dataframe which is later merged into the main dataframe.
                
            AF_lysines_df = AF_lysines_df.append(data, ignore_index=True)
        print(AF_lysines_df)

    except Exception as e:
        print("ERROR: %s"%e)
        print('Failed to obtain pLDDT data for ' + AF_code_full)
        return(AF_lysines_df)

    return(AF_lysines_df)
