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
from uniprot import get_pdbs

def search_by_technique(list_of_techniques, results_df, wanted_res):
    skip = 0
    try:

        results_df = results_df[results_df['Method Structure Obtained by'].isin(list_of_techniques)]

    except Exception as e:
        print("ERROR: %s"%e)
        print('Failed to search by method structure obtained by')
        question = [
        inquirer.List('Choice',
                            message="Continue?",
                            choices=['Yes', 'No'],
                        ),
        ]
        answer = inquirer.prompt(question)

        if answer['Choice'] == 'Yes':
            skip = 0
            pass

        elif answer['Choice'] == 'No':
            skip = 1





    try:
        resolution = float(wanted_res)

        results_df['Resolution'] = results_df['Resolution'].astype(float)

        results_df = results_df[results_df['Resolution']<= resolution]

    except Exception as e:
        print("ERROR: %s"%e)
        print('Failed to search by resolution')
        question = [
        inquirer.List('Choice',
                            message="Continue?",
                            choices=['Yes', 'No'],
                        ),
        ]
        answer = inquirer.prompt(question)

        if answer['Choice'] == 'Yes':
            skip = 0
            pass

        elif answer['Choice'] == 'No':
            skip = 1

    
    return(results_df, skip)



if __name__ == "__main__":

    try:

        columns = ['Uniprot Entry', 'PDB Code', 'Method Structure Obtained by', 'Resolution']
        df = pd.DataFrame(columns=columns)
        uniprot_code = 'P3892'
        PDB_entries = '2jfi'
        technique = 'X-ray'
        resolution = '2.9'


        data = ({'Uniprot Entry': uniprot_code, 'PDB Code': PDB_entries, 'Method Structure Obtained by': technique, 'Resolution':resolution})
        df = df.append(data, ignore_index=True)
        results_df = df

        res_raw = (search_by_technique(['X-ray'], results_df, '3.0'))
        results_df = res_raw[0]
        skip = res_raw[1]
        print(results_df)

    except Exception as e:
                print("ERROR: %s"%e)
    
    