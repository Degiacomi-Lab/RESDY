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
def find_AF_plddt(code):
    try:
        f = open("AF-O15552-F1-model_v1.pdb", "r")
        columns = ['resid', 'chain', 'plddt']
        AF_lysines_df = pd.DataFrame(columns=columns)
        for line in f:
            if re.search('CA  LYS', line):
                strline = str(line)
                data = strline.split()
                resid = data[5]
                plddt = data[10]
                chain = data[4]
                data = ({'resid': resid, 'chain': chain, 'plddt': plddt})
                AF_lysines_df = AF_lysines_df.append(data, ignore_index=True)
        print(AF_lysines_df)
    except Exception as e:
        print("ERROR: %s"%e)
        print('Failed to obtain pLDDT data for ' + code)

    return(AF_lysines_df)
    #lysine_data = re.findall('name="[\w ]*" type', line)