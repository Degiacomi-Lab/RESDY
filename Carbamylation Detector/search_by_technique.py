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
from pka_sasa import calculate_pKa_and_SASA

def search_by_technique(list_of_techniques, list_of_pdbs, wanted_res):

    list_of_pdbs_technique = list()
    wanted_res = float(wanted_res)


    for pdb in list_of_pdbs:
        try:
            print('Checking technique used to produce ' + pdb)
            web_url = 'https://www.rcsb.org/structure/' + pdb
            html = urllib.request.urlopen(web_url)
            soup = BeautifulSoup(html, 'html.parser')

            for line in soup:
                line = str(line)
                messy_resolution = re.findall('Resolution:\s</strong>[\w ()\d - . ]*', line)
                messy_techniques = re.findall('Method:\s</strong>[\w ()\d -]*</li>', line)

        except Exception as e:
            print("ERROR: %s"%e)
            print('Failed to obtain technique/resolution data for ' + pdb)
            continue

        try:
            technique = (messy_techniques[0])[17:-5]

        except Exception as e:
            print("ERROR: %s"%e)
            print('Failed to obtain technique ' + pdb + ' was obtained by.')
            continue

        try:
            resolution = (messy_resolution[0])[21:-1]
            resolution = float(resolution)

        except Exception as e:
            print("ERROR: %s"%e)
            print('Failed to obtain resolution for ' + pdb)
            continue

        try:
            if technique in list_of_techniques and resolution <= wanted_res:
                list_of_pdbs_technique.append(pdb)
                print(pdb + ' added')

        except Exception as e:
            print("ERROR: %s"%e)
            continue

        else:
            continue

    return(list_of_pdbs_technique)



if __name__ == "__main__":
    try:    
        print(search_by_technique(['X-RAY DIFFRACTION'], ['6lvn', '5h7a'], '2.80'))
    except Exception as e:
                print("ERROR: %s"%e)
    
    