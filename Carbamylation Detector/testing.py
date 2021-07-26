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


question_2 = [
inquirer.List('Choice',
                message="Do you want to only select structures obtained by certain methods?",
                choices=['Yes', 'No'],
            ),
]
answer_2 = inquirer.prompt(question_2)

list_of_pdbs = ['1d6d', '1eeg', '5h71', '167d', '6t83', '6s47', ]
list_of_pdbs_technique = list()

if answer_2['Choice'] == 'Yes':
    a = True
    list_of_techniques = list()
    while a == True:
        question_3 = [
        inquirer.List('Choice',
                    message="Which methods would you like your structures to have been obtained by?",
                    choices=['Done', 'X-RAY DIFFRACTION', 'SOLUTION NMR', 'ELECTRON MICROSCOPY', 'NEUTRON DIFFRACTION', 'ELECTRON CRYSTALLOGRAPHY', 'SOLID-STATE NMR', 'SOLUTION SCATTERING', 'FIBER DIFFRACTION', 'POWDER DIFFRACTION', 'EPR', 'THEORETICAL MODEL', 'INFRARED SPECTROSCOPY', 'FLUORESCENCE TRANSFER'],
                ),
    ]
        answer_3 = inquirer.prompt(question_3)
        if answer_3['Choice'] != 'Done':
            list_of_techniques.append(answer_3['Choice'])
            
        if answer_3['Choice'] == 'Done':
            a = False
            if len(list_of_techniques) >= 1:
                for pdb in list_of_pdbs:
                    web_url = 'https://www.rcsb.org/structure/' + pdb
                    html = urllib.request.urlopen(web_url)
                    list_UNIPROT_codes = list()
                    soup = BeautifulSoup(html, 'html.parser')

                    for line in soup: 
                        line = str(line)
                        messy_techniques = re.findall('Method:\s</strong>[\w ()\d -]*</li>', line)
                    technique = (messy_techniques[0])[17:-5]

                    if technique in list_of_techniques:
                        list_of_pdbs_technique.append(pdb)
                    else:
                        continue
        print(list_of_pdbs_technique)