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
from search_by_technique import search_by_technique
skip = 0
running = True
list_of_techniques = list()

while running:
    try:
        list_of_pdbs = list()
        print('\n' + '------------------------------------------------------------')
        print('\n' + '  Hello, press Ctrl + C at anytime to return to the start' + '\n')
        print('------------------------------------------------------------' + '\n')
        questions = [
        inquirer.List('Choice',
                        message="Do you want to search a whole organism's proteome or a single protein?",
                        choices=['Whole Proteome', 'Single Protein', 'Quit'],
                    ),
        ]
        answers = inquirer.prompt(questions)

        if answers["Choice"] == 'Whole Proteome':

            name_of_organism = input('Name of organism:')
            name_of_organism = name_of_organism.replace(' ', '+')
            code = input('Uniprot code:')


            try:
                list_of_pdbs = get_pdbs(name_of_organism, code)
            except:
                print('Try again')
                skip = 1
            if len(list_of_pdbs) == 0:
                skip = 1



            question_2 = [
            inquirer.List('Choice',
                            message="Do you want to only select structures obtained by certain methods?",
                            choices=['Yes', 'No'],
                        ),
            ]
            answer_2 = inquirer.prompt(question_2)

            question_3 = [
            inquirer.List('Choice',
                            message="Do you want to select structures based on resolution?",
                            choices=['Yes', 'No'],
                        ),
            ]
            answer_3 = inquirer.prompt(question_3)

            if answer_3['Choice'] == 'Yes':
                wanted_res = input('What maximum resolution would you like?')
            
            elif answer_3['Choice'] == 'No':
                wanted_res = 1000

            if answer_2['Choice'] == 'Yes':

                a = True

                while a == True:
                    
                    question_3 = [
                    inquirer.List('Choice',
                                message="Which methods would you like your structures to have been obtained by?",
                                choices=['Done', 'X-RAY DIFFRACTION', 'SOLUTION NMR', 'ELECTRON MICROSCOPY', 'NEUTRON DIFFRACTION', 'ELECTRON CRYSTALLOGRAPHY', 'SOLID-STATE NMR', 'SOLUTION SCATTERING', 'FIBER DIFFRACTION', 'POWDER DIFFRACTION', 'EPR', 'THEORETICAL MODEL', 'INFRARED SPECTROSCOPY', 'FLUORESCENCE TRANSFER'],
                            ),
                ]
                    answer_3 = inquirer.prompt(question_3)
                    ans = str(answer_3['Choice'])
                    
                    if ans != 'Done':
                        list_of_techniques.append(ans)

                    elif answer_3['Choice'] == 'Done':
                        a = False

                        if len(list_of_techniques) == 0:
                            list_of_techniques = ['Done', 'X-RAY DIFFRACTION', 'SOLUTION NMR', 'ELECTRON MICROSCOPY', 'NEUTRON DIFFRACTION', 'ELECTRON CRYSTALLOGRAPHY', 'SOLID-STATE NMR', 'SOLUTION SCATTERING', 'FIBER DIFFRACTION', 'POWDER DIFFRACTION', 'EPR', 'THEORETICAL MODEL', 'INFRARED SPECTROSCOPY', 'FLUORESCENCE TRANSFER']

                list_of_pdbs = search_by_technique(list_of_techniques, list_of_pdbs, wanted_res)

            elif answer_2 == 'No':
                continue




        elif answers ['Choice'] == 'Single Protein':
            PDBCODE = input('What is the PDB code?')
            list_of_pdbs.append(PDBCODE)

        elif answers ['Choice'] == 'Quit':
            skip = 1
            raise Exception

        number_of_pdbs = str(len(list_of_pdbs))

        #DOWNLOADS PDB FILES

        if skip == 0:
            for pdb in list_of_pdbs:
                subprocess.check_call("wget https://files.rcsb.org/download/" + pdb + ".pdb", shell=True)

            calculate_pKa_and_SASA(pdb)

    except KeyboardInterrupt:
            pass

    except:
        if skip == 1:
            print('\n' + 'Goodbye' + '\n')
            running = False