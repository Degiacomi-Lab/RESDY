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
                results_df = get_pdbs(name_of_organism, code)

            except:
                print('Try again')
                skip = 1
            if len(list_of_pdbs) == 0:
                skip = 1



            question_2 = [
            inquirer.List('Choice',
                            message="Do you want to select structures only obtained by certain methods?",
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
                wanted_res = input('What maximum resolution would you like? (In Angstroms)')
            
            elif answer_3['Choice'] == 'No':
                wanted_res = 1000

            if answer_2['Choice'] == 'Yes':

                a = True

                while a == True:
                    
                    question_3 = [
                    inquirer.List('Choice',
                                message="Which methods would you like your structures to have been obtained by?",
                                choices=['Done', 'X-ray', 'NMR', 'EM', 'Fiber', 'IR', 'MODEL', 'Neutron', 'Predicted'],
                            ),
                ]
                    answer_3 = inquirer.prompt(question_3)
                    ans = str(answer_3['Choice'])
                    
                    if ans != 'Done':
                        list_of_techniques.append(ans)

                    elif answer_3['Choice'] == 'Done':
                        a = False

                        if len(list_of_techniques) == 0:
                            list_of_techniques = ['Done', 'X-ray', 'NMR', 'EM', 'Fiber', 'IR', 'MODEL', 'Neutron', 'Predicted']
                
                raw_res = search_by_technique(list_of_techniques, results_df, wanted_res)
                results_df = search_by_technique[0]
                skip = search_by_technique[1]


            elif answer_2['Choice'] == 'No':
                continue
            

        elif answers ['Choice'] == 'Single Protein':

            question_4 = [
            inquirer.List('Choice',
                            message="Do you want to use a PDB or AlphaFold structure?",
                            choices=['PDB', 'AlphaFold'],
                        ),
            ]
            answer_4 = inquirer.prompt(question_4)

            if answer_4["Choice"] == 'PDB':
                try:
                    PDBCODE = input('What is the PDB code?')
                    d = {'Uniprot Entry': 'N/A', 'PDB Code': PDBCODE}
                    columns = ['Uniprot Entry', 'PDB Code']
                    results_df = pd.DataFrame(columns=columns)
                    results_df = results_df.append(d, ignore_index=True)
                    print(results_df)

                except Exception as e:
                    print("ERROR: %s"%e)

            elif answer_4['Choice'] == 'AlphaFold':
                try:
                    AF_code = input('What is the Uniprot code?')
                    AF_code = 'AF-' +AF_code + '-F1-model_v1'
                    d = {'Uniprot Entry': 'N/A', 'PDB Code': AF_code}
                    columns = ['Uniprot Entry', 'PDB Code']
                    results_df = pd.DataFrame(columns=columns)
                    results_df = results_df.append(d, ignore_index=True)
                    print(results_df)

                except Exception as e:
                    print("ERROR: %s"%e)


        elif answers ['Choice'] == 'Quit':
            skip = 1
            raise Exception

        number_of_pdbs = str(len(list_of_pdbs))




        if skip == 0:
            num = 0
            while num < len(results_df):
                code = results_df.at[num, 'PDB Code']
                if len(code) == 4:
                    try:
                        subprocess.check_call("wget https://files.rcsb.org/download/" + code + ".pdb", shell=True)
                    except Exception as e:
                        print("ERROR: %s"%e)
                else:
                    try:
                        subprocess.check_call("wget https://alphafold.ebi.ac.uk/files/" + code + ".pdb", shell=True)
                    except Exception as e:
                        print("ERROR: %s"%e)
                num = num + 1
                calculate_pKa_and_SASA(code)

    except KeyboardInterrupt:
            pass

    except:
        if skip == 1:
            print('\n' + 'Goodbye' + '\n')
            running = False