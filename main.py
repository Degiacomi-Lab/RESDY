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
                print(results_df)
            except:
                print('Try again')
                skip = 1
            if (len(results_df) == 0):
                skip = 1
                #This skip bit doesn't work...

            if skip != 1:
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
                


                if (answer_2['Choice'] == 'Yes') or (answer_3['Choice'] == 'Yes'):
                    if answer_2['Choice'] == 'Yes':
                        a = True
                        list_of_techniques = list()
                        while a == True:

                            question_5 = [
                            inquirer.List('Choice',
                                        message="Which methods would you like your structures to have been obtained by?",
                                        choices=['Done', 'X-ray', 'NMR', 'EM', 'Fiber', 'IR', 'MODEL', 'Neutron', 'Predicted'],
                                    ),
                            ]
                            answer_5 = inquirer.prompt(question_5)
                            ans = str(answer_5['Choice'])


                            if ans != 'Done':
                                list_of_techniques.append(ans)
                            if ans == 'Done':
                                a = False
                                

                    elif answer_2['Choice'] == 'No':
                        list_of_techniques = ['X-ray', 'NMR', 'EM', 'Fiber', 'IR', 'MODEL', 'Neutron', 'Predicted']

                    raw_res = search_by_technique(list_of_techniques, results_df, wanted_res)
                    results_df = raw_res[0]
                    print(results_df)
                    skip = raw_res[1]


                else:
                    pass
            

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

        columns = ['resid', 'chain', 'pKa', 'sasa', 'PDB Code']
        all_pka_sasa_res = pd.DataFrame(columns=columns)
    
        if skip == 0:
            num = 0
            results_df.index = pd.RangeIndex(len(results_df.index))

            results_df.index = range(len(results_df.index))
            while num < len(results_df):
            #while num < 30:

                code = results_df.at[num, 'PDB Code']
                if len(code) == 4:
                    try:
                        subprocess.Popen("ls", cwd="data/")
                        subprocess.check_call("wget https://files.rcsb.org/download/" + code + ".pdb", shell=True)
                    except Exception as e:
                        print("ERROR: %s"%e)
                else:
                    try:
                        subprocess.check_call("wget https://alphafold.ebi.ac.uk/files/" + code + ".pdb", shell=True)
                    except Exception as e:
                        print("ERROR: %s"%e)

                num = num + 1
                print(num)

                try:
                    pka_sasa_df = calculate_pKa_and_SASA(code)
                    all_pka_sasa_res = all_pka_sasa_res.append(pka_sasa_df)
                    print(all_pka_sasa_res)
                    all_pka_sasa_res.to_csv('results.csv')

                except Exception as e:
                    print("ERROR: %s"%e)
                    continue

                try:
                    pathway = '/Users/charliebrown/code/GitHub/carbamylation/Carbamylation_Detector/' + code + '.pdb'
                    args = ('rm', '-rf', pathway)
                    subprocess.call('%s %s %s' % args, shell=True)

                    pathway_2 = '/Users/charliebrown/code/GitHub/carbamylation/Carbamylation_Detector/' + code + '.pka'
                    args_2 = args = ('rm', '-rf', pathway_2)
                    subprocess.call('%s %s %s' % args_2, shell=True)

                except Exception as e:
                    print("ERROR: %s"%e)

                
                

        print(all_pka_sasa_res)
        


    except KeyboardInterrupt:
        pass

    except:
        if skip == 1:
            print('\n' + 'Goodbye' + '\n')
            running = False