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
import sys
import glob
from biobox.measures.calculators import sasa
from scipy.spatial.distance import _correlation_pdist_wrap
from uniprot import get_pdbs
from get_pdbs_given_uniprot import get_pdbs_given_uniprot_code
from construct_single_pdb_df import construct_single_pdb_df
from download_AF_struc import download_AF_struc
from calculate_pka import calculate_pKa
from calculate_sasa import calculate_sasa
from search_by_technique import search_by_technique
from alphafold_certainty import find_AF_plddt
from data_analysis import analyse_data
from average_by_prot import average_prot
from curate_PDB.get_data import get_data
from assemble_multimer import assemble_multimer

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
                        message="Do you want to search a organism's whole proteome or input specific PDB or Uniprot Codes?",
                        choices=['Whole Proteome', 'Input PDB Codes', 'Input Uniprot Codes', 'Quit'],
                    ),
        ]
        answers = inquirer.prompt(questions)


        #Searches Uniprot for PDB codes
        if (answers['Choice'] == 'Whole Proteome') or (answers['Choice'] == 'Input Uniprot Codes'):
            if answers["Choice"] == 'Whole Proteome':

                name_of_organism = input('Name of organism:')
                name_of_organism = name_of_organism.replace(' ', '+')
                code = input('Uniprot proteome code:')

                try:
                    results_df = get_pdbs(name_of_organism, code)
                    print(results_df)
                except:
                    print('Try again')
                    skip = 1
                if (len(results_df) == 0):
                    skip = 1
                    
            if answers['Choice'] == 'Input Uniprot Codes':

                asking_for_uniprot_codes = True
                list_of_UNIPROT_codes = list()

                while asking_for_uniprot_codes == True:

                    uniprot_code = input('Uniprot code:')
                    list_of_UNIPROT_codes.append(uniprot_code)

                    question_more_uniprot_codes = [
                    inquirer.List('Choice',
                                    message="Are there more uniprot codes to add?",
                                    choices=['Yes', 'No'],
                                ),
                    ]
                    answer_more_uniprot_codes = inquirer.prompt(question_more_uniprot_codes)

                    if answer_more_uniprot_codes['Choice'] == 'No':
                        asking_for_uniprot_codes = False
                    else:
                        pass
                
                print('Getting data for: ')
                print(list_of_UNIPROT_codes)
                results_df = get_pdbs_given_uniprot_code(list_of_UNIPROT_codes)

                if (len(results_df) == 0):
                    skip = 1




            #Allows user to chose if they only want to select structures obtained by certain techniques/of certain resolution.

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
            


        #Allows single structure to be tested (probably needs to be put in a function)
        elif answers['Choice'] == 'Input PDB Codes':
            select_strucs = True

            columns = ['Uniprot Entry', 'PDB Code', 'Method Structure Obtained by', 'Resolution', 'Chains']
            results_df = pd.DataFrame(columns=columns)

            while select_strucs == True:

                question_4 = [
                inquirer.List('Choice',
                                message="Do you want to use PDB or AlphaFold structure(s)?",
                                choices=['PDB', 'AlphaFold', 'Done'],
                            ),
                ]
                answer_4 = inquirer.prompt(question_4)

                if answer_4["Choice"] == 'PDB':
                    try:
                        PDBCODE_inpt = input('What is the PDB code?')
                        PDBCODE_inpt = PDBCODE_inpt.upper()
                        print(PDBCODE_inpt)
                        UNIPROT_code_pdb = input('What is the Uniprot Code?')
                        results_df = construct_single_pdb_df(UNIPROT_code_pdb, PDBCODE_inpt, results_df)
                        print(results_df)

                    except Exception as e:
                        print("ERROR: %s"%e)
                        print('Unable to obtain data for entry.')
                        pass

                elif answer_4['Choice'] == 'AlphaFold':
                    try:
                        AF_code = input('What is the Uniprot code?')
                        AF_code_full = 'AF-' + AF_code + '-F1-model_v1'
                        d = {'Uniprot Entry': AF_code, 'PDB Code': AF_code_full, 'Method Structure Obtained by': 'Predicted', 'Resolution': 'N/A', 'Chains': 'A'}
                        results_df = results_df.append(d, ignore_index=True)
                        print(results_df)

                    except Exception as e:
                        print("ERROR: %s"%e)
                
                if answer_4['Choice'] == 'Done':
                    select_strucs = False


        elif answers ['Choice'] == 'Quit':
            skip = 1
            raise Exception

        dict_uniprot = dict()
        if (len(results_df) != 0):
            for i in range(len(results_df)):
                pdb = results_df.at[i, 'PDB Code']
                chain = results_df.at[i, 'Chains']
                uniprot = results_df.at[i, 'Uniprot Entry']
                dict_uniprot.update({pdb: uniprot})
                print('*****************GETTING DATA FOR ' + str(pdb) + str(chain) + '*******************')
                print(chain)
                if pdb[:2]=='AF':
                    try:
                        download_AF_ystruc(pdb)
                    except:
                        continue
                else:
                    print('******* Getting data for ' + pdb + ' ' + chain + ' **********')
                    try:
                        get_data(pdb, chain)
                        for f in glob.glob("curate_PDB/raw/*"):
                            os.remove(f)
                    except Exception as e:
                        print('Error: %s'%e)
                        print('FAILED FOR: ' + pdb + chain)
                        for f in glob.glob("curate_PDB/raw/*"):
                            os.remove(f)
                        continue

        print(dict_uniprot)
        print('Assembling Multimers...')
        assemble_multimer(results_df)
        

                

        columns = ['resid', 'chain', 'pKa', 'sasa', 'PDB Code']
        all_pka_sasa_res = pd.DataFrame(columns=columns)
        num = 0
        files = list((glob.glob("assembled/*pdb")))

        if (len(files) != 0):
            results_df.index = pd.RangeIndex(len(results_df.index))

            results_df.index = range(len(results_df.index))

            

            for file in files:

                code = file[10:]


                if code[:2] == 'AF':
                    try:
                        AF_lysines_df = find_AF_plddt(code)

                    except Exception as e:
                        print("ERROR: %s"%e)
                        continue

                try:
                    pka_res_df = calculate_pKa(code)
                    pka_sasa_res_df = calculate_sasa(code, pka_res_df)
                    all_pka_sasa_res = all_pka_sasa_res.append(pka_sasa_res_df)
                    print(all_pka_sasa_res)


                except:
                    pass


                if code[:2] == 'AF':
                    df_merge = pd.merge(AF_lysines_df, all_pka_sasa_res, how='outer', on='resid')
                    df_merge.drop('chain_y', inplace=True, axis=1)
                    df_merge.rename(columns={'chain_x': 'chain'}, inplace=True)
                    df = df_merge
                    df['plddt'] = pd.to_numeric(df['plddt'],errors='coerce')
                    df = df.where(df['plddt'] > 70)
                    df = df.dropna()
                    all_pka_sasa_res = df
                else:
                    all_pka_sasa_res['plddt'] = 'N/A'

                all_pka_sasa_res['Uniprot Entry'] = all_pka_sasa_res['PDB Code'].map(dict_uniprot)
                print(all_pka_sasa_res)
                all_pka_sasa_res.to_csv('Output/results.csv')


                num = num + 1
                print(num)


                try:
                    code = code[:-4]
                    code_to_remove = code + '.pka'
                    os.remove(code_to_remove)
                    print('Files removed')

                except Exception as e:
                    print("ERROR: %s"%e)

            for f in glob.glob("assembled/*"):
                os.remove(f)
            for f in glob.glob('curate_PDB/clean/*'):
                os.remove(f)





        #Averages the pKa/sasa values for a specific residue for the PDB structures for a particular protein


        question_avgs = [
        inquirer.List('Choice',
            message="Do you want to average the data obtained from each of the PDB structures of a protein?",
                choices=['Yes', 'No'],
                    ),
        ]
        answer_avgs = inquirer.prompt(question_avgs)

        if answer_avgs['Choice'] == 'Yes':
            try:
                all_pka_sasa_res = average_prot(all_pka_sasa_res)
                print('AVERAGED DATA...')
                print(all_pka_sasa_res)
            except Exception as e:
                    print("ERROR: %s"%e)
                    print('Failed to average data.')
        
        else:
            pass






        
        #Plots pKa vs sasa on scatter plot.
        #Also allows any known carbamates to be marked.



        carbam_pdb_list = list()
        carbam_resid_list = list()
        list_of_ids = list()
        selecting_carbam_lys = True

        question_data_analysis = [
        inquirer.List('Choice',
            message="Would you like the data to be anaylsed?",
                choices=['Yes', 'No'],
                    ),
        ]

        answer_data_analysis = inquirer.prompt(question_data_analysis)


        if answer_data_analysis['Choice'] == 'Yes':
            question_carbam = [
            inquirer.List('Choice',
                        message="Are there carbamylated lysines you wish to mark?",
                        choices=['Yes', 'No'],
                        ),
            ]
            answer_carbam = inquirer.prompt(question_carbam)

            if answer_carbam['Choice'] == 'Yes':

                while selecting_carbam_lys == True:

                    list_of_pdbs = all_pka_sasa_res['PDB Code']
                    list_of_pdbs_no_dup = list()
                    for entry in list_of_pdbs:
                        if entry not in list_of_pdbs_no_dup:
                            list_of_pdbs_no_dup.append(entry)
                        else:
                            continue
                
                    list_of_resids = all_pka_sasa_res['resid']
                    list_of_resids_no_dup = list()
                    for entry in list_of_resids:
                        if entry not in list_of_resids_no_dup:
                            list_of_resids_no_dup.append(entry)
                        else:
                            continue

                    question_carbam_pdb = [
                    inquirer.List('Choice',
                                message="Which structure is the residue in?",
                                choices=list_of_pdbs_no_dup,
                                ),
                    ]
                    answer_carbam_pdb = inquirer.prompt(question_carbam_pdb)

                    carbam_pdb_list.append(answer_carbam_pdb['Choice'])
                    


                    question_carbam_resid = [
                    inquirer.List('Choice',
                                message="What is the Residue ID of the carbamylated residue?",
                                choices=list_of_resids_no_dup,
                                ),
                    ]
                    answer_carbam_resid = inquirer.prompt(question_carbam_resid)

                    carbam_resid_list.append(answer_carbam_resid['Choice'])


                    continue_q = [
                    inquirer.List('Choice',
                            message="Are there more carbamylated lysines to add?",
                            choices=['Yes', 'No'],
                            ),
                    ]
                    continue_ans = inquirer.prompt(continue_q)

                    if continue_ans['Choice'] == 'Yes':
                        continue

                    elif continue_ans['Choice'] == 'No':
                        selecting_carbam_lys = False

                
            try:
                all_pka_sasa_res = analyse_data(all_pka_sasa_res, carbam_pdb_list, carbam_resid_list)
            except Exception as e:
                    print("ERROR: %s"%e)
                
        print(all_pka_sasa_res)


    except KeyboardInterrupt:
        pass

    except Exception as e:
        print("Error: %s"%e)
        if skip == 1:
            print('\n' + 'Goodbye' + '\n')
            running = False