import re
import os
import glob
import subprocess
import urllib.request, urllib.parse, urllib.error

import inquirer
from bs4 import BeautifulSoup
import pandas as pd
import numpy as np

import uniprot_loader as ul
import alphafold_loader as af
import patcher
import measure
import postprocessing


#This is the main file from which all the different functions are called.
#The general structure is as follows:

#The user gives their input (either uniprot codes, pdb codes or organism information).
#get_initial_data is then called.
#get_initial data identifies all relevant pdb entries, then puts them through the clean_split module, which cleans them and splits them up into the various alternate conformations.
#This returns a df called pdb_codes_df, which contains every relevant pdb file, their resolution, the method obtained and the chains present.
#The files are saved in curate_PDB/conformations

#Next they go through get_data and the autopatcher which patches where needed and returns clean files for each chain in curate_PDB/clean.
#Each pdbs full structure is then reassembled by the assemble_multimer module.

#The pKa of every lysine in the multimer is then calculated.
#The sasa of every lysine in the multimer is then also calculated, using a method where only the atoms close to each lysine are considered.
#These results are then placed in the df 'pka_sasa_results'.

#Next the data is processed for all the different structures of a protein (either the most likely to form a carbamate value for each resid is taken or an average).
#Lastly the data is plotted on a graph.
skip = 0
running = True
list_of_techniques = list()

while running:
    try:
        list_of_pdbs = list()
        print('\n' + '------------------------------------------------------------')
        print('\n' + '  Hello, press Ctrl + C at anytime to return to the start' + '\n')
        print('------------------------------------------------------------' + '\n')

        #Firstly it lets the user chose how they want to input data.

        questions = [
        inquirer.List('Choice',
                        message="Do you want to search a organism's whole proteome or input specific PDB or Uniprot Codes?",
                        choices=['Whole Proteome', 'Input PDB Codes', 'Input Uniprot Codes', 'Input Codes Via .csv file', 'Quit'],
                    ),
        ]
        answers = inquirer.prompt(questions)
        columns = ['Uniprot Entry', 'PDB Code', 'Method Structure Obtained by', 'Resolution', 'Chains']
        pdb_codes_df = pd.DataFrame(columns=columns)

        if answers['Choice'] != 'Quit':

#It then asks the user if they want to wipe the log file (i.e. do you want to continue where the last left off or not).
            question_log_file = [
            inquirer.List('Choice',
                            message="Do you want to wipe the previous log file?",
                            choices=['Yes', 'No'],
                        ),
            ]
            answer_log_file = inquirer.prompt(question_log_file)

            if answer_log_file['Choice'] == 'Yes':
                try:
                    if os.path.exists('log_file.csv'):
                        os.remove('log_file.csv')
                    else:
                        print('No log file found to remove')
                except Exception as e:
                    print('Failed to wipe log file')
            

#If there isn't already a log file it writes a new one and writes the column headers.

                try:

                    if not os.path.exists('log_file.csv'):
                        f = open('log_file.csv', 'w')
                        f.write('PDB Code,Result')
                        f.close()
                    log_file_df = pd.read_csv('log_file.csv')
                    done_pdbs = log_file_df['PDB Code'].to_list()

                except Exception as e:
                    print('Error parsing log_file.csv')
                    done_pdbs = []
                    pass



        #Allows user to chose if they only want to select structures obtained by certain techniques/of certain resolution.
        #e.g. the user may want to only look at structures obtained by X-ray diffraction and with a resolution less than 3 angstroms

        #It first asks the user for their preferences.

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



    #Searches Uniprot for PDB codes
    #Firstly asks for the name of the organismn (note it has to be exactly how it is written on the uniprot page).
    #It then asks for the uniprot code (the code in the url on the proteome page - e.g. for homo sapiens the code it is UP000005640 from the URL https://www.uniprot.org/proteomes/UP000005640)
        if answers["Choice"] == 'Whole Proteome':

            name_of_organism = input('Name of organism:')
            name_of_organism = name_of_organism.replace(' ', '+')
            code = input('Uniprot proteome code:')

    #It then goes into get_pdbs and gets all the information for that organism
            try:
                pdb_codes_df = ul.get_pdbs(name_of_organism, code, pdb_codes_df, done_pdbs)
                print(pdb_codes_df)
            except:
                print('Try again')
                skip = 1
            if (len(pdb_codes_df) == 0):
                skip = 1
                
#If the user inputs the uniprot code alone it gets all the corresponding pdb codes and the AF code.

        if answers['Choice'] == 'Input Uniprot Codes':
#Firstly it makes a list of the uniprot codes from the user.
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
            
            
#Next it goes into get_data and gets relevant information/ downloads and clean all the pdb codes belonging to the uniprot codes given.

            for uniprot_code in list_of_UNIPROT_codes:
                print('Getting data for: ' + uniprot_code)
                pdb_codes_df = ul.get_pdbs_uniprot(uniprot_code, pdb_codes_df, done_pdbs)



#If the user selects input PDB codes it goes here.

        if answers['Choice'] == 'Input PDB Codes':
            select_strucs = True

#While select_strucs if true the user can input pdb codes

            while select_strucs == True:

                question_4 = [
                inquirer.List('Choice',
                                message="Do you want to use PDB or AlphaFold structure(s)?",
                                choices=['PDB', 'AlphaFold'],
                            ),
                ]
                answer_4 = inquirer.prompt(question_4)

#It asks for the pdb code and uniprot code if PDB codes are being put in or just the alphafold code if alphafold structures are being put in.
#For each pdb code it downloads and cleans the structure and constructs a df consisting of the pdb code, uniprot code, method obtained and resolution.

                if answer_4["Choice"] == 'PDB':
                    try:
                        PDBCODE_inpt = input('What is the PDB code?')
                        PDBCODE_inpt = PDBCODE_inpt.upper()
                        print(PDBCODE_inpt)
                        UNIPROT_code_pdb = input('What is the Uniprot Code?')
                        pdb_codes_df = ul.construct_single_pdb_df(UNIPROT_code_pdb, PDBCODE_inpt, pdb_codes_df, done_pdbs)
                        print(pdb_codes_df)


                    except Exception as e:
                        print("ERROR: %s"%e)
                        print('Unable to obtain data for entry.')
                        pass

                elif answer_4['Choice'] == 'AlphaFold':
                    try:
                        AF_code = input('What is the Uniprot code?')

                        keep = True
                        for entry in done_pdbs:
                            if re.search(pdb, entry):
                                keep = False

                        if keep == True:
                            AF_code_full = 'AF-' + AF_code + '-F1-model_v1'
                            d = {'Uniprot Entry': AF_code, 'PDB Code': AF_code_full, 'Method Structure Obtained by': 'Predicted', 'Resolution': 'N/A', 'Chains': 'A'}
                            pdb_codes_df = pdb_codes_df.append(d, ignore_index=True)
                            print(pdb_codes_df)

                    except Exception as e:
                        print("ERROR: %s"%e)
                
#It then asks if there are more to add, if there are it asks for another code, if there aren't it moves on.

                question_more_pdbs = [
                inquirer.List('Choice',
                                message="Are there more to add?",
                                choices=['Yes', 'No'],
                            ),
                ]
                answer_4 = inquirer.prompt(question_more_pdbs)

                if answer_4['Choice'] == 'No':
                    select_strucs = False
                else:
                    pass


#If the user wants to input through a csv file it follows the same processes as for uniprot code and pdb code inputs.
#CSV files should have the following format:
#Uniprot Entry PDB Code
#P1234         1ABC   ----- to look specifically at that pdb code
#P1234                 ----- to look at all pdbs for that uniprot entry
#P1234         AF      ---- to only look at the alphafold structure

#The csv file is then parsed and the codes fed through as if they were entered manually.
        if answers['Choice'] == 'Input Codes Via .csv file':
            csv_name = input('What is the name of the .csv file?')
            if csv_name[-4:] != '.csv':
                csv_name = csv_name + '.csv'
            pdb_codes_df = ul.from_csv_file(csv_name, pdb_codes_df, done_pdbs)
            print(pdb_codes_df)



#If the user selects quit all pdb files are wiped.

        elif answers ['Choice'] == 'Quit':
            skip = 1            
            for f in glob.glob("curate_PDB%sclean%s*"%(os.sep, os.sep)):
                try:
                    os.remove(f)
                except:
                    continue
            for f in glob.glob("curate_PDB%sconformations%s*"%(os.sep, os.sep)):
                try:
                    os.remove(f)
                except:
                    continue
            for f in glob.glob("curate_PDB%sraw%s*"%(os.sep, os.sep)):
                try:
                    os.remove(f)
                except:
                    continue
            for f in glob.glob("assembled%s*"%os.sep):
                try:
                    os.remove(f)
                except:
                    continue
            raise Exception



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

#The desired techniques are appended to a list.

#It then filters the df to remove those which don't fir the criteria set by the user.

            raw_res = ul.search_by_technique(list_of_techniques, pdb_codes_df, wanted_res)
            pdb_codes_df = raw_res[0]
            print(pdb_codes_df)
            skip = raw_res[1]


        else:
            pass
    

        #Next the dictionary dict_uniprot is created, this is used to keep track of which pdb belongs to which uniprot entry and eventually map them to each other.

        dict_uniprot = dict()

        #It then creates the pka_sasa_results df (the final df that will hold the results).

        columns = ['resid', 'chain', 'pKa', 'sasa', 'PDB Code', 'Chain_Resid']
        pka_sasa_results = pd.DataFrame(columns=columns)
        num = 0

        #A list of pdb codes from the pdb_codes_df (the df with all the pdbs in) is then created and duplicates taken out.

        list_of_pdb_codes_no_dup = list()
        list_of_pdb_codes = pdb_codes_df['PDB Code']
        for entry in list_of_pdb_codes:
            if entry not in list_of_pdb_codes_no_dup:
                list_of_pdb_codes_no_dup.append(entry)
            else:
                continue
        
        #Each one is then fed through one at a time.

        for pdb in list_of_pdb_codes_no_dup:

            try:
                fail_dict = dict()
                print(list_of_pdb_codes_no_dup)
                df_one_pdb_code = pd.DataFrame()
        #A seperate df is then created which for each pdb code.
                df_one_pdb_code = pdb_codes_df[pdb_codes_df['PDB Code'] == pdb]
                df_one_pdb_code = df_one_pdb_code.reset_index(drop=True)
                df_one_pdb_code = df_one_pdb_code[df_one_pdb_code['PDB Code'].notna()]
                list_chains = df_one_pdb_code['Chains'].to_list()

                gap_dict = dict()
        #Each entry in the df is then iterated through.

                for i in range(len(df_one_pdb_code)):
                    print(df_one_pdb_code)
                    chain = df_one_pdb_code.at[i, 'Chains']
                    uniprot = df_one_pdb_code.at[i, 'Uniprot Entry']
        #If the code is an alphacode structure it is added into the dictionary dict_uniprot here.

                    if pdb[:2] != 'AF':
                        pdb_for_dict = pdb + '_assembled'
                        dict_uniprot.update({pdb_for_dict: uniprot})
                    else:
                        dict_uniprot.update({pdb: uniprot})

                    print('***************** GETTING DATA FOR ' + str(pdb) + str(chain) + '*******************')
                    print(chain)

        #If the code is an alphacode structure it is downloaded into the assembled file here.

                    if pdb[:2]=='AF':
                        try:
                            af.download_AF_struc(pdb)
                        except:
                            continue

        #If it is a pdb code it is fed through get_data (and the autopatcher), which takes the file from curate_PDB/conformations, patches it, splits it into chains and saves it to curate_PDB/clean.
        #Note the output from get_data (gap) is used to exclude points which could potentially be in contact with a loop that is patched by the autopatcher.

                    else:
                        try:
        
                            gap = patcher.get_data(pdb, chain)
                            gap_dict[chain] = gap
                            print('****************')
                            print(gap_dict)

        #The function correct_resid corrects the resid values in the clean files so that they are all in the correct position (as sometimes autopatcher produces some which are shifted).
                            
                            patcher.correct_resid(pdb, chain)

                        except Exception as e:
                            print('Error: %s'%e)
                            print('FAILED FOR: ' + pdb + chain)
                            continue
                        try:
                            for f in glob.glob("curate_PDB%sraw%s*"%(os.sep, os.sep)):
                                os.remove(f)
                        except:
                            pass

                    print(chain)

                print(dict_uniprot)

                if pdb[:2] == 'AF':
                    pass
                else:
        #Next the multimer is assembled by taking all the individual chain's pdbs in curate_PDB/clean and assembling them into one protein (using biobox).
                    print('Assembling Multimer...')
                    chain_resid_near_failed_chain = patcher.assemble_multimer(gap_dict, pdb, list_chains)
        
#The code next cycles through the files in the assembled folder.

                files = list((glob.glob("assembled%s*pdb"%(os.sep))))

                if (len(files) != 0):
                    pdb_codes_df.index = pd.RangeIndex(len(pdb_codes_df.index))

                    pdb_codes_df.index = range(len(pdb_codes_df.index))

                    

                    for file in files:

                        code = file[10:]

#If the structure is an alphafold structure it finds the plddt value (a measure of each residue's error in the structures prediction) and puts it in a dictionary.

                        if code[:2] == 'AF':
                            try:
                                dict_plddt = af.find_AF_plddt(code)

                            except Exception as e:
                                print("ERROR: %s"%e)
                                continue

#The functions calculate_pKa and break_up_and_calculate_sasa are then used to calculate the pKa and sasa respectively.
#Each function produces a dataframe- these are later merged together to form one results df for the structure (pKa_sasa_res_df)

                        try:

#The calculate_pKa function also produces propka_lys_fails- a list of chain_resids (chain + resid e.g. A12) to be removed as they appear in the propka output file.
                            print('Obtaining data...')
                            pka_results_df, propka_lys_fails = measure.calculate_pKa(code)
                            print('Obtaining sasa data...')
                            sasa_results_df = measure.break_up_and_calculate_sasa(code)
                            print('Got pKa and sasa')

                            try:
                        #The reuslting dfs from the sasa and pka calculators are merged and any columns not needed dropped.
                        ##########POTENTIALLY COULD IMPROVE THE MERGE SO I DON'T NEED TO GET RID OF OTHER COLUMNS####################

                                pka_results_df["resid"] = pd.to_numeric(pka_results_df["resid"], downcast="float")

                                pka_sasa_res_df = pd.merge(sasa_results_df, pka_results_df, how='inner', on=['Chain_Resid', 'chain', 'resid'])
                                print(pka_sasa_res_df)
                                #pka_sasa_res_df.drop('chain_x', inplace=True, axis=1)
                                #pka_sasa_res_df.drop('resid_x', inplace=True, axis=1)
                                #pka_sasa_res_df.rename(columns={'chain_y': 'chain'}, inplace=True)
                                #pka_sasa_res_df.rename(columns={'resid_y': 'resid'}, inplace=True)

#If the structure is an AF structure the plddt data is added to the df.
#For pdb structures the value is assumed to be 0.

                                if code[:2] == 'AF':

                                    pka_sasa_res_df['plddt'] = pka_sasa_res_df['Chain_Resid'].map(dict_plddt)
                                    pka_sasa_res_df['plddt'] = pka_sasa_res_df['plddt'].fillna(0)
                                    pka_sasa_res_df['plddt'] = pd.to_numeric(pka_sasa_res_df['plddt'],errors='coerce')
                                    pka_sasa_res_df = pka_sasa_res_df.where(pka_sasa_res_df['plddt'] > 70)
                                    pka_sasa_res_df = pka_sasa_res_df[pka_sasa_res_df['plddt'].notna()]

                                else:
                                    pka_sasa_res_df['plddt'] = '0'


                            except Exception as e:
                                print('Error %s'%e)
                                print('Failed to put data into dataframe for ' + file)
                                continue
                            

                           
                            
                            if code[:2] =='AF':
                                chain_resid_near_failed_chain = list()

                            try:

#This function removes any problematic entries from pka_sasa_res_df and returns the updated df.

                                pka_sasa_res_df = postprocessing.remove_problematic(code, propka_lys_fails, chain_resid_near_failed_chain, pka_sasa_res_df)


                            except Exception as e:
                                print('Error %s'%e)
                                print('Failure removing problematic values')
                                continue


                            pka_sasa_results = pka_sasa_results.append(pka_sasa_res_df)


                        except:
                            continue

#Next the uniprot code is mapped against the pdb code into the df from dict uniprot.

                        pka_sasa_results['Uniprot Entry'] = pka_sasa_results['PDB Code'].map(dict_uniprot)
                        pka_sasa_results.to_csv('Output%sresults.csv'%os.sep)
                        print(pka_sasa_results)

#Next it is reported in the log_file that the structure passed.

                        if not os.path.exists('log_file.csv'):
                            f = open('log_file.csv', 'w')
                            f.write('PDB Code,Result')
                            f.write('\n')
                            f.write(pdb + ',Passed')
                            f.close()
                        elif os.path.exists('log_file.csv'):
                            f = open('log_file.csv', 'a')
                            f.write('\n')
                            f.write(pdb + ',Passed')
                            f.close()
                        
                        print('done')


#If the process fails, it is noted in the log_file
            except Exception as e:
                print(pdb + ' Failed')
                print('Error %s'%e)
                if not os.path.exists('log_file.csv'):
                    f = open('log_file.csv', 'w')
                    f.write('PDB Code,Result')
                    f.write('\n')
                    f.write(pdb + ', Failed')
                    f.close()
                elif os.path.exists('log_file.csv'):
                    f = open('log_file.csv', 'a')
                    f.write('\n')
                    f.write(pdb + ', Failed collecting pKa/SASA data')
                    f.close()

                continue

#Next removes .pka file and those in assembled, clean and raw 

            finally:
                try:
                    code = code[:-4]
                    code_to_remove = code + '.pka'
                    os.remove(code_to_remove)
                    print('Files removed')

                except:
                    pass


            for f in glob.glob("assembled%s*"%os.sep):
                os.remove(f)
            for f in glob.glob("curate_PDB%sclean%s*"%(os.sep, os.sep)):
                file_name = f[17:]
                if file_name[:4] == pdb:
                    os.remove(f)
            for f in glob.glob("curate_PDB%sraw%s*"%(os.sep, os.sep)):
                os.remove(f)

#Once it has done gone through all the files it removes everything from the conformations file.

        for f in glob.glob("curate_PDB%sconformations%s*"%(os.sep, os.sep)):
            os.remove(f)

#This reports what percentage of files pass the process and gives a percentage pass rate.

        percentage_passed = postprocessing.report_on_results(pdb_codes_df, pka_sasa_results)
        print('Percentage passed = ' + str(percentage_passed) + '%')



#Next it starts the data processing.

#Firstly it asks the user how they would like the data to be processed.
#Average by resid gives the mean sasa and pka values for each resid.
#Take most likely for each resid gives the entry with the lowest pKa for each resid.
#Keep data raw changes nothing.

        question_avgs = [
        inquirer.List('Choice',
            message="How would you like the data to be processed?",
                choices=['Average by resid', 'Take most likely for each resid', 'Keep data raw'],
                    ),
        ]
        answer_avgs = inquirer.prompt(question_avgs)

        if answer_avgs['Choice'] == 'Average by resid':
            try:
                print('AVERAGED DATA...')
                pka_sasa_results = postprocessing.average_prot(pka_sasa_results)
                print(pka_sasa_results)
            except Exception as e:
                    print("ERROR: %s"%e)
                    print('Failed to average data.')
                    pass

        if answer_avgs['Choice'] == 'Take most likely for each resid':

            try:
                print('TAKING MOST LIKELY FOR EACH RESID')
                pka_sasa_results = postprocessing.get_most_likely_value(pka_sasa_results)
                print(pka_sasa_results)
            except Exception as e:
                print("ERROR: %s"%e)
                print('Failed to take most likely for each resid.')
                pass
        else:
            pass

        
        #The next section plots pKa vs sasa on scatter plot.
        #Also allows any known carbamates to be marked (they will appear as a different colour on the plot).
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

                    list_of_pdbs = pka_sasa_results['PDB Code']
                    list_of_pdbs_no_dup = list()
                    for entry in list_of_pdbs:
                        if entry not in list_of_pdbs_no_dup:
                            list_of_pdbs_no_dup.append(entry)
                        else:
                            continue


                    question_carbam_pdb = [
                    inquirer.List('Choice',
                                message="Which structure is the residue in?",
                                choices=list_of_pdbs_no_dup,
                                ),
                    ]
                    answer_carbam_pdb = inquirer.prompt(question_carbam_pdb)
                    carbam_pdb = answer_carbam_pdb['Choice']
                    carbam_pdb_list.append(carbam_pdb)
                    
                    df_one_pdb = pka_sasa_results.where(pka_sasa_results['PDB Code'] == carbam_pdb)
                    list_of_resids = df_one_pdb['resid']
                    list_of_resids_no_dup = list()
                    for entry in list_of_resids:
                        if entry not in list_of_resids_no_dup:
                            list_of_resids_no_dup.append(entry)
                        else:
                            continue


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

#Next, the results are plotted.
#A scatter plot is drawn of sasa (x) vs pKa (y) with each entry as a point.
#Ideally the carbamates should have a high sasa and low pKa.

            try:
                pka_sasa_results = postprocessing.analyse_data(pka_sasa_results, carbam_pdb_list, carbam_resid_list)
            except Exception as e:
                    print("ERROR: %s"%e)
                    pass
                
        print(pka_sasa_results)

#Lastly all leftover files are removed.
        try:
            os.remove("gap_data.txt")
        except:
            pass

        try:
            os.remove("patch_data.txt")
        except:
            pass

    except KeyboardInterrupt:
        pass

    except Exception as e:
        print("Error: %s"%e)
        if skip == 1:
            print('\n' + 'Goodbye' + '\n')
            for f in glob.glob("curate_PDB%sclean%s*"%(os.sep, os.sep)):
                try:
                    os.remove(f)
                except:
                    continue
            for f in glob.glob("curate_PDB%sconformations%s*"%(os.sep, os.sep)):
                try:
                    os.remove(f)
                except:
                    continue
            for f in glob.glob("curate_PDB%sraw%s*"%(os.sep, os.sep)):
                try:
                    os.remove(f)
                except:
                    continue
            for f in glob.glob("assembled%s*"%os.sep):
                try:
                    os.remove(f)
                except:
                    continue
            running = False