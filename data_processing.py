from os import name
import numpy
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
import inquirer
import os
import re
import subprocess

#This code is called at the end to give an idea of how many of the pdb files that passed the initial extraction passed the rest of the process.
def report_on_results(pdb_codes_df, pka_sasa_results):
    try:
#Firstly it puts each pdb code from pdb_codes_df in a list then gets rid of duplicates.
        start_pdb_codes_no_dup = []
        start_pdb_codes = pdb_codes_df['PDB Code'].tolist()

        for code in start_pdb_codes:
            if code not in start_pdb_codes_no_dup:
                start_pdb_codes_no_dup.append(code)
            else:
                continue
#Next it does the same for pka_sasa_results
        end_pdb_code_no_dup = []
        end_pdb_codes = pka_sasa_results['PDB Code'].tolist()
        for code in end_pdb_codes:
            if code not in end_pdb_code_no_dup:
                end_pdb_code_no_dup.append(code)
            else:
                continue
#Lastly it gives a percentage pass rate by working out the number of entries in pka_sasa_results compared to pdb_codes_df

        percentage = (len(end_pdb_code_no_dup) / len(start_pdb_codes_no_dup))*100

    except Exception as e:
        print("Error %s"%e)

    return(percentage)




def average_prot(pka_sasa_results):
    print('Averaging Uniprot data...')
#First creates a dataframe and puts all uniprot entries in a list.

    try:

        columns = ['resid', 'chain', 'plddt', 'Uniprot Entry', 'pKa', 'PDB Code', 'sasa']
        avgd_pka_sasa = pd.DataFrame(columns=columns)
        list_of_uniprot_codes = pka_sasa_results['Uniprot Entry']
        list_of_uniprot_codes_no_dup = list()

    except Exception as e:
        print('Failed to construct dataframe')
        print("ERROR: %s"%e)
        return()

#Next removes all duplicates from the list of uniprot entries.

    try:
        for entry in list_of_uniprot_codes:
            if entry not in list_of_uniprot_codes_no_dup:
                list_of_uniprot_codes_no_dup.append(entry)
                
            else:
                continue

    except Exception as e:
        print('Failed to obtain Uniprot Entries from dataframe')
        print("ERROR: %s"%e)
        return()
#Next constructs a df for all the lysines of all the structures of a given uniprot entry.

    for uniprot_code in list_of_uniprot_codes_no_dup:
        print('Uniprot')

        try:
            df_one_uniprot_code = pka_sasa_results.where(pka_sasa_results['Uniprot Entry'] == uniprot_code)
            df_one_uniprot_code = df_one_uniprot_code[df_one_uniprot_code['Uniprot Entry'].notna()]
            

        except Exception as e:
            print('Failed to construct dataframe for' + uniprot_code)
            print("ERROR: %s"%e)
            return()

#Next it creates a new column called 'PDB' which includes just the pdb code (i.e. not the alt conformation information).
#This is so the most likely value is selected from all the alt conformations.

        try:
            df_one_uniprot_code['PDB'] = df_one_uniprot_code['PDB Code']
            df_one_uniprot_code['PDB'] = (df_one_uniprot_code['PDB Code']).str[:9]
            list_of_pdbs = df_one_uniprot_code['PDB']
            list_of_pdbs_no_dup = list()

#Next it makes a list of PDB codes from the PDB column (and gets rid of duplicates).

            for entry in list_of_pdbs:
                if entry not in list_of_pdbs_no_dup:
                    list_of_pdbs_no_dup.append(entry)
                else:
                    continue
                
        except Exception as e:
            print("ERROR: %s"%e)
            continue

#For each element in the list a new df is created with just information beloning to that pdb.

        for entry in list_of_pdbs_no_dup:
            print('pdb')
            df_one_pdb = df_one_uniprot_code.where(df_one_uniprot_code['PDB'] == entry)
            df_one_pdb = df_one_pdb[df_one_pdb['PDB'].notna()]
            
            if entry[:2] == 'AF':
                list_of_chains = df_one_pdb['chain'].to_list()
                list_of_chains_no_dup = list()
                for entry in list_of_chains:
                    if entry not in list_of_chains_no_dup:
                        list_of_chains_no_dup.append(entry)
                    else:
                        continue
                list_of_homomers = [list_of_chains_no_dup]
            else:
                list_of_homomers = check_chain_match(entry)
                print(list_of_homomers)

 

#For each group of equivalent chains it then creates a df.

            for entry in list_of_homomers:
                print('homomers')
                try:
                    print(entry)
                    boolean_series = df_one_pdb.chain.isin(entry)
                    df_only_homomers = df_one_pdb[boolean_series]
                    
                    if len(df_only_homomers) == 0:
                        print(df_only_homomers)
                        continue
                    
                    try:

#Next gets a list of resIDs from the df_only_homomers and gets rid of duplicates.

                        list_of_resids = df_only_homomers['resid']
                        list_of_resids_no_dup = list()
                        for entry in list_of_resids:
                            if entry not in list_of_resids_no_dup:
                                list_of_resids_no_dup.append(entry)
                            else:
                                continue
                            
                    except Exception as e:
                        print("ERROR: %s"%e)
                        continue

#Next makes a new df for each residue in the df_only_homomers

                    for residue in list_of_resids_no_dup:
                        print(residue)
                        try:
                            df_one_resid = df_only_homomers.where(df_only_homomers['resid'] == residue)
                            df_one_resid = df_one_resid[df_one_resid['resid'].notna()]
                        except Exception as e:
                            print('Failed to construct dataframe for residue ' + residue)
                            print("ERROR: %s"%e)
                            continue

                        try:
                            df_one_resid["sasa"] = pd.to_numeric(df_one_resid["sasa"], downcast="float")


                            avg_pka = df_one_resid['pKa'].mean()
                            avg_sasa = df_one_resid['sasa'].mean()
                            stddev_pka = df_one_resid['pKa'].std()
                            stdev_sasa = df_one_resid['sasa'].std()

                        except Exception as e:
                            print('Failed to average data for ' + residue)
                            print("ERROR: %s"%e)
                            continue

                        try:
                            avg_plddt = df_one_resid['plddt'].mean()
                   
                        except:
                            avg_plddt = 'N/A'

                        try:
                            list_of_pdbs = df_one_resid['PDB Code']
                            list_of_pdbs_no_dup = list()


                            for entry in list_of_pdbs:
                                if entry not in list_of_pdbs_no_dup:
                                    list_of_pdbs_no_dup.append(entry)
                                else:
                                    continue

                            list_chains_used_no_dup = list()
                            list_chains_used = df_one_resid['chain'].to_list()
                            for chain in list_chains_used:
                                if chain not in list_chains_used_no_dup:
                                    list_chains_used_no_dup.append(chain)
                                else:
                                    pass
                            chain_avgd = ''
                            for i in range(len(list_chains_used_no_dup)):
                                chain_avgd = chain_avgd +'/' + list_chains_used_no_dup[i]

                            PDB_codes_avgd = ''
                            for i in range(len(list_of_pdbs_no_dup)):
                                PDB_codes_avgd = PDB_codes_avgd + '/' + list_of_pdbs_no_dup[i]

                            d = {'resid': residue, 'chain': chain_avgd, 'plddt': avg_plddt, 'Uniprot Entry': uniprot_code, 'pKa': avg_pka, 'pKa stdev':stddev_pka, 'PDB Code': PDB_codes_avgd, 'sasa': avg_sasa, 'sasa stdev': stdev_sasa}
                            avgd_pka_sasa = avgd_pka_sasa.append(d, ignore_index=True)

                    
                        except Exception as e:
                            print('Failed to construct final dataframe for residue ' + str(residue))
                            print("ERROR: %s"%e)
                            continue
                except Exception as e:
                    print('Error %s'%e)
                    print('Failed to average data for residue ' + str(residue))
                    continue
            
    #Lastly appends into df and saves as Output/results_avg.csv

        try:
            avgd_pka_sasa['sasa stdev'] = avgd_pka_sasa['sasa stdev'].fillna(0)
            avgd_pka_sasa['pKa stdev'] = avgd_pka_sasa['pKa stdev'].fillna(0)

            if not os.path.exists("Output"):
                os.mkdir("Output")
            avgd_pka_sasa.to_csv('Output/results_avg.csv')
            
            print('Done')

        except Exception as e:
            print("ERROR: %s"%e)
            print('Failed to average data')
            return()


    return(avgd_pka_sasa)










#This module plots the results as a pKa vs sasa graph.
#If lysines are input it marks them as red on the graph

def analyse_data(pka_sasa_results, carbam_pdb_list, carbam_resid_list):
    print(carbam_pdb_list)
    print(carbam_resid_list)

    #Firstly it converts the resids to intergers.
    try:
        list_of_ids = list()
        pka_sasa_results = pka_sasa_results.astype({"resid": int})

    except Exception as e:
        print("ERROR: %s"%e)

#Next it creates another column where the known carbamylated lysines are marked with 'True' and the others with 'False'

    for i in range(len(carbam_resid_list)):

        try:

            pka_sasa_results = pka_sasa_results.reset_index(drop=True)

            carbam_res = pka_sasa_results.where(pka_sasa_results['resid'] == int(carbam_resid_list[i]))
            carbam_res = carbam_res.where(carbam_res['PDB Code'] == carbam_pdb_list[i])

            carbam_res_df = carbam_res[carbam_res['resid'].notna()]
            print(carbam_res_df)
            #THIS ISN"T WORKING- HAVE TO RESET INDEX

            list_ids = carbam_res_df.index.tolist()

        except Exception as e:
            print("ERROR: %s"%e)
            print('Failed to mark carbamate on resid = ' + str(i))

        try:
            list_of_ids = list_of_ids + list_ids
        except Exception as e:
                print("ERROR: %s"%e)

    try:
        pka_sasa_results['carbamylated'] = pka_sasa_results.index.isin(list_of_ids)
        pka_sasa_results.to_csv('Output/results_carbamates_marked.csv')
        print(pka_sasa_results)
    except Exception as e:
            print("ERROR: %s"%e)

#Next it plots the data.

    try:
        fig, ax = plt.subplots()
        colors = {True:'#68246D', False:'black'}
        alphas = {True:1, False:0.05}
        pka_sasa_results['sasa'] = pka_sasa_results['sasa'].astype(float)
        pka_sasa_results['pKa'] = pka_sasa_results['pKa'].astype(float)
        plt.scatter(pka_sasa_results['sasa'], pka_sasa_results['pKa'], c=pka_sasa_results['carbamylated'].map(colors), alpha=pka_sasa_results['carbamylated'].map(alphas))
        plt.title('pKa vs sasa')
        plt.xlabel('sasa')
        plt.ylabel('pKa')
        plt.show()
        pka_sasa_results.to_csv('RESULTS_marked.csv')
        print('done')

    except Exception as e:
            print("ERROR: %s"%e)
            print('Error plotting data.')

    return(pka_sasa_results)





#This returns a df where for each residue the most result that it most likely to be carbamylated is given.

def get_most_likely_value(pka_sasa_results):

#First creates a dataframe and puts all uniprot entries in a list.

    try:

        columns = ['resid', 'chain', 'plddt', 'Uniprot Entry', 'pKa', 'PDB Code', 'sasa']
        low_pka_sasa = pd.DataFrame(columns=columns)
        list_of_uniprot_codes = pka_sasa_results['Uniprot Entry']
        list_of_uniprot_codes_no_dup = list()

    except Exception as e:
        print('Failed to construct dataframe')
        print("ERROR: %s"%e)
        return()

#Next removes all duplicates from the list of uniprot entries.

    try:
        for entry in list_of_uniprot_codes:
            if entry not in list_of_uniprot_codes_no_dup:
                list_of_uniprot_codes_no_dup.append(entry)
                
            else:
                continue

    except Exception as e:
        print('Failed to obtain Uniprot Entries from dataframe')
        print("ERROR: %s"%e)
        return()
#Next constructs a df for all the lysines of all the structures of a given uniprot entry.

    for uniprot_code in list_of_uniprot_codes_no_dup:

        try:
            df_one_uniprot_code = pka_sasa_results.where(pka_sasa_results['Uniprot Entry'] == uniprot_code)
            df_one_uniprot_code = df_one_uniprot_code[df_one_uniprot_code['Uniprot Entry'].notna()]
            

        except Exception as e:
            print('Failed to construct dataframe for' + uniprot_code)
            print("ERROR: %s"%e)
            return()

#Next it creates a new column called 'PDB' which includes just the pdb code (i.e. not the alt conformation information).
#This is so the most likely value is selected from all the alt conformations.

        try:
            df_one_uniprot_code['PDB'] = df_one_uniprot_code['PDB Code']
            df_one_uniprot_code['PDB'] = (df_one_uniprot_code['PDB Code']).str[:9]
            list_of_pdbs = df_one_uniprot_code['PDB']
            list_of_pdbs_no_dup = list()

#Next it makes a list of PDB codes from the PDB column (and gets rid of duplicates).

            for entry in list_of_pdbs:
                if entry not in list_of_pdbs_no_dup:
                    list_of_pdbs_no_dup.append(entry)
                else:
                    continue
                
        except Exception as e:
            print("ERROR: %s"%e)
            continue

#For each element in the list a new df is created with just information beloning to that pdb.

        for entry in list_of_pdbs_no_dup:
            df_one_pdb = df_one_uniprot_code.where(df_one_uniprot_code['PDB'] == entry)
            df_one_pdb = df_one_pdb[df_one_pdb['PDB'].notna()]
            
            if entry[:2] == 'AF':
                list_of_chains = df_one_pdb['chain'].to_list()
                list_of_chains_no_dup = list()
                for entry in list_of_chains:
                    if entry not in list_of_chains_no_dup:
                        list_of_chains_no_dup.append(entry)
                    else:
                        continue
                list_of_homomers = [list_of_chains_no_dup]
            else:
                list_of_homomers = check_chain_match(entry)

 

#For each group of equivalent chains it then creates a df.

            for entry in list_of_homomers:
                try:
                    print(entry)
                    boolean_series = df_one_pdb.chain.isin(entry)
                    df_only_homomers = df_one_pdb[boolean_series]
                    print(df_only_homomers)
                    
                    try:

#Next gets a list of resIDs from the df_only_homomers and gets rid of duplicates.

                        list_of_resids = df_only_homomers['resid']
                        list_of_resids_no_dup = list()
                        for entry in list_of_resids:
                            if entry not in list_of_resids_no_dup:
                                list_of_resids_no_dup.append(entry)
                            else:
                                continue
                            
                    except Exception as e:
                        print("ERROR: %s"%e)
                        continue

#Next makes a new df for each residue in the df_only_homomers

                    for residue in list_of_resids_no_dup:
                        
                        try:
                            df_one_resid = df_only_homomers.where(df_only_homomers['resid'] == residue)
                            df_one_resid = df_one_resid[df_one_resid['resid'].notna()]
                        except Exception as e:
                            print('Failed to construct dataframe for residue ' + residue)
                            print("ERROR: %s"%e)
                            continue

#Next it sorts the data by pKa and selects the pKa and sasa data corresponding to the lowest pKa entry for each resID.

                        try:
                            residue = str(residue)
                            df_one_resid["pKa"] = pd.to_numeric(df_one_resid["pKa"], downcast="float")
                            sorted_pka_df = df_one_resid.sort_values(by=['pKa'], ascending=True)
                            sorted_pka_df = sorted_pka_df.reset_index(drop=True)
                            lowest_pka = sorted_pka_df.at[0, 'pKa']
                            chain = sorted_pka_df.at[0, 'chain']
                            lowest_pka = str(lowest_pka)
                            sasa = sorted_pka_df.at[0, 'sasa']
                            sasa = str(sasa)

                        except Exception as e:
                            print('Failed to average data for ' + residue)
                            print("ERROR: %s"%e)
                            continue

                        try:
                            avg_plddt = df_one_resid['plddt'].mean()
                        except:
                            avg_plddt = 'N/A'

#Next it gets a list of all the pdb codes that the entry includes.

                        try:
                            list_of_pdbs = df_one_resid['PDB Code']
                            list_of_pdbs_no_dup = list()

                            for entry in list_of_pdbs:
                                if entry not in list_of_pdbs_no_dup:
                                    list_of_pdbs_no_dup.append(entry)
                                else:
                                    continue

#It then concatonates them together to create a long string consisting of for example: /1ABC-alt-1A/1ABC-alt-1A/

                            PDB_codes_avgd = ''
                            for i in range(len(list_of_pdbs_no_dup)):
                                PDB_codes_avgd = PDB_codes_avgd + '/' + list_of_pdbs_no_dup[i]

#Lastly it appends all the information to a df and saves as in the output folder as results_likely.csv

                            d = {'resid': residue, 'chain': chain, 'plddt': avg_plddt, 'Uniprot Entry': uniprot_code, 'pKa': lowest_pka, 'PDB Code': PDB_codes_avgd, 'sasa': sasa}
                            low_pka_sasa = low_pka_sasa.append(d, ignore_index=True)
                            
                    
                        except Exception as e:
                            print('Failed to construct final dataframe for residue ' + residue)
                            print("ERROR: %s"%e)
                            continue
                except Exception as e:
                    print('Error %s'%e)
                    print('Failure finding most likely')
                    return
                



            low_pka_sasa.to_csv('Output/results_likely.csv')


        
        print('Done')
        print(low_pka_sasa)
    return(low_pka_sasa)



def highest_sasa_lowest_pka(pka_sasa_results):
#First creates a dataframe and puts all uniprot entries in a list.

    try:

        columns = ['resid', 'chain', 'plddt', 'Uniprot Entry', 'pKa', 'PDB Code', 'sasa']
        low_pka_sasa = pd.DataFrame(columns=columns)
        list_of_uniprot_codes = pka_sasa_results['Uniprot Entry']
        list_of_uniprot_codes_no_dup = list()

    except Exception as e:
        print('Failed to construct dataframe')
        print("ERROR: %s"%e)
        return()

#Next removes all duplicates from the list of uniprot entries.

    try:
        for entry in list_of_uniprot_codes:
            if entry not in list_of_uniprot_codes_no_dup:
                list_of_uniprot_codes_no_dup.append(entry)
                
            else:
                continue

    except Exception as e:
        print('Failed to obtain Uniprot Entries from dataframe')
        print("ERROR: %s"%e)
        return()
#Next constructs a df for all the lysines of all the structures of a given uniprot entry.

    for uniprot_code in list_of_uniprot_codes_no_dup:

        try:
            df_one_uniprot_code = pka_sasa_results.where(pka_sasa_results['Uniprot Entry'] == uniprot_code)
            df_one_uniprot_code = df_one_uniprot_code[df_one_uniprot_code['Uniprot Entry'].notna()]
            

        except Exception as e:
            print('Failed to construct dataframe for' + uniprot_code)
            print("ERROR: %s"%e)
            return()

#Next it creates a new column called 'PDB' which includes just the pdb code (i.e. not the alt conformation information).
#This is so the most likely value is selected from all the alt conformations.

        try:
            df_one_uniprot_code['PDB'] = df_one_uniprot_code['PDB Code']
            df_one_uniprot_code['PDB'] = (df_one_uniprot_code['PDB Code']).str[:9]
            list_of_pdbs = df_one_uniprot_code['PDB']
            list_of_pdbs_no_dup = list()

#Next it makes a list of PDB codes from the PDB column (and gets rid of duplicates).

            for entry in list_of_pdbs:
                if entry not in list_of_pdbs_no_dup:
                    list_of_pdbs_no_dup.append(entry)
                else:
                    continue
                
        except Exception as e:
            print("ERROR: %s"%e)
            continue

#For each element in the list a new df is created with just information beloning to that pdb.

        for entry in list_of_pdbs_no_dup:
            df_one_pdb = df_one_uniprot_code.where(df_one_uniprot_code['PDB'] == entry)
            df_one_pdb = df_one_pdb[df_one_pdb['PDB'].notna()]
            
            if entry[:2] == 'AF':
                list_of_chains = df_one_pdb['chain'].to_list()
                list_of_chains_no_dup = list()
                for entry in list_of_chains:
                    if entry not in list_of_chains_no_dup:
                        list_of_chains_no_dup.append(entry)
                    else:
                        continue
                list_of_homomers = [list_of_chains_no_dup]
            else:
                list_of_homomers = check_chain_match(entry)

 

#For each group of equivalent chains it then creates a df.

            for entry in list_of_homomers:
                try:
                    print(entry)
                    boolean_series = df_one_pdb.chain.isin(entry)
                    df_only_homomers = df_one_pdb[boolean_series]
                    print(df_only_homomers)
                    
                    try:

#Next gets a list of resIDs from the df_only_homomers and gets rid of duplicates.

                        list_of_resids = df_only_homomers['resid']
                        list_of_resids_no_dup = list()
                        for entry in list_of_resids:
                            if entry not in list_of_resids_no_dup:
                                list_of_resids_no_dup.append(entry)
                            else:
                                continue
                            
                    except Exception as e:
                        print("ERROR: %s"%e)
                        continue

#Next makes a new df for each residue in the df_only_homomers

                    for residue in list_of_resids_no_dup:
                        
                        try:
                            df_one_resid = df_only_homomers.where(df_only_homomers['resid'] == residue)
                            df_one_resid = df_one_resid[df_one_resid['resid'].notna()]
                        except Exception as e:
                            print('Failed to construct dataframe for residue ' + residue)
                            print("ERROR: %s"%e)
                            continue

                        try:

#It then selects the lowest pKa found and highest sasa found for all the entries in df_one_resid.

                            residue = str(residue)
                            df_one_resid["pKa"] = pd.to_numeric(df_one_resid["pKa"], downcast="float")
                            sorted_pka_df = df_one_resid.sort_values(by=['pKa'], ascending=True)
                            sorted_pka_df = sorted_pka_df.reset_index(drop=True)
                            lowest_pka = sorted_pka_df.at[0, 'pKa']

                            df_one_resid["sasa"] = pd.to_numeric(df_one_resid["sasa"], downcast="float")
                            sorted_sasa_df = df_one_resid.sort_values(by=['sasa'], ascending=False)
                            sorted_sasa_df = sorted_sasa_df.reset_index(drop=True)
                            highest_sasa = sorted_sasa_df.at[0, 'sasa']

                            chain = sorted_pka_df.at[0, 'chain']
                            lowest_pka = str(lowest_pka)
                            highest_sasa = str(highest_sasa)

                        except Exception as e:
                            print('Failed to average data for ' + residue)
                            print("ERROR: %s"%e)
                            continue

                        try:
                            avg_plddt = df_one_resid['plddt'].mean()
                        except:
                            avg_plddt = 'N/A'

                        try:
                            list_of_pdbs = df_one_resid['PDB Code']
                            list_of_pdbs_no_dup = list()

                            for entry in list_of_pdbs:
                                if entry not in list_of_pdbs_no_dup:
                                    list_of_pdbs_no_dup.append(entry)
                                else:
                                    continue

#It then concatonates them together to create a long string consisting of for example: /1ABC-alt-1A/1ABC-alt-1A/.

                            PDB_codes_avgd = ''
                            for i in range(len(list_of_pdbs_no_dup)):
                                PDB_codes_avgd = PDB_codes_avgd + '/' + list_of_pdbs_no_dup[i]

                            d = {'resid': residue, 'chain': chain, 'plddt': avg_plddt, 'Uniprot Entry': uniprot_code, 'pKa': lowest_pka, 'PDB Code': PDB_codes_avgd, 'sasa': highest_sasa}
                            low_pka_sasa = low_pka_sasa.append(d, ignore_index=True)
                            
                    
                        except Exception as e:
                            print('Failed to construct final dataframe for residue ' + residue)
                            print("ERROR: %s"%e)
                            continue
                        
    #Lastly appends into df and saves as Output/results.csv


                        low_pka_sasa.to_csv('Output/results_most_likely.csv')

                    #low_pka_sasa.to_csv('Output/results.csv')
                    
                    print('Done')
                    print(low_pka_sasa)

                except Exception as e:
                    print('Error %s'%e)
    
    return(low_pka_sasa)




#This function is important for the data processesing section
#It checks which chains in a protein are equivalent and tells the program to treat them as eqivalent (i.e. to take the lowest values for each resid from only the equivalent chains).

def check_chain_match(pdb):
    pdb = pdb[:4]
    try:
#Firstly it downloads the fasta file from online.
        web_url = "https://www.rcsb.org/fasta/entry/" + pdb + '/download'
        file_name_fasta = pdb + '.fasta'
        subprocess.check_call("wget -O " + file_name_fasta + " " + web_url, shell=True)
    except Exception as e:
        print('Error %s'%e)
        print('Failed to download fasta sequence')

#Next it opens the file, parses it and puts all the equivalent chains into a list of lists.
#For example if A and B are equivalent and so are C and D the list will read [[A, B], [C, D]]
#This is then fed back into the get_most_likely function.

    try:
        f = open(file_name_fasta, 'r')
        entries = list()
        for line in f:
            m = re.findall('Chain[ a-z , A-Z \[\]]*|', line)
            for entry in m:
                if len(entry) > 1:
                    entries.append(entry)
    except Exception as e:
        print('Error %s'%e)
        print('Failed to get chain information.')

    list_of_homomers = list()
    for entry in entries:
        try:
            clean_entry = list()
            words = entry.split(' ')
            
            for word in words:
                
                if word[:5] == 'Chain':
                    words.remove(word)
            for word in words:
                if len(word) > 1:
                    if word[1] == ',':
                        word = word[0]
                        clean_entry.append(word)
                if len(word) > 2:
                    if word[1] == '[':
                        word = word[0]
                        clean_entry.append(word)
                elif len(word) == 1:
                    clean_entry.append(word)

                clean_entry_no_dup = list()
                for entry in clean_entry:
                    if entry not in clean_entry_no_dup:
                        clean_entry_no_dup.append(entry)
                    else:
                        pass
            list_of_homomers.append(clean_entry_no_dup)
        except Exception as e:
            print('Error %s'%e)
            continue
        try:
            os.remove(file_name_fasta)
        except:
            pass
    print(list_of_homomers)
                    

    return list_of_homomers










if __name__ == "__main__":




    try:
        question = [
        inquirer.List('Choice',
                        message="How would you like the data to be processed?",
                        choices=['Average for each resid', 'Most likely to form carbamate for each resid', 'No processing needed', 'Plot input csv file', 'Only plot carbamates', 'Plot with errors', 'Only PDBs', 'Highest sasa lowest pka', 'Mark online structures', 'Show carbamates and mark af'],
                    ),
        ]
        answer = inquirer.prompt(question)

        if answer['Choice'] == 'Average for each resid':
            pka_sasa_results = pd.read_csv('Output/results/results.csv')
            pka_sasa_results = pka_sasa_results.drop(['Unnamed: 0'], axis=1)
            carbam_pdb_list = []
            carbam_resid_list = []
            average_prot(pka_sasa_results)


        if answer['Choice'] == 'Most likely to form carbamate for each resid':
            pka_sasa_results = pd.read_csv('Output/results.csv')
            pka_sasa_results = pka_sasa_results.drop(['Unnamed: 0'], axis=1)
            get_most_likely_value(pka_sasa_results)
        
        
        if answer['Choice'] == 'Plot input csv file':
            csv_file_name = input['Name of CSV file:']
            if csv_file_name[-4:] != '.csv':
                csv_file_name = csv_file_name + '.pdb'
            pka_sasa_results = pd.read_csv('Output/results_all_together/results_all.csv')
            pka_sasa_results = pka_sasa_results.drop(['Unnamed: 0'], axis=1)
            fig, ax = plt.subplots()
            colors = {True:'#68246D', False:'black'}
            alphas = {True:1, False:0.05}
            pka_sasa_results['sasa'] = pka_sasa_results['sasa'].astype(float)
            pka_sasa_results['pKa'] = pka_sasa_results['pKa'].astype(float)
            plt.scatter(pka_sasa_results['sasa'], pka_sasa_results['pKa'], c=pka_sasa_results['carbamylated'].map(colors), alpha=pka_sasa_results['carbamylated'].map(alphas))
            plt.title('pKa vs sasa')
            plt.xlabel('sasa')
            plt.ylabel('pKa')
            plt.show()
        
        
    except Exception as e:
        print("ERROR: %s"%e)
