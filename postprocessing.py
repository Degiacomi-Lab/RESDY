import os
import re
import subprocess

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker

import biobox as bb

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
        percentage = (len(end_pdb_code_no_dup) / float(len(start_pdb_codes_no_dup)))*100.0

    except Exception as e:
        print("Error: %s"%e)
        percentage = 0.0

    return percentage


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
        return

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
        return
    
#Next constructs a df for all the lysines of all the structures of a given uniprot entry.
    for uniprot_code in list_of_uniprot_codes_no_dup:
        print('Uniprot')

        try:
            df_one_uniprot_code = pka_sasa_results.where(pka_sasa_results['Uniprot Entry'] == uniprot_code)
            df_one_uniprot_code = df_one_uniprot_code[df_one_uniprot_code['Uniprot Entry'].notna()]
            
        except Exception as e:
            print('Failed to construct dataframe for' + uniprot_code)
            print("ERROR: %s"%e)
            return

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
                #print(list_of_homomers)


#For each group of equivalent chains it then creates a df.
            for entry in list_of_homomers:
                #print('homomers')
                try:
                    #print(entry)
                    boolean_series = df_one_pdb.chain.isin(entry)
                    df_only_homomers = df_one_pdb[boolean_series]
                    
                    if len(df_only_homomers) == 0:
                        #print(df_only_homomers)
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
                        #print(residue)
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
                
            avgd_pka_sasa.to_csv(os.path.join("Output", "results_avg.csv"))

        except Exception as e:
            print("ERROR: %s"%e)
            print('Failed to average data')
            return

    return avgd_pka_sasa


#This function plots the results as a pKa vs SASA graph.
#If lysines are input it marks them as red on the graph
def analyse_data(pka_sasa_results, carbam_pdb_list, carbam_resid_list):
    #print(carbam_pdb_list)
    #print(carbam_resid_list)

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
            #print(carbam_res_df)
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
        pka_sasa_results.to_csv(os.path.join('Output','results_carbamates_marked.csv'))
        #print(pka_sasa_results)
        
    except Exception as e:
            print("ERROR: %s"%e)

    #Next it plots the data.
    try:
        fig, ax = plt.subplots()
        colors = {True:'#68246D', False:'black'}
        alphas = {True:1, False:0.05}
        pka_sasa_results['sasa'] = pka_sasa_results['sasa'].astype(float)
        pka_sasa_results['pKa'] = pka_sasa_results['pKa'].astype(float)
        plt.scatter(pka_sasa_results['sasa'], pka_sasa_results['pKa'], c=pka_sasa_results['carbamylated'].map(colors))#, alpha=pka_sasa_results['carbamylated'].map(alphas))
        plt.title('pKa vs sasa')
        plt.xlabel('sasa')
        plt.ylabel('pKa')
        plt.show()
        pka_sasa_results.to_csv('RESULTS_marked.csv')
        #print('done')

    except Exception as e:
            print("ERROR: %s"%e)
            print('Error plotting data.')

    return pka_sasa_results


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
        return
    
#Next constructs a df for all the lysines of all the structures of a given uniprot entry.
    for uniprot_code in list_of_uniprot_codes_no_dup:

        try:
            df_one_uniprot_code = pka_sasa_results.where(pka_sasa_results['Uniprot Entry'] == uniprot_code)
            df_one_uniprot_code = df_one_uniprot_code[df_one_uniprot_code['Uniprot Entry'].notna()]
            
        except Exception as e:
            print('Failed to construct dataframe for' + uniprot_code)
            print("ERROR: %s"%e)
            return

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
                
            low_pka_sasa.to_csv(os.path.join('Output', 'results_likely.csv'))

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
        return
    
#Next constructs a df for all the lysines of all the structures of a given uniprot entry
    for uniprot_code in list_of_uniprot_codes_no_dup:

        try:
            df_one_uniprot_code = pka_sasa_results.where(pka_sasa_results['Uniprot Entry'] == uniprot_code)
            df_one_uniprot_code = df_one_uniprot_code[df_one_uniprot_code['Uniprot Entry'].notna()]
            

        except Exception as e:
            print('Failed to construct dataframe for' + uniprot_code)
            print("ERROR: %s"%e)
            return

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
                        low_pka_sasa.to_csv(os.path.join('Output', 'results_most_likely.csv'))

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

        if sys.platform == "win32":
            line = "curl -O " + file_name_fasta + " " + web_url
        else:
            line = "wget -O " + file_name_fasta + " " + web_url
            
        subprocess.check_call(line, shell=True)
        
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


def remove_problematic(code, propka_lys_fails, chain_resid_near_failed_chain, pka_sasa_res_df):

    def report_lys_fail(code, propka_lys_fails, chain_resid_clash, chain_resid_near_failed_chain):

        if not os.path.exists('Failed_Lysines'):
            os.mkdir('Failed_Lysines')

        code_no_pdb = code[:-4]
        
        path = os.path.join("Failed_Lysines", "%s_fails.txt"%code_no_pdb)
        #path = "Failed_Lysines%s%s_fails.txt"%(os.sep, code_no_pdb)        
        f = open(path, 'w')
        for entry in propka_lys_fails:
            f.write(entry)
            f.write('-Failed propka- check propka log for explanation')
            f.write('\n')

        for entry in chain_resid_clash:
            f.write(entry)
            f.write('-Clash error- likely due to autopatch adding a section which causes clashing')
            f.write('\n')

        for entry in chain_resid_near_failed_chain:
            f.write(entry)
            f.write('-Near chain error- likely due to autopatch adding a section which causes chains to overlap')
            f.write('\n')
        return


    def remove_failures(code, propka_lys_fails, chain_resid_near_failed_chain, pka_sasa_res_df):
    #The program then removes residues that have been selected as problematic along the line.
        print('Removing problematic residues...')
        print('\n')


        #Firstly, those which were selected as problematic in the propka report are removed.
        try:
            inverse_boolean_series = ~pka_sasa_res_df.Chain_Resid.isin(propka_lys_fails)
            pka_sasa_res_df = pka_sasa_res_df[inverse_boolean_series]
        except Exception as e:
            print('Error %s'%e)
            print('Error removing problematic propka residues')
            pka_sasa_res_df = pka_sasa_res_df[0:0]
            return(pka_sasa_res_df)

        #Next chain clash is investigated, if two chains clash the relevant residues are removed
        if code[:2] != 'AF':

            try:
                chain_resid_clash = check_clash(code)
                print(chain_resid_clash)
                inverse_boolean_series = ~pka_sasa_res_df.Chain_Resid.isin(chain_resid_clash)
                pka_sasa_res_df = pka_sasa_res_df[inverse_boolean_series]
                print(pka_sasa_res_df)

            except Exception as e:
                print('Error %s'%e)
                print('Error removing clashing atoms')
                pka_sasa_res_df = pka_sasa_res_df[0:0]
                return(pka_sasa_res_df)

            #Next residues are removed if they are exposed to any chain that failed the autopatch (as without the chain they are usually exposed to the results are not reliable).
            try:

                inverse_boolean_series = ~pka_sasa_res_df.Chain_Resid.isin(chain_resid_near_failed_chain)
                pka_sasa_res_df = pka_sasa_res_df[inverse_boolean_series]

            except Exception as e:
                print('Error %s'%e)
                print('Error removing files near failed chain')
                pka_sasa_res_df = pka_sasa_res_df[0:0]
                return(pka_sasa_res_df)
        else:
            chain_resid_clash = list()
        
        print('Residues removed:')
        print(propka_lys_fails)
        print(chain_resid_clash)
        print(chain_resid_near_failed_chain)

        return chain_resid_clash, pka_sasa_res_df


    #This module checks for chain clash.
    #Chain clash occurs when there is an issue during the assembly which means that two chains can partially overlap.
    #Lysines in clashing regions are ignored.
    def check_clash(pdb_code):
        list_of_chains = list()
        list_of_resid = list()

        #Firstly a df is constructed to store results in.
        columns_2 = ['PDB Code', 'Lysine Index', 'Score']
        score_df = pd.DataFrame(columns=columns_2)
        print('BREAKING UP MOLECULE')

        #Next the assembled file is opened in biobox and the data is stored in a df.
        try:
            M = bb.Molecule()
            path = os.path.join("assembled", pdb_code)

            M.import_pdb(path, include_hetatm=True)
            df = M.data

            #Next the corrdinates and index for every lysine is found.
            lys_coords, lys_idx = M.atomselect('*','LYS', 'CA', use_resname=True, get_index=True)

            #Next the chain and resid for each lysine are put into lists.
            for entry in lys_idx:
                chain = df.at[entry, 'chain']
                list_of_chains.append(chain)

            for entry in lys_idx:
                resid = df.at[entry, 'resid']
                list_of_resid.append(resid)

            #Next the coordinates and index for every atom are found.
            all_coords, idx = M.atomselect('*','*','*', get_index=True)

        except Exception as e:
            print('Error %s'%e)
            print('Failure checking chain clash')
            return []

        #The program then cycles through each lysine and finds the distnace to every atom.
        for j in range(len(lys_coords)):
            try:
                lys_chain = list_of_chains[j]
                score = 0

                try:          
                    for i in range(len(all_coords)):
                        x_dist = (((lys_coords[j])[0] - (all_coords[i])[0])**2)
                        y_dist = (((lys_coords[j])[1] - (all_coords[i])[1])**2)
                        z_dist = (((lys_coords[j])[2] - (all_coords[i])[2])**2)
                        distance = np.sqrt(x_dist + y_dist + z_dist)
                        
                        #If the distance to an atom on another chain is found the be less than 4.5 an extra point is added to the score.

                        if distance < 4.5:
                            index_of_aa = idx[i]
                            test_chain = df.at[index_of_aa, 'chain']

                            if test_chain != lys_chain:             
                                score = score + 1

                except Exception as e:
                    print('Error %s'%e)
                
                #The score for each lysine is then added to a df.
                chain_resid = str(list_of_chains[j]) + str(list_of_resid[j])

                data = ({'PDB Code': pdb_code, 'Lysine Index': lys_idx[j], 'Resid': list_of_resid[j], 'Chain': list_of_chains[j], 'Chain_Resid': chain_resid,'Score': score})
                score_df = score_df.append(data, ignore_index=True)

            except Exception as e:
                print('Error %s'%e)
                print('Failure chekcing chain clash')
                score = 0
                data = ({'PDB Code': pdb_code, 'Lysine Index': lys_idx[j], 'Resid': list_of_resid[j], 'Chain': list_of_chains[j], 'Chain_Resid': chain_resid,'Score': score})
                score_df = score_df.append(data, ignore_index=True)
            
                continue

            #If the score is greater than 1 the lysine is removed (actual removal happens in main) from the results.
            try:
                score_df = score_df.where(score_df['Score'] > 1)
                score_df = score_df.dropna()

                chain_resid_clash = score_df['Chain_Resid'].tolist()
            except Exception as e:
                print('Error %s'%e)
                print('Failure chekcing chain clash')
                return []

            return(chain_resid_clash)
    
    chain_resid_clash, pka_sasa_res_df = remove_failures(code, propka_lys_fails, chain_resid_near_failed_chain, pka_sasa_res_df)
    report_lys_fail(code, propka_lys_fails, chain_resid_clash, chain_resid_near_failed_chain)

    return(pka_sasa_res_df)

###############################################################

if __name__ == "__main__":

    import inquirer

    try:
        question = [
        inquirer.List('Choice',
                        message="How would you like the data to be processed?",
                        choices=['Average for each resid', 'Most likely to form carbamate for each resid', 'No processing needed', 'Plot input csv file', 'Only plot carbamates', 'Plot with errors', 'Only PDBs', 'Highest sasa lowest pka', 'Mark online structures', 'Show carbamates and mark af'],
                    ),
        ]
        answer = inquirer.prompt(question)

        if answer['Choice'] == 'Average for each resid':
            pka_sasa_results = pd.read_csv(os.path.join('Output','results','results.csv'))
            pka_sasa_results = pka_sasa_results.drop(['Unnamed: 0'], axis=1)
            carbam_pdb_list = []
            carbam_resid_list = []
            average_prot(pka_sasa_results)

        if answer['Choice'] == 'Most likely to form carbamate for each resid':
            pka_sasa_results = pd.read_csv(os.path.join('Output','results.csv'))
            pka_sasa_results = pka_sasa_results.drop(['Unnamed: 0'], axis=1)
            get_most_likely_value(pka_sasa_results)
        
        if answer['Choice'] == 'Plot input csv file':
            csv_file_name = input['Name of CSV file:']
            if csv_file_name[-4:] != '.csv':
                csv_file_name = csv_file_name + '.pdb'
            pka_sasa_results = pd.read_csv(os.path.join('Output','results_all_together','results_all.csv'))
            pka_sasa_results = pka_sasa_results.drop(['Unnamed: 0'], axis=1)
            fig, ax = plt.subplots()
            colors = {True:'#68246D', False:'black'}
            alphas = {True:1, False:0.05}
            pka_sasa_results['sasa'] = pka_sasa_results['sasa'].astype(float)
            pka_sasa_results['pKa'] = pka_sasa_results['pKa'].astype(float)
            plt.scatter(pka_sasa_results['sasa'], pka_sasa_results['pKa'], c=pka_sasa_results['carbamylated'].map(colors))#, alpha=pka_sasa_results['carbamylated'].map(alphas))
            plt.title('pKa vs sasa')
            plt.xlabel('sasa')
            plt.ylabel('pKa')
            plt.show()
            
    except Exception as e:
        print("ERROR: %s"%e)