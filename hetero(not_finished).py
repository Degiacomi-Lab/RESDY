import biobox as bb
import numpy as np
import pandas as pd
import urllib.request, urllib.parse, urllib.error
import re
from bs4 import BeautifulSoup

all_pka_sasa_res = pd.read_csv('Output/results_up_to_2508/results.csv')
all_pka_sasa_res = all_pka_sasa_res.drop(['Unnamed: 0'], axis=1)

#Need to make this better...
#This function produces a df containing the most likely to form a carbamate for each residue.
#If the protein is a homomer it takes the lowest from all chains, whereasif it is a heteromer it doesn't.

try:

    columns = ['resid', 'chain', 'plddt', 'Uniprot Entry', 'pKa', 'PDB Code', 'sasa']
    low_pka_sasa = pd.DataFrame(columns=columns)
    list_of_uniprot_codes = all_pka_sasa_res['Uniprot Entry']
    list_of_uniprot_codes_no_dup = list()


except Exception as e:
    print('Failed to construct dataframe')
    print("ERROR: %s"%e)


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

#Next constructs a df for all the lysines of all the structures of a given uniprot entry.

for uniprot_code in list_of_uniprot_codes_no_dup:

    try:
        df_one_uniprot_code = all_pka_sasa_res.where(all_pka_sasa_res['Uniprot Entry'] == uniprot_code)
        df_one_uniprot_code = df_one_uniprot_code[df_one_uniprot_code['Uniprot Entry'].notna()]
        list_chains = df_one_uniprot_code['chain'].to_list()
        list_ch_re = df_one_uniprot_code['Chain_Resid'].to_list()
        list_ch_re_no_dup = []

        for entry in list_ch_re:
            if entry not in list_ch_re_no_dup:
                list_ch_re_no_dup.append(entry)
                
            else:
                continue

        list_chains_no_dup = []

        for entry in list_chains:
            if entry not in list_chains_no_dup:
                list_chains_no_dup.append(entry)
            else:
                continue

#Here it checks to see if the protein is a homomer or heteromer

        x = True
        for j in range(len(list_chains_no_dup)):

                code = list_ch_re_no_dup[0]
                code_no_let = code[1:]
                test = list_chains_no_dup[j] + code_no_let

                code_2 = list_ch_re_no_dup[1]
                code_no_let_2 = code[1:]
                test_2 = list_chains_no_dup[j] + code_no_let_2

                if (test in list_ch_re) and (test_2 in list_ch_re):
                    continue



                else:
                    print(uniprot_code)
                    print(test)
                    print(test_2)
                    print(list_ch_re)
                    #If it is a hetereomer it does not take the lowest pKa value.
                    x = False
                    low_pka_sasa = low_pka_sasa.append(df_one_uniprot_code)
                    print('Breaking')
                    break


#If it is a homomer it take the entry corresponding to the lowest pKa value for each residue.
        if x == True:
            try:

                list_of_resids = df_one_uniprot_code['resid']
                list_of_resids_no_dup = list()
                for entry in list_of_resids:
                    if entry not in list_of_resids_no_dup:
                        list_of_resids_no_dup.append(entry)
                    else:
                        continue
                    
            except Exception as e:
                print("ERROR: %s"%e)
                continue

            for residue in list_of_resids_no_dup:
                if x == False:
                    print(residue)
                    print(list_of_resids_no_dup)
                try:
                    df_one_resid = df_one_uniprot_code.where(df_one_uniprot_code['resid'] == residue)
                    df_one_resid = df_one_resid[df_one_resid['resid'].notna()]
                except Exception as e:
                    print('Failed to construct dataframe for residue ' + residue)
                    print("ERROR: %s"%e)
                    continue
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

                try:
                    list_of_pdbs = df_one_resid['PDB Code']
                    list_of_pdbs_no_dup = list()

                    for entry in list_of_pdbs:
                        if entry not in list_of_pdbs_no_dup:
                            list_of_pdbs_no_dup.append(entry)
                        else:
                            continue

                    PDB_codes_avgd = ''
                    for i in range(len(list_of_pdbs_no_dup)):
                        PDB_codes_avgd = PDB_codes_avgd + '/' + list_of_pdbs_no_dup[i]

                    d = {'resid': residue, 'chain': chain, 'plddt': avg_plddt, 'Uniprot Entry': uniprot_code, 'pKa': lowest_pka, 'PDB Code': PDB_codes_avgd, 'sasa': sasa}
                    low_pka_sasa = low_pka_sasa.append(d, ignore_index=True)
                        
                
                except Exception as e:
                    print('Failed to construct final dataframe for residue ' + residue)
                    print("ERROR: %s"%e)
                    continue
            
#Lastly appends into df and saves as Output/results.csv


        low_pka_sasa.to_csv('Output/results_likely_test.csv')
        


    except Exception as e:
        print('Failed to construct dataframe for' + uniprot_code)
        print("ERROR: %s"%e)
    
