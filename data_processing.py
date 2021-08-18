from os import name
import numpy
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
import inquirer
import os

#This code calculates the mean and standard deviation for the sasa and pKa values from all the structures of a given protein.
def report_on_results(results_df, all_pka_sasa_res):
    try:
        start_pdb_codes_no_dup = []
        start_pdb_codes = results_df['PDB Code'].tolist()

        for code in start_pdb_codes:
            if code not in start_pdb_codes_no_dup:
                start_pdb_codes_no_dup.append(code)
            else:
                continue

        end_pdb_code_no_dup = []
        end_pdb_codes = all_pka_sasa_res['PDB Code'].tolist()
        for code in end_pdb_codes:
            if code not in end_pdb_code_no_dup:
                end_pdb_code_no_dup.append(code)
            else:
                continue

        percentage = (len(end_pdb_code_no_dup) / len(start_pdb_codes_no_dup))*100

    except Exception as e:
        print("Error %s"%e)

    return(percentage)




def average_prot(all_pka_sasa_res):
    print('Averaging Uniprot data...')

#First creates a dataframe and puts all uniprot entries in a list.

    try:
        columns = ['resid', 'chain', 'plddt', 'Uniprot Entry', 'pKa', 'pKa stdev', 'PDB Code', 'sasa', 'sasa stdev']
        avgd_pka_sasa = pd.DataFrame(columns=columns)
        list_of_uniprot_codes = all_pka_sasa_res['Uniprot Entry']
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
            df_one_uniprot_code = all_pka_sasa_res.where(all_pka_sasa_res['Uniprot Entry'] == uniprot_code)
            df_one_uniprot_code = df_one_uniprot_code[df_one_uniprot_code['Uniprot Entry'].notna()]
            list_chains = df_one_uniprot_code['chain']
            list_chains_no_dup = list()

#Next makes a list of all the chains in the structure

            for entry in list_chains:
                if entry not in list_chains_no_dup:
                    list_chains_no_dup.append(entry)
                else:
                    continue

        except Exception as e:
            print('Failed to construct dataframe for' + uniprot_code)
            print("ERROR: %s"%e)
            return()


#Next makes a df including all the lysines in a specific chain in all uniprot entries.

        for chain in list_chains_no_dup:

            try:
                df_one_chain = df_one_uniprot_code.where(df_one_uniprot_code['chain'] == chain)
                df_one_chain= df_one_chain[df_one_chain['chain'].notna()]

#Next makes a list of the lysine resids in that chain and gets rid of all duplicates.

                list_of_resids = df_one_chain['resid']
                list_of_resids_no_dup = list()
                for entry in list_of_resids:
                    if entry not in list_of_resids_no_dup:
                        list_of_resids_no_dup.append(entry)
                    else:
                        continue
                    
            except Exception as e:
                print('Failed to construct dataframe for chain ' + chain)
                print("ERROR: %s"%e)
                continue
#Next calculates the mean and standard deviation for each individual lysine's pKa/sasa in all the available structures and puts them in a df.

            for residue in list_of_resids_no_dup:
                try:
                    df_one_resid = df_one_chain.where(df_one_chain['resid'] == residue)
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

                    PDB_codes_avgd = ''
                    for i in range(len(list_of_pdbs_no_dup)):
                        PDB_codes_avgd = PDB_codes_avgd + '/' + list_of_pdbs_no_dup[i]
    
                    d = {'resid': residue, 'chain': chain, 'plddt': avg_plddt, 'Uniprot Entry': uniprot_code, 'pKa': avg_pka, 'pKa stdev':stddev_pka, 'PDB Code': PDB_codes_avgd, 'sasa': avg_sasa, 'sasa stdev': stdev_sasa}
                    avgd_pka_sasa = avgd_pka_sasa.append(d, ignore_index=True)
    
            
                except Exception as e:
                    print('Failed to construct final dataframe for residue ' + residue)
                    print("ERROR: %s"%e)
                    continue
                
    #Lastly appends into df and saves as Output/results.csv

        try:
            avgd_pka_sasa['sasa stdev'] = avgd_pka_sasa['sasa stdev'].fillna(0)
            avgd_pka_sasa['pKa stdev'] = avgd_pka_sasa['pKa stdev'].fillna(0)

            if not os.path.exists("Output"):
                os.mkdir("Output")

            avgd_pka_sasa.to_csv('Output/results.csv')
            
            print('Done')

        except Exception as e:
            print("ERROR: %s"%e)
            print('Failed to average data')
            return()


    return(avgd_pka_sasa)










#This module plots the results as a pKa vs sasa graph.
#If lysines are input it marks them as red on the graph

def analyse_data(all_pka_sasa_res, carbam_pdb_list, carbam_resid_list):
    print(carbam_pdb_list)
    print(carbam_resid_list)

    #Firstly it converts the resids to intergers.
    try:
        list_of_ids = list()
        all_pka_sasa_res = all_pka_sasa_res.astype({"resid": int})

    except Exception as e:
        print("ERROR: %s"%e)

#Next it creates another column where the known carbamylated lysines are marked with 'True' and the others with 'False'

    for i in range(len(carbam_resid_list)):

        try:
            carbam_res = all_pka_sasa_res.where(all_pka_sasa_res['resid'] == int(carbam_resid_list[i]))
            carbam_res = carbam_res.where(carbam_res['PDB Code'] == carbam_pdb_list[i])

            carbam_res_df = carbam_res[carbam_res['resid'].notna()]
            print(carbam_res_df)
            list_ids = carbam_res_df.index.tolist()

        except Exception as e:
            print("ERROR: %s"%e)
            print('Failed to mark carbamate on resid = ' + str(i))

        try:
            list_of_ids = list_of_ids + list_ids
        except Exception as e:
                print("ERROR: %s"%e)

    try:
        all_pka_sasa_res['carbamylated'] = all_pka_sasa_res.index.isin(list_of_ids)
        print(all_pka_sasa_res)
    except Exception as e:
            print("ERROR: %s"%e)

#Next it plots the data.

    try:
        fig, ax = plt.subplots()
        colors = {True:'red', False:'black'}
        all_pka_sasa_res['sasa'] = all_pka_sasa_res['sasa'].astype(float)
        plt.scatter(all_pka_sasa_res['sasa'], all_pka_sasa_res['pKa'], c=all_pka_sasa_res['carbamylated'].map(colors))
        plt.title('pKa vs sasa')
        plt.xlabel('sasa')
        plt.ylabel('pKa')
        plt.show()

    except Exception as e:
            print("ERROR: %s"%e)
            print('Error plotting data.')

    return(all_pka_sasa_res)


def get_most_likely_value(all_pka_sasa_res):

#First creates a dataframe and puts all uniprot entries in a list.

    try:

        columns = ['resid', 'chain', 'plddt', 'Uniprot Entry', 'pKa', 'PDB Code', 'sasa']
        low_pka_sasa = pd.DataFrame(columns=columns)
        list_of_uniprot_codes = all_pka_sasa_res['Uniprot Entry']
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
            df_one_uniprot_code = all_pka_sasa_res.where(all_pka_sasa_res['Uniprot Entry'] == uniprot_code)
            df_one_uniprot_code = df_one_uniprot_code[df_one_uniprot_code['Uniprot Entry'].notna()]
            list_chains = df_one_uniprot_code['chain']
            list_chains_no_dup = list()

#Next makes a list of all the chains in the structure

            for entry in list_chains:
                if entry not in list_chains_no_dup:
                    list_chains_no_dup.append(entry)
                else:
                    continue

        except Exception as e:
            print('Failed to construct dataframe for' + uniprot_code)
            print("ERROR: %s"%e)
            return()


#Next makes a df including all the lysines in a specific chain in all uniprot entries.

        for chain in list_chains_no_dup:

            try:
                df_one_chain = df_one_uniprot_code.where(df_one_uniprot_code['chain'] == chain)
                df_one_chain= df_one_chain[df_one_chain['chain'].notna()]

#Next makes a list of the lysine resids in that chain and gets rid of all duplicates.

                list_of_resids = df_one_chain['resid']
                list_of_resids_no_dup = list()
                for entry in list_of_resids:
                    if entry not in list_of_resids_no_dup:
                        list_of_resids_no_dup.append(entry)
                    else:
                        continue
                    
            except Exception as e:
                print('Failed to construct dataframe for chain ' + chain)
                print("ERROR: %s"%e)
                continue
#Next calculates the mean and standard deviation for each individual lysine's pKa/sasa in all the available structures and puts them in a df.

            for residue in list_of_resids_no_dup:
                
                try:
                    df_one_resid = df_one_chain.where(df_one_chain['resid'] == residue)
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



            if not os.path.exists("Output"):
                os.mkdir("Output")

            #low_pka_sasa.to_csv('Output/results.csv')
            
            print('Done')
            print(low_pka_sasa)
    return(low_pka_sasa)







if __name__ == "__main__":




    try:
        all_pka_sasa_res = pd.read_csv('Output/results.csv')
        all_pka_sasa_res = all_pka_sasa_res.drop(['Unnamed: 0'], axis=1)
        carbam_pdb_list = []
        carbam_resid_list = []
        get_most_likely_value(all_pka_sasa_res)
        
    except Exception as e:
        print("ERROR: %s"%e)