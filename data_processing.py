from os import name
import numpy
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
import inquirer
import os
import re

#This code is called at the end to give an idea of how many of the pdb files that passed the initial extraction passed the rest of the process.
def report_on_results(results_df, all_pka_sasa_res):
    try:
#Firstly it puts each pdb code from results_df in a list then gets rid of duplicates.
        start_pdb_codes_no_dup = []
        start_pdb_codes = results_df['PDB Code'].tolist()

        for code in start_pdb_codes:
            if code not in start_pdb_codes_no_dup:
                start_pdb_codes_no_dup.append(code)
            else:
                continue
#Next it does the same for all_pka_sasa_res
        end_pdb_code_no_dup = []
        end_pdb_codes = all_pka_sasa_res['PDB Code'].tolist()
        for code in end_pdb_codes:
            if code not in end_pdb_code_no_dup:
                end_pdb_code_no_dup.append(code)
            else:
                continue
#Lastly it gives a percentage pass rate by working out the number of entries in all_pka_sasa_res compared to results_df

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


        except Exception as e:
            print('Failed to construct dataframe for' + uniprot_code)
            print("ERROR: %s"%e)
            return()





#Next makes a list of the lysine resids in that uniprot entry and gets rid of all duplicates.
        try:
            list_of_resids = df_one_uniprot_code['resid']
            list_of_resids_no_dup = list()
            for entry in list_of_resids:
                if entry not in list_of_resids_no_dup:
                    list_of_resids_no_dup.append(entry)
                else:
                    continue
                
        except Exception as e:
            print('Failed to construct dataframe for chain ' + uniprot_code)
            print("ERROR: %s"%e)
            continue
#Next calculates the mean and standard deviation for each individual lysine's pKa/sasa in all the available structures and puts them in a df.

        for residue in list_of_resids_no_dup:
            try:
                df_one_resid = df_one_uniprot_code.where(df_one_uniprot_code['resid'] == residue)
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

                d = {'resid': residue, 'chain': 'all', 'plddt': avg_plddt, 'Uniprot Entry': uniprot_code, 'pKa': avg_pka, 'pKa stdev':stddev_pka, 'PDB Code': PDB_codes_avgd, 'sasa': avg_sasa, 'sasa stdev': stdev_sasa}
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
            avgd_pka_sasa.to_csv('Output/results_avg.csv')
            
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
        all_pka_sasa_res.to_csv('Output/results_carbamates_marked.csv')
        print(all_pka_sasa_res)
    except Exception as e:
            print("ERROR: %s"%e)

#Next it plots the data.

    try:
        fig, ax = plt.subplots()
        colors = {True:'#68246D', False:'black'}
        alphas = {True:1, False:0.05}
        all_pka_sasa_res['sasa'] = all_pka_sasa_res['sasa'].astype(float)
        all_pka_sasa_res['pKa'] = all_pka_sasa_res['pKa'].astype(float)
        plt.scatter(all_pka_sasa_res['sasa'], all_pka_sasa_res['pKa'], c=all_pka_sasa_res['carbamylated'].map(colors), alpha=all_pka_sasa_res['carbamylated'].map(alphas))
        plt.title('pKa vs sasa')
        plt.xlabel('sasa')
        plt.ylabel('pKa')
        plt.show()
        all_pka_sasa_res.to_csv('RESULTS_marked.csv')
        print('done')

    except Exception as e:
            print("ERROR: %s"%e)
            print('Error plotting data.')

    return(all_pka_sasa_res)





#This returns a df where for each residue the most result that it most likely to be carbamylated is given.

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
            

        except Exception as e:
            print('Failed to construct dataframe for' + uniprot_code)
            print("ERROR: %s"%e)
            return()




#Next makes a list of the lysine resids in that chain and gets rid of all duplicates.

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

#Next calculates the mean and standard deviation for each individual lysine's pKa/sasa in all the available structures and puts them in a df.

        for residue in list_of_resids_no_dup:
            
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


            low_pka_sasa.to_csv('Output/results_likely.csv')

        #low_pka_sasa.to_csv('Output/results.csv')
        
        print('Done')
        print(low_pka_sasa)
    return(low_pka_sasa)



def highest_sasa_lowest_pka(all_pka_sasa_res):
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
            

        except Exception as e:
            print('Failed to construct dataframe for' + uniprot_code)
            print("ERROR: %s"%e)
            return()




#Next makes a list of the lysine resids in that chain and gets rid of all duplicates.

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

#Next calculates the mean and standard deviation for each individual lysine's pKa/sasa in all the available structures and puts them in a df.

        for residue in list_of_resids_no_dup:
            
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
    return(low_pka_sasa)











if __name__ == "__main__":




    try:
        question = [
        inquirer.List('Choice',
                        message="How would you like the data to be processed?",
                        choices=['Average for each resid', 'Most likely to form carbamate for each resid', 'No processing needed', 'Plot input csv file', 'Only plot carbamates', 'Plot with errors', 'Only PDBs', 'Highest sasa lowest pka'],
                    ),
        ]
        answer = inquirer.prompt(question)

        if answer['Choice'] == 'Average for each resid':
            all_pka_sasa_res = pd.read_csv('Output/results_2608/results.csv')
            all_pka_sasa_res = all_pka_sasa_res.drop(['Unnamed: 0'], axis=1)
            carbam_pdb_list = []
            carbam_resid_list = []
            average_prot(all_pka_sasa_res)
            all_pka_sasa_res = pd.read_csv('Output/results_avg.csv')
            all_pka_sasa_res = all_pka_sasa_res.drop(['Unnamed: 0'], axis=1)
            analyse_data(all_pka_sasa_res, carbam_pdb_list, carbam_resid_list)

        if answer['Choice'] == 'Most likely to form carbamate for each resid':
            all_pka_sasa_res = pd.read_csv('Output/results.csv')
            all_pka_sasa_res = all_pka_sasa_res.drop(['Unnamed: 0'], axis=1)
            #carbam_pdb_list = ['/1U8F-alt1A_assembled/1U8F-alt1B_assembled', '/2LJ5-alt-78_assembled/2LJ5-alt-100_assembled/2LJ5-alt-114_assembled/2LJ5-alt-44_assembled/2LJ5-alt-50_assembled/2LJ5-alt-128_assembled/2LJ5-alt-87_assembled/2LJ5-alt-93_assembled/2LJ5-alt-276_assembled/2LJ5-alt-262_assembled/2LJ5-alt-289_assembled/2LJ5-alt-288_assembled/2LJ5-alt-263_assembled/2LJ5-alt-277_assembled/2LJ5-alt-92_assembled/2LJ5-alt-86_assembled/2LJ5-alt-129_assembled/2LJ5-alt-51_assembled/2LJ5-alt-45_assembled/2LJ5-alt-115_assembled/2LJ5-alt-101_assembled/2LJ5-alt-79_assembled/2LJ5-alt-117_assembled/2LJ5-alt-103_assembled/2LJ5-alt-53_assembled/2LJ5-alt-47_assembled/2LJ5-alt-90_assembled/2LJ5-alt-84_assembled/2LJ5-alt-301_assembled/2LJ5-alt-261_assembled/2LJ5-alt-275_assembled/2LJ5-alt-249_assembled/2LJ5-alt-248_assembled/2LJ5-alt-274_assembled/2LJ5-alt-260_assembled/2LJ5-alt-300_assembled/2LJ5-alt-85_assembled/2LJ5-alt-91_assembled/2LJ5-alt-46_assembled/2LJ5-alt-52_assembled/2LJ5-alt-102_assembled/2LJ5-alt-116_assembled/2LJ5-alt-56_assembled/2LJ5-alt-42_assembled/2LJ5-alt-112_assembled/2LJ5-alt-106_assembled/2LJ5-alt-95_assembled/2LJ5-alt-81_assembled/2LJ5-alt-258_assembled/2LJ5-alt-264_assembled/2LJ5-alt-270_assembled/2LJ5-alt-271_assembled/2LJ5-alt-265_assembled/2LJ5-alt-259_assembled/2LJ5-alt-80_assembled/2LJ5-alt-94_assembled/2LJ5-alt-107_assembled/2LJ5-alt-113_assembled/2LJ5-alt-43_assembled/2LJ5-alt-57_assembled/2LJ5-alt-139_assembled/2LJ5-alt-41_assembled/2LJ5-alt-55_assembled/2LJ5-alt-105_assembled/2LJ5-alt-111_assembled/2LJ5-alt-69_assembled/2LJ5-alt-82_assembled/2LJ5-alt-96_assembled/2LJ5-alt-273_assembled/2LJ5-alt-267_assembled/2LJ5-alt-298_assembled/2LJ5-alt-299_assembled/2LJ5-alt-266_assembled/2LJ5-alt-272_assembled/2LJ5-alt-97_assembled/2LJ5-alt-83_assembled/2LJ5-alt-68_assembled/2LJ5-alt-110_assembled/2LJ5-alt-104_assembled/2LJ5-alt-54_assembled/2LJ5-alt-40_assembled/2LJ5-alt-138_assembled/2LJ5-alt-163_assembled/2LJ5-alt-177_assembled/2LJ5-alt-27_assembled/2LJ5-alt-33_assembled/2LJ5-alt-188_assembled/2LJ5-alt-215_assembled/2LJ5-alt-201_assembled/2LJ5-alt-229_assembled/2LJ5-alt-3_assembled/2LJ5-alt-2_assembled/2LJ5-alt-228_assembled/2LJ5-alt-200_assembled/2LJ5-alt-214_assembled/2LJ5-alt-189_assembled/2LJ5-alt-32_assembled/2LJ5-alt-26_assembled/2LJ5-alt-176_assembled/2LJ5-alt-162_assembled/2LJ5-alt-174_assembled/2LJ5-alt-160_assembled/2LJ5-alt-18_assembled/2LJ5-alt-148_assembled/2LJ5-alt-30_assembled/2LJ5-alt-24_assembled/2LJ5-alt-202_assembled/2LJ5-alt-216_assembled/2LJ5-alt-1_assembled/2LJ5-alt-217_assembled/2LJ5-alt-203_assembled/2LJ5-alt-25_assembled/2LJ5-alt-31_assembled/2LJ5-alt-149_assembled/2LJ5-alt-19_assembled/2LJ5-alt-161_assembled/2LJ5-alt-175_assembled/2LJ5-alt-35_assembled/2LJ5-alt-21_assembled/2LJ5-alt-159_assembled/2LJ5-alt-171_assembled/2LJ5-alt-165_assembled/2LJ5-alt-207_assembled/2LJ5-alt-213_assembled/2LJ5-alt-5_assembled/2LJ5-alt-4_assembled/2LJ5-alt-212_assembled/2LJ5-alt-206_assembled/2LJ5-alt-164_assembled/2LJ5-alt-170_assembled/2LJ5-alt-158_assembled/2LJ5-alt-20_assembled/2LJ5-alt-34_assembled/2LJ5-alt-22_assembled/2LJ5-alt-36_assembled/2LJ5-alt-166_assembled/2LJ5-alt-172_assembled/2LJ5-alt-199_assembled/2LJ5-alt-238_assembled/2LJ5-alt-210_assembled/2LJ5-alt-204_assembled/2LJ5-alt-6_assembled/2LJ5-alt-7_assembled/2LJ5-alt-205_assembled/2LJ5-alt-211_assembled/2LJ5-alt-239_assembled/2LJ5-alt-198_assembled/2LJ5-alt-173_assembled/2LJ5-alt-167_assembled/2LJ5-alt-37_assembled/2LJ5-alt-23_assembled/2LJ5-alt-142_assembled/2LJ5-alt-156_assembled/2LJ5-alt-12_assembled/2LJ5-alt-181_assembled/2LJ5-alt-195_assembled/2LJ5-alt-234_assembled/2LJ5-alt-220_assembled/2LJ5-alt-208_assembled/2LJ5-alt-209_assembled/2LJ5-alt-221_assembled/2LJ5-alt-235_assembled/2LJ5-alt-194_assembled/2LJ5-alt-180_assembled/2LJ5-alt-13_assembled/2LJ5-alt-157_assembled/2LJ5-alt-143_assembled/2LJ5-alt-155_assembled/2LJ5-alt-39_assembled/2LJ5-alt-141_assembled/2LJ5-alt-11_assembled/2LJ5-alt-169_assembled/2LJ5-alt-196_assembled/2LJ5-alt-182_assembled/2LJ5-alt-223_assembled/2LJ5-alt-237_assembled/2LJ5-alt-9_assembled/2LJ5-alt-8_assembled/2LJ5-alt-236_assembled/2LJ5-alt-222_assembled/2LJ5-alt-183_assembled/2LJ5-alt-197_assembled/2LJ5-alt-168_assembled/2LJ5-alt-10_assembled/2LJ5-alt-140_assembled/2LJ5-alt-38_assembled/2LJ5-alt-154_assembled/2LJ5-alt-14_assembled/2LJ5-alt-178_assembled/2LJ5-alt-150_assembled/2LJ5-alt-28_assembled/2LJ5-alt-144_assembled/2LJ5-alt-193_assembled/2LJ5-alt-187_assembled/2LJ5-alt-226_assembled/2LJ5-alt-232_assembled/2LJ5-alt-233_assembled/2LJ5-alt-227_assembled/2LJ5-alt-186_assembled/2LJ5-alt-192_assembled/2LJ5-alt-145_assembled/2LJ5-alt-29_assembled/2LJ5-alt-151_assembled/2LJ5-alt-179_assembled/2LJ5-alt-15_assembled/2LJ5-alt-17_assembled/2LJ5-alt-147_assembled/2LJ5-alt-153_assembled/2LJ5-alt-184_assembled/2LJ5-alt-190_assembled/2LJ5-alt-219_assembled/2LJ5-alt-231_assembled/2LJ5-alt-225_assembled/2LJ5-alt-224_assembled/2LJ5-alt-230_assembled/2LJ5-alt-218_assembled/2LJ5-alt-191_assembled/2LJ5-alt-185_assembled/2LJ5-alt-152_assembled/2LJ5-alt-146_assembled/2LJ5-alt-16_assembled/2LJ5-alt-121_assembled/2LJ5-alt-59_assembled/2LJ5-alt-135_assembled/2LJ5-alt-65_assembled/2LJ5-alt-109_assembled/2LJ5-alt-71_assembled/2LJ5-alt-257_assembled/2LJ5-alt-243_assembled/2LJ5-alt-294_assembled/2LJ5-alt-280_assembled/2LJ5-alt-281_assembled/2LJ5-alt-295_assembled/2LJ5-alt-242_assembled/2LJ5-alt-256_assembled/2LJ5-alt-70_assembled/2LJ5-alt-108_assembled/2LJ5-alt-64_assembled/2LJ5-alt-134_assembled/2LJ5-alt-58_assembled/2LJ5-alt-120_assembled/2LJ5-alt-136_assembled/2LJ5-alt-122_assembled/2LJ5-alt-72_assembled/2LJ5-alt-66_assembled/2LJ5-alt-99_assembled/2LJ5-alt-240_assembled/2LJ5-alt-254_assembled/2LJ5-alt-268_assembled/2LJ5-alt-283_assembled/2LJ5-alt-297_assembled/2LJ5-alt-296_assembled/2LJ5-alt-282_assembled/2LJ5-alt-269_assembled/2LJ5-alt-255_assembled/2LJ5-alt-241_assembled/2LJ5-alt-98_assembled/2LJ5-alt-67_assembled/2LJ5-alt-73_assembled/2LJ5-alt-123_assembled/2LJ5-alt-137_assembled/2LJ5-alt-77_assembled/2LJ5-alt-63_assembled/2LJ5-alt-133_assembled/2LJ5-alt-127_assembled/2LJ5-alt-88_assembled/2LJ5-alt-279_assembled/2LJ5-alt-245_assembled/2LJ5-alt-251_assembled/2LJ5-alt-286_assembled/2LJ5-alt-292_assembled/2LJ5-alt-293_assembled/2LJ5-alt-287_assembled/2LJ5-alt-250_assembled/2LJ5-alt-244_assembled/2LJ5-alt-278_assembled/2LJ5-alt-89_assembled/2LJ5-alt-126_assembled/2LJ5-alt-132_assembled/2LJ5-alt-62_assembled/2LJ5-alt-76_assembled/2LJ5-alt-60_assembled/2LJ5-alt-118_assembled/2LJ5-alt-74_assembled/2LJ5-alt-124_assembled/2LJ5-alt-48_assembled/2LJ5-alt-130_assembled/2LJ5-alt-252_assembled/2LJ5-alt-246_assembled/2LJ5-alt-291_assembled/2LJ5-alt-285_assembled/2LJ5-alt-284_assembled/2LJ5-alt-290_assembled/2LJ5-alt-247_assembled/2LJ5-alt-253_assembled/2LJ5-alt-131_assembled/2LJ5-alt-49_assembled/2LJ5-alt-125_assembled/2LJ5-alt-75_assembled/2LJ5-alt-119_assembled/2LJ5-alt-61_assembled/2MBH-alt-20_assembled/2MBH-alt-19_assembled/2MBH-alt-18_assembled/2MBH-alt-9_assembled/2MBH-alt-8_assembled/2MBH-alt-1_assembled/2MBH-alt-3_assembled/2MBH-alt-2_assembled/2MBH-alt-6_assembled/2MBH-alt-7_assembled/2MBH-alt-5_assembled/2MBH-alt-4_assembled/2MBH-alt-15_assembled/2MBH-alt-14_assembled/2MBH-alt-16_assembled/2MBH-alt-17_assembled/2MBH-alt-13_assembled/2MBH-alt-12_assembled/2MBH-alt-10_assembled/2MBH-alt-11_assembled','/AF-P0DMV9-F1-model_v1','/AF-P0DMV8-F1-model_v1','/AF-P10809-F1-model_v1','/AF-P10809-F1-model_v1','/AF-P63261-F1-model_v1','/AF-Q99832-F1-model_v1','/AF-P62807-F1-model_v1','/AF-P62807-F1-model_v1','/3BYH-alt-1_assembled','/3B6R-alt-1_assembled/3DRB-alt-1_assembled','/5W8J-alt1B_assembled/5W8J-alt1A_assembled','/6SBV-alt1B_assembled/6SBV-alt1A_assembled','/AF-Q12931-F1-model_v1','/AF-Q8N257-F1-model_v1','/AF-Q6PEY2-F1-model_v1','/AF-P68032-F1-model_v1','/AF-P63267-F1-model_v1','/AF-O60506-F1-model_v1','/AF-P62750-F1-model_v1','/1NI2-alt-1_assembled','/1QCK-alt-3_assembled/1QCK-alt-10_assembled/1QCK-alt-2_assembled/1QCK-alt-12_assembled/1QCK-alt-13_assembled/1QCK-alt-1_assembled/1QCK-alt-5_assembled/1QCK-alt-17_assembled/1QCK-alt-16_assembled/1QCK-alt-4_assembled/1QCK-alt-6_assembled/1QCK-alt-14_assembled/1QCK-alt-7_assembled/1QCK-alt-18_assembled/1QCK-alt-19_assembled/1QCK-alt-9_assembled/1QCK-alt-8_assembled/1QCK-alt-21_assembled/1QCK-alt-20_assembled','/AF-P60981-F1-model_v1','/AF-Q562R1-F1-model_v1','/AF-P62736-F1-model_v1','/AF-P48668-F1-model_v1','/3TZD-alt1C_assembled/3TZD-alt1B_assembled/3TZD-alt1A_assembled','/AF-P02538-F1-model_v1','/AF-O00567-F1-model_v1','/4XZ2-alt-1_assembled/4XYJ-alt-1_assembled','/AF-Q9ULV4-F1-model_v1','/1PL8-alt1B_assembled/1PL8-alt1A_assembled','/AF-Q53GS9-F1-model_v1','/AF-Q9Y2R4-F1-model_v1','/AF-P55265-F1-model_v1','/AF-Q6S8J3-F1-model_v1','/2P5X-alt-1_assembled/AF-O95671-F1-model_v1','/AF-Q8N142-F1-model_v1','/AF-O00471-F1-model_v1','/AF-Q68EM7-F1-model_v1','/AF-P37802-F1-model_v1','/AF-P35527-F1-model_v1','/AF-P07900-F1-model_v1','/AF-P61224-F1-model_v1','/5EXW-alt1B_assembled/5EXW-alt1A_assembled/5F0X-alt1B_assembled/5F0X-alt1A_assembled','/AF-P08133-F1-model_v1','/AF-P39023-F1-model_v1','/AF-O43390-F1-model_v1','/AF-P13861-F1-model_v1','/AF-Q9Y285-F1-model_v1','/6ZJE-alt1A_assembled/6ZJE-alt1C_assembled/6ZJE-alt1B_assembled','/5DZZ-alt-1_assembled','/AF-Q969X6-F1-model_v1','/AF-Q63ZY6-F1-model_v1','/3BBF-alt1B_assembled/3BBF-alt1A_assembled/3BBC-alt1B_assembled/3BBC-alt1A_assembled','/1UCN-alt-1_assembled/2HVD-alt-1_assembled','/AF-P62805-F1-model_v1','/6Y1E-alt1A_assembled/6Y1E-alt1B_assembled/2A2R-alt1B_assembled/2A2R-alt1A_assembled','/4XBJ-alt1B_assembled/4XBJ-alt1A_assembled/4WR3-alt1B_assembled/4WR3-alt1A_assembled','/1BXN-alt-1_assembled','/1BD0-alt-1_assembled/1SFT-alt-1_assembled','/1E3U-alt1B_assembled/1E3U-alt1A_assembled/1EWZ-alt-1_assembled','/1E9Y-alt-1_assembled/1E9Z-alt-1_assembled','/2OEJ-alt-1_assembled','/1RQE-alt1B_assembled/1RQE-alt1A_assembled','/AF-P05020-F1-model_v1','/7B53-alt1A_assembled/7B53-alt1B_assembled/7B60-alt1A_assembled/7B60-alt1B_assembled','/1YBQ-alt1B_assembled/1YBQ-alt1A_assembled','/3BG9-alt-1_assembled','/4C6B-alt1A_assembled/4C6B-alt1B_assembled']
            #carbam_resid_list = ['66', '11', '539', '539', '218', '233', '328', '401', '6', '58', '328', '156', '41', '42', '577', '58', '394', '330', '329', '256', '88', '83', '18', '44', '329', '330', '426', '34', '426', '230', '156', '205', '5', '459', '583', '564', '295', '100', '158', '297', '295', '88', '321', '419', '104', '125', '191', '39', '374', '135', '117', '146', '2026', '281', '56', '124', '12', '92', '127', '122', '203', '129', '70', '219', '173', '184', '103', '225', '162', '741', '1556']
            get_most_likely_value(all_pka_sasa_res)
        
        if answer['Choice'] == 'No processing needed':
            all_pka_sasa_res = pd.read_csv('Output/results_2608/results_carbamates_marked.csv')
            all_pka_sasa_res = all_pka_sasa_res.drop(['Unnamed: 0'], axis=1)
            carbam_pdb_list = ['4XBJ-alt1B', '4XBJ-alt1A']
            carbam_resid_list = ['122', '122']
            analyse_data(all_pka_sasa_res, carbam_pdb_list, carbam_resid_list)
        
        if answer['Choice'] == 'Plot input csv file':

            all_pka_sasa_res = pd.read_csv('Output/results_2708_marked.csv')
            all_pka_sasa_res = all_pka_sasa_res.drop(['Unnamed: 0'], axis=1)
            fig, ax = plt.subplots()
            colors = {True:'#68246D', False:'black'}
            alphas = {True:1, False:0.05}
            all_pka_sasa_res['sasa'] = all_pka_sasa_res['sasa'].astype(float)
            all_pka_sasa_res['pKa'] = all_pka_sasa_res['pKa'].astype(float)
            plt.scatter(all_pka_sasa_res['sasa'], all_pka_sasa_res['pKa'], c=all_pka_sasa_res['carbamylated'].map(colors), alpha=all_pka_sasa_res['carbamylated'].map(alphas))
            plt.title('pKa vs sasa')
            plt.xlabel('sasa')
            plt.ylabel('pKa')
            plt.show()
        
        if answer['Choice'] == 'Only plot carbamates':
            all_pka_sasa_res = pd.read_csv('Output/results_2608/results_avg_carbam_marked.csv')
            all_pka_sasa_res = all_pka_sasa_res.drop(['Unnamed: 0'], axis=1)
            df_only_carbs = all_pka_sasa_res.where(all_pka_sasa_res['carbamylated'] == True)
            df_only_carbs = df_only_carbs[df_only_carbs['carbamylated'].notna()]
            df_only_carbs['AF'] = df_only_carbs['PDB Code'].str.startswith('/AF')
            print(df_only_carbs)
            colors = {True:'red', False:'black'}
            plt.scatter(df_only_carbs['sasa'], df_only_carbs['pKa'], c=df_only_carbs['AF'].map(colors))
            plt.title('Graph Showing Which Points Come From Alphafold Structures')
            plt.xlabel('sasa')
            plt.ylabel('pKa')
            plt.show()

        if answer['Choice'] == 'Plot with errors':
            all_pka_sasa_res = pd.read_csv('Output/results_2608/results_avg_carbam_marked.csv')
            all_pka_sasa_res = all_pka_sasa_res.drop(['Unnamed: 0'], axis=1)
            df_only_carbs = all_pka_sasa_res.where(all_pka_sasa_res['carbamylated'] == True)
            df_only_carbs = df_only_carbs[df_only_carbs['carbamylated'].notna()]
            plt.errorbar(df_only_carbs['sasa'], df_only_carbs['pKa'], yerr=df_only_carbs['pKa stdev'], xerr=df_only_carbs['sasa stdev'], fmt='o', ecolor='black', mfc='red', mec = 'red', lw = 0.5)
            plt.title('Graph Showing Standard Deviations')
            plt.xlabel('sasa')
            plt.ylabel('pKa')
            plt.show()

        if answer['Choice'] == 'Only PDBs':
            all_pka_sasa_res = pd.read_csv('Output/results_2608/results_carbamates_marked.csv')
            all_pka_sasa_res['AF'] = all_pka_sasa_res['PDB Code'].str.startswith('/AF')
            df_only_pdbs = all_pka_sasa_res.where(all_pka_sasa_res['AF'] == False)
            df_only_pdbs = df_only_pdbs[df_only_pdbs['carbamylated'].notna()]
            colors = {True:'#68246D', False:'black'}
            alphas = {True:1, False:0.05}
            plt.scatter(df_only_pdbs['sasa'], df_only_pdbs['pKa'], c=df_only_pdbs['carbamylated'].map(colors), alpha=df_only_pdbs['carbamylated'].map(alphas))
            plt.show()
            print(df_only_pdbs)

        if answer['Choice'] == 'Highest sasa lowest pka':
            all_pka_sasa_res = pd.read_csv('Output/results_2608/results.csv')
            all_pka_sasa_res = all_pka_sasa_res.drop(['Unnamed: 0'], axis=1)
            low_pka = highest_sasa_lowest_pka(all_pka_sasa_res)

            

        
    except Exception as e:
        print("ERROR: %s"%e)
carbam_pdb_list = ['/1U8F-alt1A_assembled/1U8F-alt1B_assembled', 
'/2LJ5-alt-78_assembled/2LJ5-alt-100_assembled/2LJ5-alt-114_assembled/2LJ5-alt-44_assembled/2LJ5-alt-50_assembled/2LJ5-alt-128_assembled/2LJ5-alt-87_assembled/2LJ5-alt-93_assembled/2LJ5-alt-276_assembled/2LJ5-alt-262_assembled/2LJ5-alt-289_assembled/2LJ5-alt-288_assembled/2LJ5-alt-263_assembled/2LJ5-alt-277_assembled/2LJ5-alt-92_assembled/2LJ5-alt-86_assembled/2LJ5-alt-129_assembled/2LJ5-alt-51_assembled/2LJ5-alt-45_assembled/2LJ5-alt-115_assembled/2LJ5-alt-101_assembled/2LJ5-alt-79_assembled/2LJ5-alt-117_assembled/2LJ5-alt-103_assembled/2LJ5-alt-53_assembled/2LJ5-alt-47_assembled/2LJ5-alt-90_assembled/2LJ5-alt-84_assembled/2LJ5-alt-301_assembled/2LJ5-alt-261_assembled/2LJ5-alt-275_assembled/2LJ5-alt-249_assembled/2LJ5-alt-248_assembled/2LJ5-alt-274_assembled/2LJ5-alt-260_assembled/2LJ5-alt-300_assembled/2LJ5-alt-85_assembled/2LJ5-alt-91_assembled/2LJ5-alt-46_assembled/2LJ5-alt-52_assembled/2LJ5-alt-102_assembled/2LJ5-alt-116_assembled/2LJ5-alt-56_assembled/2LJ5-alt-42_assembled/2LJ5-alt-112_assembled/2LJ5-alt-106_assembled/2LJ5-alt-95_assembled/2LJ5-alt-81_assembled/2LJ5-alt-258_assembled/2LJ5-alt-264_assembled/2LJ5-alt-270_assembled/2LJ5-alt-271_assembled/2LJ5-alt-265_assembled/2LJ5-alt-259_assembled/2LJ5-alt-80_assembled/2LJ5-alt-94_assembled/2LJ5-alt-107_assembled/2LJ5-alt-113_assembled/2LJ5-alt-43_assembled/2LJ5-alt-57_assembled/2LJ5-alt-139_assembled/2LJ5-alt-41_assembled/2LJ5-alt-55_assembled/2LJ5-alt-105_assembled/2LJ5-alt-111_assembled/2LJ5-alt-69_assembled/2LJ5-alt-82_assembled/2LJ5-alt-96_assembled/2LJ5-alt-273_assembled/2LJ5-alt-267_assembled/2LJ5-alt-298_assembled/2LJ5-alt-299_assembled/2LJ5-alt-266_assembled/2LJ5-alt-272_assembled/2LJ5-alt-97_assembled/2LJ5-alt-83_assembled/2LJ5-alt-68_assembled/2LJ5-alt-110_assembled/2LJ5-alt-104_assembled/2LJ5-alt-54_assembled/2LJ5-alt-40_assembled/2LJ5-alt-138_assembled/2LJ5-alt-163_assembled/2LJ5-alt-177_assembled/2LJ5-alt-27_assembled/2LJ5-alt-33_assembled/2LJ5-alt-188_assembled/2LJ5-alt-215_assembled/2LJ5-alt-201_assembled/2LJ5-alt-229_assembled/2LJ5-alt-3_assembled/2LJ5-alt-2_assembled/2LJ5-alt-228_assembled/2LJ5-alt-200_assembled/2LJ5-alt-214_assembled/2LJ5-alt-189_assembled/2LJ5-alt-32_assembled/2LJ5-alt-26_assembled/2LJ5-alt-176_assembled/2LJ5-alt-162_assembled/2LJ5-alt-174_assembled/2LJ5-alt-160_assembled/2LJ5-alt-18_assembled/2LJ5-alt-148_assembled/2LJ5-alt-30_assembled/2LJ5-alt-24_assembled/2LJ5-alt-202_assembled/2LJ5-alt-216_assembled/2LJ5-alt-1_assembled/2LJ5-alt-217_assembled/2LJ5-alt-203_assembled/2LJ5-alt-25_assembled/2LJ5-alt-31_assembled/2LJ5-alt-149_assembled/2LJ5-alt-19_assembled/2LJ5-alt-161_assembled/2LJ5-alt-175_assembled/2LJ5-alt-35_assembled/2LJ5-alt-21_assembled/2LJ5-alt-159_assembled/2LJ5-alt-171_assembled/2LJ5-alt-165_assembled/2LJ5-alt-207_assembled/2LJ5-alt-213_assembled/2LJ5-alt-5_assembled/2LJ5-alt-4_assembled/2LJ5-alt-212_assembled/2LJ5-alt-206_assembled/2LJ5-alt-164_assembled/2LJ5-alt-170_assembled/2LJ5-alt-158_assembled/2LJ5-alt-20_assembled/2LJ5-alt-34_assembled/2LJ5-alt-22_assembled/2LJ5-alt-36_assembled/2LJ5-alt-166_assembled/2LJ5-alt-172_assembled/2LJ5-alt-199_assembled/2LJ5-alt-238_assembled/2LJ5-alt-210_assembled/2LJ5-alt-204_assembled/2LJ5-alt-6_assembled/2LJ5-alt-7_assembled/2LJ5-alt-205_assembled/2LJ5-alt-211_assembled/2LJ5-alt-239_assembled/2LJ5-alt-198_assembled/2LJ5-alt-173_assembled/2LJ5-alt-167_assembled/2LJ5-alt-37_assembled/2LJ5-alt-23_assembled/2LJ5-alt-142_assembled/2LJ5-alt-156_assembled/2LJ5-alt-12_assembled/2LJ5-alt-181_assembled/2LJ5-alt-195_assembled/2LJ5-alt-234_assembled/2LJ5-alt-220_assembled/2LJ5-alt-208_assembled/2LJ5-alt-209_assembled/2LJ5-alt-221_assembled/2LJ5-alt-235_assembled/2LJ5-alt-194_assembled/2LJ5-alt-180_assembled/2LJ5-alt-13_assembled/2LJ5-alt-157_assembled/2LJ5-alt-143_assembled/2LJ5-alt-155_assembled/2LJ5-alt-39_assembled/2LJ5-alt-141_assembled/2LJ5-alt-11_assembled/2LJ5-alt-169_assembled/2LJ5-alt-196_assembled/2LJ5-alt-182_assembled/2LJ5-alt-223_assembled/2LJ5-alt-237_assembled/2LJ5-alt-9_assembled/2LJ5-alt-8_assembled/2LJ5-alt-236_assembled/2LJ5-alt-222_assembled/2LJ5-alt-183_assembled/2LJ5-alt-197_assembled/2LJ5-alt-168_assembled/2LJ5-alt-10_assembled/2LJ5-alt-140_assembled/2LJ5-alt-38_assembled/2LJ5-alt-154_assembled/2LJ5-alt-14_assembled/2LJ5-alt-178_assembled/2LJ5-alt-150_assembled/2LJ5-alt-28_assembled/2LJ5-alt-144_assembled/2LJ5-alt-193_assembled/2LJ5-alt-187_assembled/2LJ5-alt-226_assembled/2LJ5-alt-232_assembled/2LJ5-alt-233_assembled/2LJ5-alt-227_assembled/2LJ5-alt-186_assembled/2LJ5-alt-192_assembled/2LJ5-alt-145_assembled/2LJ5-alt-29_assembled/2LJ5-alt-151_assembled/2LJ5-alt-179_assembled/2LJ5-alt-15_assembled/2LJ5-alt-17_assembled/2LJ5-alt-147_assembled/2LJ5-alt-153_assembled/2LJ5-alt-184_assembled/2LJ5-alt-190_assembled/2LJ5-alt-219_assembled/2LJ5-alt-231_assembled/2LJ5-alt-225_assembled/2LJ5-alt-224_assembled/2LJ5-alt-230_assembled/2LJ5-alt-218_assembled/2LJ5-alt-191_assembled/2LJ5-alt-185_assembled/2LJ5-alt-152_assembled/2LJ5-alt-146_assembled/2LJ5-alt-16_assembled/2LJ5-alt-121_assembled/2LJ5-alt-59_assembled/2LJ5-alt-135_assembled/2LJ5-alt-65_assembled/2LJ5-alt-109_assembled/2LJ5-alt-71_assembled/2LJ5-alt-257_assembled/2LJ5-alt-243_assembled/2LJ5-alt-294_assembled/2LJ5-alt-280_assembled/2LJ5-alt-281_assembled/2LJ5-alt-295_assembled/2LJ5-alt-242_assembled/2LJ5-alt-256_assembled/2LJ5-alt-70_assembled/2LJ5-alt-108_assembled/2LJ5-alt-64_assembled/2LJ5-alt-134_assembled/2LJ5-alt-58_assembled/2LJ5-alt-120_assembled/2LJ5-alt-136_assembled/2LJ5-alt-122_assembled/2LJ5-alt-72_assembled/2LJ5-alt-66_assembled/2LJ5-alt-99_assembled/2LJ5-alt-240_assembled/2LJ5-alt-254_assembled/2LJ5-alt-268_assembled/2LJ5-alt-283_assembled/2LJ5-alt-297_assembled/2LJ5-alt-296_assembled/2LJ5-alt-282_assembled/2LJ5-alt-269_assembled/2LJ5-alt-255_assembled/2LJ5-alt-241_assembled/2LJ5-alt-98_assembled/2LJ5-alt-67_assembled/2LJ5-alt-73_assembled/2LJ5-alt-123_assembled/2LJ5-alt-137_assembled/2LJ5-alt-77_assembled/2LJ5-alt-63_assembled/2LJ5-alt-133_assembled/2LJ5-alt-127_assembled/2LJ5-alt-88_assembled/2LJ5-alt-279_assembled/2LJ5-alt-245_assembled/2LJ5-alt-251_assembled/2LJ5-alt-286_assembled/2LJ5-alt-292_assembled/2LJ5-alt-293_assembled/2LJ5-alt-287_assembled/2LJ5-alt-250_assembled/2LJ5-alt-244_assembled/2LJ5-alt-278_assembled/2LJ5-alt-89_assembled/2LJ5-alt-126_assembled/2LJ5-alt-132_assembled/2LJ5-alt-62_assembled/2LJ5-alt-76_assembled/2LJ5-alt-60_assembled/2LJ5-alt-118_assembled/2LJ5-alt-74_assembled/2LJ5-alt-124_assembled/2LJ5-alt-48_assembled/2LJ5-alt-130_assembled/2LJ5-alt-252_assembled/2LJ5-alt-246_assembled/2LJ5-alt-291_assembled/2LJ5-alt-285_assembled/2LJ5-alt-284_assembled/2LJ5-alt-290_assembled/2LJ5-alt-247_assembled/2LJ5-alt-253_assembled/2LJ5-alt-131_assembled/2LJ5-alt-49_assembled/2LJ5-alt-125_assembled/2LJ5-alt-75_assembled/2LJ5-alt-119_assembled/2LJ5-alt-61_assembled/2MBH-alt-20_assembled/2MBH-alt-19_assembled/2MBH-alt-18_assembled/2MBH-alt-9_assembled/2MBH-alt-8_assembled/2MBH-alt-1_assembled/2MBH-alt-3_assembled/2MBH-alt-2_assembled/2MBH-alt-6_assembled/2MBH-alt-7_assembled/2MBH-alt-5_assembled/2MBH-alt-4_assembled/2MBH-alt-15_assembled/2MBH-alt-14_assembled/2MBH-alt-16_assembled/2MBH-alt-17_assembled/2MBH-alt-13_assembled/2MBH-alt-12_assembled/2MBH-alt-10_assembled/2MBH-alt-11_assembled',
'/AF-P0DMV9-F1-model_v1',
'/AF-P0DMV8-F1-model_v1',
'/AF-P10809-F1-model_v1',
'/AF-P10809-F1-model_v1',
'/AF-P63261-F1-model_v1',
'/AF-Q99832-F1-model_v1',
'/AF-P62807-F1-model_v1',
'/AF-P62807-F1-model_v1',
'/3BYH-alt-1_assembled',
'/3B6R-alt-1_assembled/3DRB-alt-1_assembled',
'/5W8J-alt1B_assembled/5W8J-alt1A_assembled',
'/6SBV-alt1B_assembled/6SBV-alt1A_assembled',
'/AF-Q12931-F1-model_v1'
'/AF-Q8N257-F1-model_v1',
'/AF-Q6PEY2-F1-model_v1',
'/AF-P68032-F1-model_v1',
'/AF-P63267-F1-model_v1',
'/AF-O60506-F1-model_v1',
'/AF-P62750-F1-model_v1',
'/1NI2-alt-1_assembled',
'/1QCK-alt-3_assembled/1QCK-alt-10_assembled/1QCK-alt-2_assembled/1QCK-alt-12_assembled/1QCK-alt-13_assembled/1QCK-alt-1_assembled/1QCK-alt-5_assembled/1QCK-alt-17_assembled/1QCK-alt-16_assembled/1QCK-alt-4_assembled/1QCK-alt-6_assembled/1QCK-alt-14_assembled/1QCK-alt-7_assembled/1QCK-alt-18_assembled/1QCK-alt-19_assembled/1QCK-alt-9_assembled/1QCK-alt-8_assembled/1QCK-alt-21_assembled/1QCK-alt-20_assembled',
'/AF-P60981-F1-model_v1',
'/AF-Q562R1-F1-model_v1',
'/AF-P62736-F1-model_v1',
'/AF-P48668-F1-model_v1',
'/3TZD-alt1C_assembled/3TZD-alt1B_assembled/3TZD-alt1A_assembled',
'/AF-P02538-F1-model_v1',
'/AF-O00567-F1-model_v1',
'/4XZ2-alt-1_assembled/4XYJ-alt-1_assembled',
'/AF-Q9ULV4-F1-model_v1',
'/1PL8-alt1B_assembled/1PL8-alt1A_assembled',
'/AF-Q53GS9-F1-model_v1',
'/AF-Q9Y2R4-F1-model_v1',
'/AF-P55265-F1-model_v1',
'/AF-Q6S8J3-F1-model_v1',
'/2P5X-alt-1_assembled/AF-O95671-F1-model_v1',
'/AF-Q8N142-F1-model_v1',
'/AF-O00471-F1-model_v1',
'/AF-Q68EM7-F1-model_v1',
'/AF-P37802-F1-model_v1',
'/AF-P35527-F1-model_v1',
'/AF-P07900-F1-model_v1',
'/AF-P61224-F1-model_v1',
'/5EXW-alt1B_assembled/5EXW-alt1A_assembled/5F0X-alt1B_assembled/5F0X-alt1A_assembled',
'/AF-P08133-F1-model_v1',
'/AF-P39023-F1-model_v1',
'/AF-O43390-F1-model_v1',
'/AF-P13861-F1-model_v1',
'/AF-Q9Y285-F1-model_v1',
'/6ZJE-alt1A_assembled/6ZJE-alt1C_assembled/6ZJE-alt1B_assembled',
'/5DZZ-alt-1_assembled',
'/AF-Q969X6-F1-model_v1',
'/AF-Q63ZY6-F1-model_v1',
'/3BBF-alt1B_assembled/3BBF-alt1A_assembled/3BBC-alt1B_assembled/3BBC-alt1A_assembled',
'/1UCN-alt-1_assembled/2HVD-alt-1_assembled',
'/AF-P62805-F1-model_v1',
'/6Y1E-alt1A_assembled/6Y1E-alt1B_assembled/2A2R-alt1B_assembled/2A2R-alt1A_assembled',
'/4XBJ-alt1B_assembled/4XBJ-alt1A_assembled/4WR3-alt1B_assembled/4WR3-alt1A_assembled',
'/1BXN-alt-1_assembled',
'/1BD0-alt-1_assembled/1SFT-alt-1_assembled',
'/1E3U-alt1B_assembled/1E3U-alt1A_assembled/1EWZ-alt-1_assembled',
'/1E9Y-alt-1_assembled/1E9Z-alt-1_assembled',
'/2OEJ-alt-1_assembled',
'/1RQE-alt1B_assembled/1RQE-alt1A_assembled',
'/AF-P05020-F1-model_v1',
'/7B53-alt1A_assembled/7B53-alt1B_assembled/7B60-alt1A_assembled/7B60-alt1B_assembled',
'/1YBQ-alt1B_assembled/1YBQ-alt1A_assembled',
'/3BG9-alt-1_assembled',
'/4C6B-alt1A_assembled/4C6B-alt1B_assembled']

carbam_resid_list = ['66', '11', '539', '539', '218', '233', '328', '401', '6', '58', '328', '156', '41', '42', '577', '58', '394', '330', '329', '256', '88', '83', '18', '44', '329', '330', '426', '34', '426', '230', '156', '205', '5', '459', '583', '564', '295', '100', '158', '297', 
'295', '88', '321', '419', '104', '125', '191', '39', '374', '135', '117', '146', '2026', '281', '56', '124', '12', '92', '127', '122', '203', '129', '70', '219', '173', '184', '103', '225', '162', '741', '1556']