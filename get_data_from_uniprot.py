import inquirer
import urllib.request, urllib.parse, urllib.error
import re
from bs4 import BeautifulSoup
import csv
import pandas as pd
import os
import subprocess
import numpy as np
import sys
import biobox as bb



#This function obtains all the pdb codes given an organism
#Note the name has to be exactly that used on the uniprot website and the code needs to be the code in the URL for the proteome
def get_pdbs(name_of_organism, code):

    #Firstly a df is constructed for our results to go in

    columns = ['Uniprot Entry', 'PDB Code', 'Method Structure Obtained by', 'Resolution', 'Chains']
    df = pd.DataFrame(columns=columns)

    list_of_entries = list()
    list_of_pdbs = list()
    list_clean = list()
    starting_number = 0
    end = 0
    list_UNIPROT_codes = list()
    try:
#Next it finds the names of all the chromosomes the organism has (this information is needed later)

        web_url = 'https://www.uniprot.org/proteomes/' + code
        html = urllib.request.urlopen(web_url)
        soup = BeautifulSoup(html, 'html.parser')
        table = soup.find_all('table')
        for line in table:
            line = str(line)
            ID = re.findall('name="[\w ]*" type', line)

            IDS = list()

            for chromosome in ID:
                chromosome = chromosome[6:-6]
                chromosome = chromosome.replace(' ', '+')
                IDS.append(chromosome)

    except Exception as e:
        print("ERROR: %s"%e)
        print('Uniprot code is invalid')
        return()

    for chromosome in IDS:
        try:
#Next the program obtains the other uniprot code which is needed in the URL later on.

            print(chromosome)
            chromosome = chromosome.lower()
            web_url = 'https://www.uniprot.org/uniprot/?query=proteome:' + code + '+AND+proteomecomponent:%22' + chromosome + '%22&sort=score'
            html = html = urllib.request.urlopen(web_url)
            soup = BeautifulSoup(html, 'html.parser')

            for line in soup:
                line = str(line)
                code_other = re.findall('%[\w ()\d -]*%5d%22', line)

            for other_uniprot_code in code_other:
                other_uniprot_code = (code_other[0])[1:-6]

#Next the code gets the species name for the first three entries and checks they correspond to the species name entered.
#This check is important as if incorrect data is added by the user the code will produce incorrect information.

            if len(IDS) > 1:
                web_url = 'https://www.uniprot.org/uniprot/?query=proteomecomponent%3a%22' + chromosome + '%22&fil=organism%3a%22' + name_of_organism + '+%' + other_uniprot_code + '%5d%22+AND+proteome%3a' + code + '&offset=0&sort=score&columns=id%2centry+name%2creviewed%2cprotein+names%2cgenes%2corganism%2clength'
                html = urllib.request.urlopen(web_url)
                soup = BeautifulSoup(html, 'html.parser')

            if len(IDS) == 1:
                web_url = 'https://www.uniprot.org/uniprot/?query=proteomecomponent%3a' + chromosome + '&fil=organism%3a%22' + name_of_organism + '+%' + other_uniprot_code + '%5d%22+AND+proteome%3a' + code + '&offset=0&sort=score&columns=id%2centry+name%2creviewed%2cprotein+names%2cgenes%2corganism%2clength'
                html = urllib.request.urlopen(web_url)
                soup = BeautifulSoup(html, 'html.parser')

            for line in soup:
                line = str(line)
                list_of_species = re.findall('>[\w ()\d -]*</a></td><td class="number', line)
                for species in list_of_species:
                    species = species[1:-26]
                    species = species.replace(' ', '+')
                    list_clean.append(species)


            if list_clean[0] != name_of_organism:
                end = 1
                raise Exception('Likely due to incorrect information added')
            if list_clean[1] != name_of_organism:
                end = 1
                raise Exception('Likely due to incorrect information added')
            if list_clean[2] != name_of_organism:
                end = 1
                raise Exception('Likely due to incorrect information added')

#This code identifies the number of uniprot codes in the chromosome currently being scanned.

            for line in soup:
                line = str(line)
                number_of_prot = re.findall('var resultsize = .*;', line)
            number_of_prot_2 = number_of_prot[0]
            number_of_prot_2 = int(number_of_prot_2[17:-1])
            print('Number of Uniprot Entries Identified: ', number_of_prot_2)
        

        except Exception as e:
            print("ERROR: %s"%e)
            print('Failed obtaining Uniprot codes for ' + chromosome)
            if end == 1:
                break
            else:
                continue

#This code produces a list going up in increments of 25 (0, 25, 50... [number of entries]).
#This list is needed so for the URL later on so that the code cycles through the uniprot website 25 proteins at a time until the end.

        while starting_number < number_of_prot_2:

            list_of_entries.append(starting_number)
            starting_number = starting_number + 25


        list_of_entries.append(number_of_prot_2)

#Next the code goes through uniprot 25 entries at a time and appends uniprot codes of proteins to a list.
        
        for number in list_of_entries:
            try:
                strnum = str(number)
                web_url = 'https://www.uniprot.org/uniprot/?query=proteomecomponent%3a%22' + chromosome + '%22&fil=organism%3a%22' + name_of_organism + '+%' + other_uniprot_code + '%5d%22+AND+proteome%3a' + code + '&offset=' + strnum + '&sort=score&columns=id%2centry+name%2creviewed%2cprotein+names%2cgenes%2corganism%2clength'
                html = urllib.request.urlopen(web_url)
                soup = BeautifulSoup(html, 'html.parser')

                for line in soup:
                    line = str(line)
                    messy_list_of_uniprot_codes = re.findall('id="\w*"><td', line)

                for messy_protein_code in messy_list_of_uniprot_codes:
                    protein_code = messy_protein_code[4:-5]
                    list_UNIPROT_codes.append(protein_code)

            except Exception as e:
                print("ERROR: %s"%e)
                print('Failed to obtain Uniprot codes on:' + web_url)
                continue
                
            if len(list_UNIPROT_codes) == 0:
                print('No Uniprot codes were found on: ' + web_url)
                continue
#For each uniprot code identified it then parses through the .txt file and appends the Uniprot code, PDB code, method structure obtained by, resolution and chain information to a df for each pdb code.
            for protein_code_clean in list_UNIPROT_codes:

                try:
                    url_2 = 'https://www.uniprot.org/uniprot/' + protein_code_clean + '.txt'
                                            
                    html_2 = urllib.request.urlopen(url_2)

                    for line in html_2:
                        line = str(line)
                        messy_entry = re.findall('PDB; [\w -. ; \d /]*=', line)


                        for entry in messy_entry:
                            try:
                                words = entry.split()
                                chain_ent_num = (len(words) - 1)
                                chain_info = words[chain_ent_num]
                                chain_info = chain_info[:-1]
                                chain_info = chain_info.split('/')
                                PDBCODE = (words[1])[:-1]
                                method_obtained = (words[2])[:-1]
                                resolution = (words[3])[:-1]

#Here the code calls the get_chains function which gets chain information for the protein from the PDB.

                                unique_values_chain = get_chains(PDBCODE)
                                for i in range(len(unique_values_chain)):
                                    data = ({'Uniprot Entry': protein_code_clean, 'PDB Code': PDBCODE, 'Method Structure Obtained by': method_obtained, 'Resolution': resolution, 'Chains': unique_values_chain[i]})
                                    df = df.append(data, ignore_index=True)
                            except Exception as e:
                                print("Error %s"%e)
                                continue

#If for any uniprot code there are no entries it instead appends the alphafold code- if there isn't an alphafold entry it is picked up later but for now it always assumes there is one.

                    if len(messy_entry) == 0:
                        AF_code = 'AF-' + protein_code_clean + 'F1-model_v1'
                        data = ({'Uniprot Entry': protein_code_clean, 'PDB Code': AF_code, 'Method Structure Obtained by': 'Predicted', 'Resolution': 'N/A', 'Chains': 'N/A'})
                        df = df.append(data, ignore_index=True)


                except Exception as e:
                    print("ERROR: %s"%e)
                    print('Failed to obtain PDB codes for ' + protein_code_clean)
                    continue


            number = number + 25

            if number > number_of_prot_2:
                print('Number of Uniprot Entries Searched: ' + str(number_of_prot_2) + '\n')
            else:
                print('Number of Uniprot Entries Searched: ' + str(number) + '\n')
            if number > number_of_prot_2:
                number = number_of_prot_2
            print((str((number/number_of_prot_2)*100))[0:3] + '%')
            

    df.to_csv('temp_res.csv')
    return df




#Given a list of uniprot codes thos function goes to the .txt URL and finds all corresponding .pdb files and all the relevant information.
def get_pdbs_given_uniprot_code(list_UNIPROT_codes):

    columns = ['Uniprot Entry', 'PDB Code', 'Method Structure Obtained by', 'Resolution', 'Chains']
    df = pd.DataFrame(columns=columns)

    for protein_code_clean in list_UNIPROT_codes:
        try:
            url_2 = 'https://www.uniprot.org/uniprot/' + protein_code_clean + '.txt'
            html_2 = urllib.request.urlopen(url_2)
        except Exception as e:
            print('Error %s'%e)
            print('Failed to obtain data for Uniprot entry: ' + protein_code_clean)
            continue

    
        for line in html_2:
            try:
                line = str(line)
                messy_entry = re.findall('PDB; [\w -. ; \d /]*=', line)


                for entry in messy_entry:
                    words = entry.split()

                    PDBCODE = (words[1])[:-1]
                    method_obtained = (words[2])[:-1]
                    resolution = (words[3])[:-1]

                    #As before it calls in the get_chains function to get chain info.

                    unique_values_chain = get_chains(PDBCODE)
                    for i in range(len(unique_values_chain)):
                        data = ({'Uniprot Entry': protein_code_clean, 'PDB Code': PDBCODE, 'Method Structure Obtained by': method_obtained, 'Resolution': resolution, 'Chains': unique_values_chain[i]})
                        df = df.append(data, ignore_index=True)

                        if len(messy_entry) == 0:
                            AF_code = 'AF-' + protein_code_clean + 'F1-model_v1'
                            data = ({'Uniprot Entry': protein_code_clean, 'PDB Code': AF_code, 'Method Structure Obtained by': 'Predicted', 'Resolution': 'N/A', 'Chains': 'N/A'})
                            df = df.append(data, ignore_index=True)

            except Exception as e:
                print('Error %s'%e)
                print('Failed to obtain data for Uniprot entry ' + protein_code_clean)
    print(df)

    return(df)



#Gets chain info given PDB code.
def get_chains(PDBCODE_inpt):
    try:
        try:
            #Fisrt downloads .pdb file and opens in biobox
            subprocess.check_call("wget https://files.rcsb.org/download/" + PDBCODE_inpt + ".pdb", shell=True)
            M = bb.Molecule()
            M.import_pdb(PDBCODE_inpt + '.pdb')
            df = M.data
        except:
            print('LARGE PROTEIN- USING .cif file')
            subprocess.check_call("wget https://files.rcsb.org/download/" + PDBCODE_inpt + ".cif", shell=True)
            M = bb.Molecule()
            M.import_pqr(PDBCODE_inpt + ".cif")
            df = M.data
        
        #Next puts the chains in a lis and makes sure there are no duplicates.

        chain_column = df['chain'].tolist()
        unique_values_chain = list()
        for i in range(len(chain_column)):
            if chain_column[i] not in unique_values_chain:
                unique_values_chain.append(chain_column[i])

        os.remove(PDBCODE_inpt + '.pdb')

    except Exception as e:
        print('Failed to obtain chain information for ' + PDBCODE_inpt)
        print("Error %s"%e)
        print("This may be due to the fact that the protein is too large and therefore a .pdb structure doesn't exist.")
        return()

    return(unique_values_chain)




#Given a list of techniques from the user and the desired resolution this function removes pdb entries from results_df that don't fit the criteeria
def search_by_technique(list_of_techniques, results_df, wanted_res):
    skip = 0
    try:

        results_df = results_df[results_df['Method Structure Obtained by'].isin(list_of_techniques)]

    except Exception as e:
        print("ERROR: %s"%e)
        print('Failed to search by method structure obtained by')
        question = [
        inquirer.List('Choice',
                            message="Continue?",
                            choices=['Yes', 'No'],
                        ),
        ]
        answer = inquirer.prompt(question)

        if answer['Choice'] == 'Yes':
            skip = 0
            pass

        elif answer['Choice'] == 'No':
            skip = 1

    try:
        resolution = float(wanted_res)

        results_df['Resolution'] = results_df['Resolution'].astype(float)

        results_df = results_df[results_df['Resolution']<= resolution]

    except Exception as e:
        print("ERROR: %s"%e)
        print('Failed to search by resolution')
        question = [
        inquirer.List('Choice',
                            message="Continue?",
                            choices=['Yes', 'No'],
                        ),
        ]
        answer = inquirer.prompt(question)

        if answer['Choice'] == 'Yes':
            skip = 0
            pass

        elif answer['Choice'] == 'No':
            skip = 1

    
    return(results_df, skip)






if __name__ == "__main__":

    try:
        #print(get_pdbs('Oryctolagus+cuniculus+(Rabbit)', 'UP000001811'))
        print(get_chains('6YAM'))
        #columns = ['Uniprot Entry', 'PDB Code', 'Method Structure Obtained by', 'Resolution']
        #df = pd.DataFrame(columns=columns)
        #uniprot_code = 'P3892'
        #PDB_entries = '2jfi'
        #technique = 'X-ray'
        #resolution = '2.9'


        #data = ({'Uniprot Entry': uniprot_code, 'PDB Code': PDB_entries, 'Method Structure Obtained by': technique, 'Resolution':resolution})
        #df = df.append(data, ignore_index=True)
        #results_df = df

        #res_raw = (search_by_technique(['X-ray'], results_df, '3.0'))
        #results_df = res_raw[0]
        #skip = res_raw[1]
        #print(results_df)

        #list_UNIPROT_codes = ['P50897']
        #print(get_pdbs_given_uniprot_code(list_UNIPROT_codes))
        
    except Exception as e:
                print("ERROR: %s"%e)


