import glob
import re
import os

import inquirer
import urllib.request, urllib.parse, urllib.error
from bs4 import BeautifulSoup

import pandas as pd
import numpy as np
import biobox as bb

import pdb_loader as pl


def create_empty_dataframe():
    '''
    generate an empty DataFrame according to predefined structure (columns)
    '''
    columns = ['Uniprot Entry', 'PDB Code', 'Method', 'Resolution', 'Chains']
    return pd.DataFrame(columns=columns)
            

def checks(pdb, done_pdbs):
    '''
    Check whether a pdb code is the log_file (the file saying what has already been done).
    '''
    for entry in done_pdbs:

        if re.search(pdb, entry):
            
            print('> %s is already the log_file, skipping...'%pdb)
            if not os.path.exists('log_file.csv'):
                f = open('log_file.csv', 'w')
                f.write('PDB Code,Result\n')
                f.write(pdb + ',Failed as pdb file is already in log_file- remove to continue.')
                f.close()

            elif os.path.exists('log_file.csv'):
                f = open('log_file.csv', 'a')
                f.write('\n')
                f.write(pdb + ',Failed as pdb file is already in log_file- remove to continue.')
                f.close()
                
            return False
          
    return True


def get_organism_proteins(name_of_organism, code, df=[], done_pdbs=[]):
    '''
    Obtain all the PDB codes belonging to an organism.
    The name has to be exactly that used on the uniprot website and the code needs to be the code in the URL for the proteome.    
    '''
    
    #if no DataFrame is provide, build an empty one
    if len(df) == 0:
        df = create_empty_dataframe()


    list_of_entries = list()
    list_clean = list()
    starting_number = 0
    end = 0
    list_UNIPROT_codes = list()
    
    try:
        
        #find the names of all the chromosomes the organism has
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
        raise Exception('Uniprot code is invalid %s'%e)
    
    #iterate over all cromosomes to get their proteins
    for chromosome in IDS:
        
        try:
            
            #obtain the other uniprot code which is needed in the URL later on.
            print("chromosome: %s"%chromosome)
            chromosome = chromosome.lower()
            web_url = 'https://www.uniprot.org/uniprot/?query=proteome:' + code + '+AND+proteomecomponent:%22' + chromosome + '%22&sort=score'
            html = html = urllib.request.urlopen(web_url)
            soup = BeautifulSoup(html, 'html.parser')

            for line in soup:
                line = str(line)
                code_other = re.findall('%[\w ()\d -]*%5d%22', line)

            for other_uniprot_code in code_other:
                other_uniprot_code = (code_other[0])[1:-6]

            #get the species name for the first three entries and checks they correspond to the species name entered.
            #This check is important as if incorrect data is added by the user, the code will produce incorrect information.
            if len(IDS) > 1:
                web_url = 'https://www.uniprot.org/uniprot/?query=proteomecomponent%3a%22' + chromosome + '%22&fil=organism%3a%22' + name_of_organism + '+%' + other_uniprot_code + '%5d%22+AND+proteome%3a' + code + '&offset=0&sort=score&columns=id%2centry+name%2creviewed%2cprotein+names%2cgenes%2corganism%2clength'          
            else:
                #one cromosome only
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

            #identify the number of uniprot codes in the chromosome currently being scanned.
            for line in soup:
                line = str(line)
                number_of_prot = re.findall('var resultsize = .*;', line)
            number_of_prot_2 = number_of_prot[0]
            number_of_prot_2 = int(number_of_prot_2[17:-1])
            print('> Number of Uniprot Entries Identified: ', number_of_prot_2)

        except Exception as e:
            print("ERROR: %s"%e)
            print('Failed obtaining Uniprot codes for ' + chromosome)
            if end == 1:
                break
            
            else:
                continue

        #produce a list going up in increments of 25 (0, 25, 50... [number of entries]).
        #This list is needed so for the URL later on so that the code cycles through the uniprot website
        #25 proteins at a time until the end.
        while starting_number < number_of_prot_2:
            list_of_entries.append(starting_number)
            starting_number = starting_number + 25

        list_of_entries.append(number_of_prot_2)

        #go through uniprot 25 entries at a time and appends uniprot codes of proteins to a list.       
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
                print('Failed to obtain Uniprot codes. %s'%e)
                continue
                
            if len(list_UNIPROT_codes) == 0:
                print('No Uniprot codes were found on: ' + web_url)
                continue

            #For each uniprot code identified it then parses through the .txt file
            #and appends the Uniprot code, PDB code, method structure obtained by, resolution and chain information
            #to a df for each PDB code.
            for protein_code_clean in list_UNIPROT_codes:

                try:
                    df = get_protein_structures(protein_code_clean, df=df)
                except Exception as e:
                    print("%s", e)
                    continue
            
            #Here the code calls the get_chains function which gets chain information for the protein from the PDB.
            number = number + 25

            if number > number_of_prot_2:
                print('Number of Uniprot Entries Searched: ' + str(number_of_prot_2) + '\n')
            else:
                print('Number of Uniprot Entries Searched: ' + str(number) + '\n')
                
            if number > number_of_prot_2:
                number = number_of_prot_2
                
            print((str((number/number_of_prot_2)*100))[0:3] + '%')
            
    return df


def get_protein_structures(uniprot_code, pdb_code_target="", df=[], done_pdbs=[]):
    '''
    Download protein structures associated with a given UNIPROT code.
    The protein is then cleaned (keep only protein atoms, remove hydrogens, MSE and KCX amino acids, split alternative conformations in multiple PDBs)
    If a PDB code is also provided, only that PDB will be downloaded (e.g. useful for consistency check between UNIPROT and PDB)
    If a DataFrame df is provided, extracted structures will be appended to it
    If a list pdb previously analysed is provided in done_pdbs, only PDB not present in the list will be processed.
    '''    
    
    #if no DataFrame is provided, build an empty one
    if len(df) == 0:
        df = create_empty_dataframe()

    #check if there is uniprot information available for the protein
    try:
        url_2 = 'https://www.uniprot.org/uniprot/' + uniprot_code + '.txt'
        html_2 = urllib.request.urlopen(url_2)

    except Exception as e:
        raise Exception('Failed to obtain UNIPROT data. %s'%e)
    
    if pdb_code_target != "":
    
        #appends the AF structure to the df (if this isn't present it will be removed later).
        try:
    
            AF_code = 'AF-' + uniprot_code + '-F1-model_v1'
            
            if AF_code in done_pdbs:
                pass
            
            else:
                data = ({'Uniprot Entry': uniprot_code, 'PDB Code': AF_code, 'Method': 'Predicted', 'Resolution': 'N/A', 'Chains': 'N/A'})
                df = pd.concat([df, pd.DataFrame.from_records(data, index=[0])], ignore_index=True)
    
        #search for available PDB structures
        except Exception as e:
            print('Error %s'%e)

    for line in html_2:
        try:
            line = str(line)
            messy_entry = re.findall('PDB; [\w -. ; \d /]*=', line)

            for entry in messy_entry:
                words = entry.split()
                
                PDBCODE = (words[1])[:-1]
                
                if pdb_code_target != "" and PDBCODE != pdb_code_target:
                    continue
     
                #check if the pdb is in the log file    
                keep = checks(PDBCODE, done_pdbs)
                if keep == False:
                    continue
     
                method_obtained = (words[2])[:-1]
                resolution = (words[3])[:-1]
     
                # if the PDB is not in the log file
                # load, clean, and split it in alternate conformations
                pl.clean_and_split_alt_conformations(PDBCODE)
              
                files = np.array(glob.glob(os.path.join("curate_PDB", "conformations", "*%s*pdb"%PDBCODE)))
                for f in files:
                    conf = f[25:-4]
          
                    #Get chains obtains relevant chain information.
                    unique_values_chain = get_chains(f)

                    #append the data to the df
                    for i in range(len(unique_values_chain)):
                        data = ({'Uniprot Entry': uniprot_code, 'PDB Code': conf, 'Method': method_obtained, 'Resolution': resolution, 'Chains': unique_values_chain[i]})
                        df = pd.concat([df, pd.DataFrame.from_records(data, index=[0])], ignore_index=True)

        except Exception as e:
            print('Error %s'%e)
            continue

    return df


def get_chains(pdb):
    '''
    Get list of chains in a given PDB code
    '''

    try:
        M = bb.Molecule()
        M.import_pdb(pdb, include_hetatm=True)
        unique_values_chain = list(set(list(M.data['chain'])))

    except Exception as e:
        print("Failed to obtain chain information.\n A single .pdb structure may not exist because the protein is too large.")
        raise Exception("%s"%e)

    return unique_values_chain


def search_by_technique(list_of_techniques, wanted_res, pdb_codes_df=[]):
    '''
    Given a list of techniques from the user and the desired resolution
    this function removes pdb entries from pdb_codes_df that don't fit the criteeria
    the second returned parameters defined whether, upon failure, main program should continue (if False, it stops)
    '''    

    if len(pdb_codes_df) == 0:
        pdb_codes_df = create_empty_dataframe()

    try:
        pdb_codes_df = pdb_codes_df[pdb_codes_df['Method'].isin(list_of_techniques)]

    except Exception as e:
        print('Failed to search by method %s'%e)
        question = [
        inquirer.List('Choice',
                            message="Continue?",
                            choices=['Yes', 'No'],
                        ),]
        answer = inquirer.prompt(question)

        if answer['Choice'] == 'Yes':
            return True
        elif answer['Choice'] == 'No':
            return False

    try:
        resolution = float(wanted_res)

        pdb_codes_df['Resolution'] = pdb_codes_df['Resolution'].astype(float)
        pdb_codes_df = pdb_codes_df[pdb_codes_df['Resolution']<= resolution]

    except Exception as e:
        print('Failed to search by resolution. %s'%e)
        question = [
        inquirer.List('Choice',
                            message="Continue?",
                            choices=['Yes', 'No'],
                        ),
        ]
        answer = inquirer.prompt(question)

        if answer['Choice'] == 'Yes':
            return True
        elif answer['Choice'] == 'No':
            return False
            

    return pdb_codes_df, True


def from_csv_file(csv_file, pdb_codes_df=[], done_pdbs=[]):
    '''
    Parse a .csv file to find uniprot and pdb codes to pass into the pipeline.
    If you want to just input uniprot codes put them in the first column and leave the second empty
    If you want to input pdb codes put the uniprot code in the first column and pdb code in the second.
    '''
    
    if len(pdb_codes_df) == 0:
        pdb_codes_df = create_empty_dataframe()

    #read .csv file.
    try:
        csv_df = pd.read_csv(csv_file)
        print('.csv file successfully opened')
        
    except Exception as e:
        raise Exception('Failed to read %s. %s'%(csv_file, e))

    #Next abstract column names
    try:
        column_names = list(csv_df.columns)
        csv_df[column_names[1]] = csv_df[column_names[1]].fillna(0)

    except Exception as e:
        raise Exception('Failed to get data from .csv file. %s'%e)
        
    #get the data about the protein and append it to the pdb_codes_df dataframe.
    for i in range(len(csv_df)):
        try:
            uniprot_code = csv_df.at[i, column_names[0]]
            pdb_code = csv_df.at[i, column_names[1]]
            
            print("\n", uniprot_code, pdb_code)
            
            if pdb_code == 0:
                pdb_codes_df = get_protein_structures(uniprot_code, pdb_codes_df, done_pdbs)

            if pdb_code == 'AF':
                pdb_code = 'AF-' + uniprot_code + '-F1-model_v1'

                keep = True
                if pdb_code in done_pdbs:
                    keep = False

                if keep == True:
                    d = {'Uniprot Entry': uniprot_code, 'PDB Code': pdb_code, 'Method': 'Predicted', 'Resolution': 'N/A', 'Chains': 'A'}
                    pdb_codes_df = pd.concat([pdb_codes_df, pd.DataFrame.from_records(d, index=[0])], ignore_index=True)

            else:
                pdb_codes_df = get_protein_structures(uniprot_code, pdb_code, pdb_codes_df, done_pdbs=done_pdbs)
                
        except Exception as e:
            print('> Error %s'%e)
            continue

    return pdb_codes_df


########################################################

if __name__ == "__main__":
    
    print(get_organism_proteins('Oryctolagus+cuniculus+(Rabbit)', 'UP000001811'))

    #print(get_chains('6YAM'))

    #res_raw = search_by_technique(['X-ray'], pdb_codes_df, '3.0')
    #print(pdb_codes_df)
   
    #pdb_codes_df = get_protein_structures("P09167", "1PRE")
    #print(pdb_codes_df)
    
    #test = from_csv_file("inputs\\input_codes_4.csv")
    #print(test)