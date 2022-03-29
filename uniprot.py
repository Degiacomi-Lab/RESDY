import re
import urllib.request, urllib.parse, urllib.error
from bs4 import BeautifulSoup

import pandas as pd
import numpy as np

class Uniprot(object):
    
    def __init__(self, done_pdbs=[]):

        columns = ['Uniprot Entry', 'PDB Code', 'Method', 'Resolution', 'Chains']
        self.df = pd.DataFrame(columns=columns)
                        
   
    def get_organism_proteins(self, name_of_organism, code):
        '''
        Obtain all the PDB codes belonging to an organism.
        The name has to be exactly that used on the uniprot website and the code needs to be the code in the URL for the proteome.    
        '''
           
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
                        self.get_protein_data(protein_code_clean)
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
    
    
    def get_protein_data(self, uniprot_code, pdb_code_target="", chain_target=""):
        '''
        Download protein structures associated with a given UNIPROT code.
        The protein is then cleaned (keep only protein atoms, remove hydrogens, MSE and KCX amino acids, split alternative conformations in multiple PDBs)
        If a PDB code is also provided, only that PDB will be downloaded (e.g. useful for consistency check between UNIPROT and PDB)
        If a DataFrame df is provided, extracted structures will be appended to it
        '''    
    
        
        #check if there is uniprot information available for the protein
        try:
            url_2 = 'https://www.uniprot.org/uniprot/' + uniprot_code + '.txt'
            html_2 = urllib.request.urlopen(url_2)
    
        except Exception as e:
            raise Exception('Failed to obtain UNIPROT data. %s'%e)
       
        if pdb_code_target == "":
        
            #appends the AF structure to the df (if this isn't present it will be removed later).
            try:
        
                AF_code = 'AF-' + uniprot_code + '-F1-model_v1'
              
                data = ({'Uniprot Entry': uniprot_code, 'PDB Code': AF_code, 'Method': 'Predicted', 'Resolution': np.nan, 'Chains': np.nan})
                self.df = pd.concat([self.df, pd.DataFrame.from_records(data, index=[0])], ignore_index=True)
        
            #search for available PDB structures
            except Exception as e:
                print('Error %s'%e)
    
        for line in html_2:
            try:
                line = str(line)
                messy_entry = re.findall('PDB; [\w -. ; \d /]*=', line)
    
                if len(messy_entry)>0:
    
                    words = messy_entry[0].split("; ")
                    PDBCODE = (words[1])
                    
                    if pdb_code_target != "" and PDBCODE != pdb_code_target:
                        continue

         
                    method_obtained = (words[2])
                    
                    try:
                        resolution = (words[3].split()[0])
                    except Exception:
                        resolution = np.nan
                    chains = words[-1][:-1]
    
                    if chain_target != "":      
                        for j, c in enumerate(chains.split('/')):
                            if chain_target == c:     
                                data = {'Uniprot Entry': uniprot_code, 'PDB Code': PDBCODE, 'Method': method_obtained, 'Resolution': resolution, 'Chains' : c}
                    else:
                        data = {'Uniprot Entry': uniprot_code, 'PDB Code': PDBCODE, 'Method': method_obtained, 'Resolution': resolution, 'Chains' : chains}
                     
                    
                    self.df = pd.concat([self.df, pd.DataFrame.from_records(data, index=[0])], ignore_index=True)
    
            except Exception as e:
                print('Error %s'%e)
                continue
    
    
    def from_csv_file(self, csv_file):
        '''
        Parse a .csv file to find uniprot and pdb codes to pass into the pipeline.
        If you want to just input uniprot codes put them in the first column and leave the second empty
        If you want to input pdb codes put the uniprot code in the first column and pdb code in the second.
        '''
    
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
    
        #get the data about the protein and append it to the dataframe.
        for i in range(len(csv_df)):
            try:
                uniprot_code = csv_df.at[i, column_names[0]]
                pdb_code = csv_df.at[i, column_names[1]]
                 
                if pdb_code == 0:
                    self.get_protein_data(uniprot_code)
    
                if pdb_code == 'AF':
                    pdb_code = 'AF-' + uniprot_code + '-F1-model_v1'
    
                    d = {'Uniprot Entry': uniprot_code, 'PDB Code': pdb_code, 'Method': 'Predicted', 'Resolution': np.nan, 'Chains': np.nan}
                    self.df = pd.concat([self.df, pd.DataFrame.from_records(d, index=[0])], ignore_index=True)
    
                else:
                    self.get_protein_data(uniprot_code, pdb_code)
                    
            except Exception as e:
                print('> Error %s'%e)
                continue


    def filter_by_technique(self, list_of_techniques):
        '''
        return the subset of entries obtained by certain methods (x-ray, NMR, EM, Predicted)
        '''
        return self.df[self.df['Method'].isin(list_of_techniques)]
    


########################################################

if __name__ == "__main__":
       
    UP = Uniprot()
 
    if False:
        UP.get_organism_proteins('Oryctolagus+cuniculus+(Rabbit)', 'UP000001811')
    
    
    if True:
        UP.get_protein_data("P09167")

    if True:
        UP.from_csv_file("inputs\\input_codes_4.csv")
    
    print(UP.df)