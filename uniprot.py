import re
import urllib.request, urllib.parse, urllib.error
import requests
from requests.adapters import HTTPAdapter, Retry

import pandas as pd
import numpy as np


class Uniprot(object):
    
    def __init__(self, done_pdbs=[]):

        columns = ['Uniprot_Entry', 'PDB Code', 'Method', 'Resolution', 'Chains']
        self.df = pd.DataFrame(columns=columns)
  
        
    def count_organism_proteins(self, code, reviewed_only=False):
    
       if reviewed_only:
           url = f'https://rest.uniprot.org/uniprotkb/search?format=list&query=%28%28proteome%3A{code}%29%29%20AND%20%28reviewed%3Atrue%29&size=500'
       else:
           url = f'https://rest.uniprot.org/uniprotkb/search?format=list&query=%28%28proteome%3A{code}%29%29&size=500'
       
       retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
       session = requests.Session()
       session.mount("https://", HTTPAdapter(max_retries=retries))
       response = session.get(url)
       response.raise_for_status()
       total = response.headers["x-total-results"]
       
       return total  
    
    
    def get_organism_proteins(self, code, reviewed_only=False):
    
        if reviewed_only:
            url = f'https://rest.uniprot.org/uniprotkb/search?format=list&query=%28%28proteome%3A{code}%29%29%20AND%20%28reviewed%3Atrue%29&size=500'
        else:
            url = f'https://rest.uniprot.org/uniprotkb/search?format=list&query=%28%28proteome%3A{code}%29%29&size=500'
    
        re_next_link = re.compile(r'<(.+)>; rel="next"')
        retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
        session = requests.Session()
        session.mount("https://", HTTPAdapter(max_retries=retries))
    
        def _get_next_link(headers):
            if "Link" in headers:
                match = re_next_link.match(headers["Link"])
                if match:
                    return match.group(1)
    
        def _get_batch(batch_url):
            while batch_url:
                response = session.get(batch_url)
                response.raise_for_status()
                total = response.headers["x-total-results"]
                yield response, total
                batch_url = _get_next_link(response.headers)
    
        codes = []
        for batch, total in _get_batch(url):
            for line in batch.text.splitlines()[1:]:
                codes.append(line)
            print(f'{len(codes)} / {total}')
        
        return codes
    

    
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
        
                AF_code = 'AF-' + uniprot_code + '-F1-model_v3'
              
                data = ({'Uniprot_Entry': uniprot_code, 'PDB Code': AF_code, 'Method': 'Predicted', 'Resolution': np.nan, 'Chains': np.nan})
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
                                data = {'Uniprot_Entry': uniprot_code, 'PDB Code': PDBCODE, 'Method': method_obtained, 'Resolution': resolution, 'Chains' : c}
                    else:
                        data = {'Uniprot_Entry': uniprot_code, 'PDB Code': PDBCODE, 'Method': method_obtained, 'Resolution': resolution, 'Chains' : chains}
                     
                    
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
    
                    d = {'Uniprot_Entry': uniprot_code, 'PDB Code': pdb_code, 'Method': 'Predicted', 'Resolution': np.nan, 'Chains': np.nan}
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
        UP.get_organism_proteins('UP000001811')
    
    if True:
        UP.get_protein_data("P09167")

    if True:
        UP.from_csv_file("inputs\\input_codes_4.csv")
    
    print(UP.df)
