import os
import re
import threading
import concurrent.futures
import json
import urllib.request
import requests
from requests.adapters import HTTPAdapter, Retry
import pandas as pd
import numpy as np


class Uniprot(object):
    '''
    Class to handle the first step in the overall carbamylation prediction pipeline. This will take
    in a spreadsheet with the column headers: 'Uniprot_Entry' and 'PDB_Code' (the latter is used to
    specify specific PDB structures, can be left blank to extract all available structures for the
    given UNIPROT entry). The class will extract information on the protein and return a list of
    structures associated to take onto the next step: protein.py.
    '''

    def __init__(self, outdir='data'):
        '''
        Initialise the Uniprot class. This class provides the methods to take a list of UNIPROT
        codes and return a list of structures that can be extracted for use in the overall model.
        
        :param outdir: The name of the output directory. Note this is currently only used for the
            GO term work here.
        '''
        columns = ['Uniprot_Entry', 'PDB_Code', 'Method', 'Resolution', 'Chains']
        self.df = pd.DataFrame(columns=columns)
        self.outdir = outdir
        pd.reset_option('display.max_rows')
        pd.set_option("display.max_columns", None)

        self.GO_dict = {}  # code as the key
        self.uniprot_to_GO = {}  # uniprot as key
        self.name_to_code = {}
        self.code_to_name = None


    def count_organism_proteins(self, code, reviewed_only=False):
        '''
        Contact the uniprot api to enquire about the number of proteins within a given organism
        code. Possible to filter this to be just reviewed proteins.

        :param code: The uniprot code for the organism of interest
        :type code: str
        :param reviewed_only: Tunable option to limit the proteins returned to just those that are
            reviewed
        :type reviewed_only: bool
        :returns: The total number of proteins associated with the organism code
        :rtype: int
        '''
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
        '''
        Contacts the uniprot api to obtain the lists of proteins associated with the given uniprot
        organism code. Possible to filter this to be just the reviewed proteins here.

        :param code: The uniprot code for the organism of interest
        :type code: str
        :param reviewed_only: Tunable option to limit the proteins returned to just those that are
            reviewed
        :type reviewed_only: bool
        :returns: List of proteins associated with the organism code
        :rtype: list
        '''
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
            for line in batch.text.splitlines():
                line = line.strip()
                if line and line != 'Entry':
                    codes.append(line)
            print(f'{len(codes)} / {total}')

        return codes



    def get_protein_data(self, uniprot_code, pdb_code_target="", chain_target=""):
        '''
        Download protein structures associated with a given UNIPROT code. If a PDB code is also
        provided, only that PDB will be downloaded (e.g. useful for consistency check between
        UNIPROT and PDB) If a DataFrame df is provided, extracted structures will be appended to it.
        Extracts the name for AF structure, contacts AF API and finds latest version, if any
        problems in finding the latest version, uses v6.

        :param uniprot_code: The uniprot code for the organism of interest
        :type uniprot_code: str
        :param pdb_code_target: If a specific PDB is required, can be specified here, else all PDBs
            will be obtained
        :type pdb_code_target: str
        :param chain_target: If a specific chain is required, can be specified here, else all chains
            will be obtained
        :type chain_target: str
        '''

        #check if there is uniprot information available for the protein
        try:
            url_2 = 'https://www.uniprot.org/uniprot/' + uniprot_code + '.txt'
            html_2 = urllib.request.urlopen(url_2)

        except Exception as e:
            raise Exception(f'Failed to obtain UNIPROT data for code: {uniprot_code}. {e}') from e

        if pdb_code_target == "":

            #appends the AF structure to the df (if this isn't present it will be removed later).
            try:
                # get latest AF structure details for the uniprot code
                try:
                    af_details_response = requests.get(f'https://alphafold.ebi.ac.uk/api/prediction/{uniprot_code}?include_complexes=false', timeout=10)
                    if not af_details_response.ok:
                        af_details_response.raise_for_status()
                    af_data = af_details_response.json()
                    latest_version = str(af_data[0]['latestVersion'])
                    af_code = f'AF-{uniprot_code}-F1-model_v{latest_version}'
                except Exception as e:
                    print(f'>> Failed to obtain information about AF structure for uniprot code {uniprot_code} with error {e}, using v6 for AF code.')
                    af_code = f'AF-{uniprot_code}-F1-model_v6'

                data = ({'Uniprot_Entry': uniprot_code, 'PDB_Code': af_code, 'Method': 'Predicted', 'Resolution': np.nan, 'Chains': np.nan})
                self.df = pd.concat([self.df, pd.DataFrame.from_records(data, index=[0])], ignore_index=True)

            except Exception as e:
                print(f'Error {e}')

        for line in html_2:
            try:
                line = str(line)
                messy_entry = re.findall(r'PDB; [\w -. ; \d /]*=', line)

                if len(messy_entry)>0:

                    words = messy_entry[0].split("; ")
                    pdb_code = words[1]

                    if pdb_code_target != "" and pdb_code != pdb_code_target:
                        continue

                    method_obtained = words[2]

                    try:
                        resolution = words[3].split()[0]
                    except Exception:
                        resolution = np.nan
                    chains = words[-1][:-1]

                    data = None
                    if chain_target != "":
                        for c in chains.split('/'):
                            if chain_target == c:
                                data = {'Uniprot_Entry': uniprot_code, 'PDB_Code': pdb_code,
                                        'Method': method_obtained, 'Resolution': resolution, 'Chains' : c}
                                break
                    else:
                        data = {'Uniprot_Entry': uniprot_code, 'PDB_Code': pdb_code,
                                'Method': method_obtained, 'Resolution': resolution, 'Chains' : chains}

                    if data is not None:
                        self.df = pd.concat([self.df, pd.DataFrame.from_records(data, index=[0])], ignore_index=True)

            except Exception as e:
                print(f'Error {e}')
                continue


    def from_organism(self, code, reviewed_only=False):
        '''
        Function to go from an organism to dataframe of proteins ready to be curated in the protein
        class. This will take an organism code, find all the proteins associated with this proteome,
        then iterate over the returned list of structures to append to the class output dataframe
        containing full details

        :param code: Code for the organism
        :type code: str
        :param reviewed_only: Toggle to return only reviewed proteins or all
        :type reviewed_only: bool
        '''
        prot_list = self.get_organism_proteins(code=code, reviewed_only=reviewed_only)
        for i, prot in enumerate(prot_list):
            print(f'Extracting protein information for protein {prot} in organism '
                  f'{code}; {i}/{len(prot_list)}')
            self.get_protein_data(uniprot_code=prot)
        self.df = self.df.drop_duplicates()


    def from_csv_file(self, csv_file):
        '''
        Parse a .csv file to find uniprot and pdb codes to pass into the pipeline. If you want to
        just input uniprot codes put them in the first column and leave the second empty. If you
        want to input pdb codes put the uniprot code in the first column and pdb code in the second.

        :param csv_file: Path to the csv file containing the desired uniprot codes and optional pdb
            codes
        :type csv_file: str
        '''

        try:
            csv_df = pd.read_csv(csv_file)
            print('.csv input file of Uniprot codes successfully opened')
        except Exception as e:
            raise Exception(f'Failed to read {csv_file}: {e}') from e

        try:
            csv_df['PDB_Code'] = csv_df['PDB_Code'].fillna(0)
        except Exception as e:
            raise Exception(f'Failed to get data from .csv file. {e}') from e

        for i, r in csv_df.iterrows():
            try:
                uniprot_code = r['Uniprot_Entry']
                pdb_code = r['PDB_Code']

                if pdb_code == 0:
                    self.get_protein_data(uniprot_code)

                elif pdb_code == 'AF':
                    pdb_code = f'AF-{uniprot_code}-F1-model_v6'

                    d = {'Uniprot_Entry': uniprot_code, 'PDB_Code': pdb_code, 'Method': 'Predicted', 'Resolution': np.nan, 'Chains': np.nan}
                    self.df = pd.concat([self.df, pd.DataFrame.from_records(d, index=[0])], ignore_index=True)

                else:
                    self.get_protein_data(uniprot_code, pdb_code)

            except Exception as e:
                print(f'> Error {e}')
                continue

        self.df = self.df.drop_duplicates()


    def filter_by_technique(self, list_of_techniques):
        '''
        Return a subset of the proteins obtained from the search based on the method in which the
        structure was obtained; options: x-ray, NMR, EM, Predicted

        :param list_of_techniques: The list of techniques to filter the proteins obtained by
        :type list_of_techniques: list
        :returns: Dataframe containing the subset of proteins
        '''
        return self.df[self.df['Method'].isin(list_of_techniques)]


    def save_state(self, outname='potential_proteins.csv'):
        '''
        Quick function to aid saving a list of protein structures that have been found ready to be
        used in the protein class.

        :param outname: Name of the file to write the UP.df dataframe to
        :type outname: str
        '''
        self.df.to_csv(outname, index_label=False, index=False)


    def GO_search_term(self, df, code = '', name = ''):
        '''
        List the subset of UNIPROT codes associated with a GO Term.

        :param df: Dataframe to get matches to GO term from
        :type  df: pandas.DataFrame
        :param code: Code of the Go term to analyse; format: 7 numbers only
        :type code: str
        :param name: Name corresponding to the GO term
        :type name: str
        '''
        if not self.GO_dict:
            print('>> No GO terms have been searched for and as such this would return '
                  'nothing. Please call GO_get_data() on a list of uniprot entries and '
                  'then recall this function as needed.')
            return pd.DataFrame()

        if code == '' and name == '':
            return 'Insufficient input!'
        elif code == '' and name != '':
            code = self.name_to_code[name]
        elif code != '' and name != '':
            # check if they match
            if code != self.name_to_code[name]:
                return f'Unmatched GO term code and name, wrong input code {code}.'

        uni_list = self.GO_dict.get(code)
        df_out = pd.DataFrame()
        for uni in uni_list:
            cdf = df[df['Uniprot_Entry'] == uni]
            df_out = pd.concat([df_out, cdf], ignore_index=True)

        return df_out


    def GO_search_protein(self, uniprot_entry):
        '''
        List all the GO Terms associated with a specific UNIPROT code based
        on previous work done by calling GO_get_data() to extract all the information
        on the dataframe for GO terms. If GO_get_data() has not been already called,
        this function will call it before searching the GO dictionary to give a list
        of GO terms associated.
        
        :param uniprot_code: The uniprot code of interest to do an API call on to get
            the GO terms associated with this code.
        :type uniprot_code: str
        '''

        if self.GO_dict == {}:
            self.GO_get_data(uniprot=uniprot_entry)
        GO_list = []
        for code, uni_list in self.GO_dict.items():
            if uniprot_entry in uni_list:
                GO_list.append(code)
        GO_list = [a for a in [self.code_to_name.get(code) for code in GO_list] if a is not None]
        return GO_list


    def _GO_get_data(self, uniprot_code, lock, index, total):
        '''
        Worker function of the self.GO_get_data method. This performs and API call to
        the UNIPROT website to obtain the details of the GO codes associated with the
        specific UNIPROT code of interest. After details have been found, it assigns the
        GO terms to the GO_dict key with the GO term as the key and the UNIPROT code
        in a list as the value. If already in the dict, it just adds the UNIPROT code to
        the list associated. It also adds it into the mapping between the GO term and
        the GO name associated.
        
        :param uniprot_code: The uniprot code of interest to do an API call on to get
            the GO terms associated with this code.
        :type uniprot_code: str
        :param lock: Lock provided by the multithreading instance for adding data out to
            the dicionaries
        :type lock: multithreading process lock
        :param index: Number of the protein to investigate in the dataframe
        :type index: int
        :param total: Total number of the proteins to investigate in the dataframe
        :type total: int
        '''

        print(f'Searching for {index}/{total} protein.')
        try:
            url_2 = f'https://www.uniprot.org/uniprot/{uniprot_code}.txt'
            with urllib.request.urlopen(url_2, timeout=10) as response:
                html_2 = response.read().decode('utf-8')
        except Exception as e:
            print(f'Failed to obtain UNIPROT data for {uniprot_code}. {e}')

        try:
            for line in html_2.splitlines():
                line = str(line)
                find_go = re.findall('GO; GO', line)
                if len(find_go) > 0:
                    # GO; GO:0030089; C:phycobilisome; IEA:UniProtKB-KW.
                    # GO; GO:0102834; F:1-18:1-2-16:0-monogalactosyldiacylglycerol acyl-lipid omega-6 desaturase activity; IEA:UniProtKB-EC.
                    GO_code = line.split('; ')[1].split(':')[1]
                    idx = line.split('; ')[2].index(':') + 1
                    GO_word = line.split('; ')[2][idx:]
                    with lock:
                        if GO_code not in self.GO_dict.keys():
                            self.GO_dict[GO_code] = [uniprot_code]
                        else:
                            self.GO_dict[GO_code].append(uniprot_code)

                        if uniprot_code not in self.uniprot_to_GO.keys():
                            self.uniprot_to_GO[uniprot_code] = [GO_code]
                        else:
                            self.uniprot_to_GO[uniprot_code].append(GO_code)

                        if GO_word not in self.name_to_code.keys():
                            self.name_to_code[GO_word] = GO_code

        except Exception as e:
            print(f'Error analysis GO data for Uniprot: {uniprot_code} with error: {e}')


    def GO_get_data(self, uniprot):
        '''
        Go over the measures dataframe given as input and extract all possible GO terms
        that have associated with the Uniprot codes within the dataframe. This calls a
        worked function _GO_get_data iteratively to do this automatically and will create
        class dictionaries for converting a GO code to a name and name to a GO code.
        
        :param uniprot: Input either a uniprot code or list of uniprot codes and then the
            function will go through the codes and obtain the GO terms associated.
        :type uniprot: str, list
        '''

        if isinstance(uniprot, str):
            uniprot = [uniprot]
        uniprot_codes = list(set(uniprot))
        num = len(uniprot_codes)
        locks = [threading.Lock()]*num
        with concurrent.futures.ThreadPoolExecutor() as executor:
            executor.map(self._GO_get_data, uniprot_codes, locks, range(num), [num]*num)

        self.code_to_name = {v: k for k, v in self.name_to_code.items()}


    def save_GO_metadata(self, outname='GO_metadata.json'):
        '''
        Save the GO term metadata that has been extracted to a json file. This will save the
        UNIPROT code to GO term data, which can be loaded in to use later on or loaded in here
        to avoid searching across a full dataset again. It will also save the GO terms to GO
        names dict to allow full recovery after
        '''

        if not self.GO_dict:
            print('>> No UNIPROT codes have had GO terms extracted for and as such there is nothing'
                  'to write to a file. Please collect some GO data before saving.')
            return

        metadata = {k: {'uniprot': self.GO_dict[k], 'name': self.code_to_name[k]}
                    for k in self.GO_dict}
        json_GO_metadata = json.loads(json.dumps(metadata, default=str))

        with open(os.path.join(self.outdir, outname), 'w') as fh:
                json.dump(json_GO_metadata, fh, indent=2)


    def load_GO_metadata(self, metadata_name='GO_metadata.json'):
        '''
        Load the GO metadata saved from previous investigation into this. Will setup the two GO term
        dictionaries for use, however will not setup the 
        '''

        GO_metadata_path = os.path.join(self.outdir, metadata_name)
        if not os.path.exists(GO_metadata_path):
            print(f'>> Path given to the GO metadata json file does not exist at location:'
                  f'{GO_metadata_path}, could not read anything.')
            return

        with open(GO_metadata_path) as fh:
            json_go_data = json.load(fh)

        for GO_code in json_go_data.keys():
            go_metadata = json_go_data[GO_code]
            name = go_metadata['name']
            uniprot = go_metadata['uniprot']
            if isinstance(uniprot, str):
                uniprot = [uniprot]

            if GO_code not in self.GO_dict.keys():
                self.GO_dict[GO_code] = uniprot
            else:
                for uni in uniprot:
                    self.GO_dict[GO_code].append(uni)
                self.GO_dict[GO_code] = list(set(self.GO_dict[GO_code]))

            if GO_code not in self.code_to_name.keys():
                self.code_to_name[GO_code] = name

        self.name_to_code = {v:k for k,v in self.code_to_name.items()}



if __name__ == "__main__":

    UP = Uniprot()

    #print(UP.get_organism_proteins('UP000001806'))
    UP.from_organism('UP000000625')
    #UP.get_protein_data("P09167")
    #UP.from_csv_file("inputs\\input_codes_4.csv")
    print(UP.df)
