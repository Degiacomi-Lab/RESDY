import re
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

    def __init__(self):
        '''
        Initialise the Uniprot class. This class provides the methods to take a list of uniprot
        codes and return a list of structures that can be extracted for use in the overall model.
        '''
        columns = ['Uniprot_Entry', 'PDB_Code', 'Method', 'Resolution', 'Chains']
        self.df = pd.DataFrame(columns=columns)
        pd.reset_option('display.max_rows')
        pd.set_option("display.max_columns", None)


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

            #search for available PDB structures
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

        #read .csv file.
        try:
            csv_df = pd.read_csv(csv_file)
            print('.csv input file of Uniprot codes successfully opened')
        except Exception as e:
            raise Exception(f'Failed to read {csv_file}: {e}') from e

        #Next abstract column names
        try:
            csv_df['PDB_Code'] = csv_df['PDB_Code'].fillna(0)
        except Exception as e:
            raise Exception(f'Failed to get data from .csv file. {e}') from e

        #get the data about the protein and append it to the dataframe.
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


if __name__ == "__main__":

    UP = Uniprot()

    #print(UP.get_organism_proteins('UP000001806'))
    UP.from_organism('UP000000625')
    #UP.get_protein_data("P09167")
    #UP.from_csv_file("inputs\\input_codes_4.csv")
    print(UP.df)
