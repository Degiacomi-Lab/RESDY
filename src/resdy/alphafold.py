import csv
import os
import json
import requests
import pandas as pd


def download_AF_struc(pdb, outfolder="result"):
    '''
    Download AlphaFold2 structures into the assembled folder.

    .. rubric:: Method

    - Check that curated directory exists in outdir, otherwise create one
    - Use the appropriate download tool to download the AF structure for the code given

    :param pdb: The AF code for the structure to extract the PLDDT values from.
    :type pdb: str
    :param outfolder: The output directory used to know where the downloaded files should be
        written to. This will then write to a folder named 'curated' inside this.
    :type outfolder: str
    '''

    download_path = os.path.join(outfolder, "curated")
    if not os.path.exists(download_path):
        os.makedirs(download_path)

    if os.path.exists(os.path.join(download_path, f'{pdb}.pdb')):
        print(f'>> AF structure for {pdb} is already downloaded, using previous version')
        return

    print(f"> downloading AlphaFold structure {pdb}")

    try:
        response = requests.get(url=f'https://alphafold.ebi.ac.uk/files/{pdb}.pdb',
                                timeout=20)
        response.raise_for_status()
        open(os.path.join(download_path, f'{pdb}.pdb'), 'wb').write(response.content)

    except Exception as e:
        print(f'>> AF structure not found for {pdb}, error: {e}')
        raise Exception(f'AF structure not found for {pdb}: {e}') from e

    return


def find_af_plddt(af_code_full, outfolder='result', resnames=['LYS']):
    '''
    Obtain PLDDT (a measure of certainty where 100 is high and 70 low) value for each residue
    of interest in an alphafold structure. This will be reported into a json file per structure.
    After the curation has been completed, all af_plddt json files will be extracted and collated
    into the overall AF_PLDDT_Output.csv file.

    .. rubric:: Method

    - Search the curated AF structure pdb file to extract the PLDDT values.
    - Write a .json file for the af_code that is PLDDT data is extracted from

    :param af_code_full: The AF code for the structure to extract the PLDDT values from
    :type af_code_full: str
    :param outfolder: The output directory used to check that the PLDDT out file is in the correct
        place and allow appending to this
    :type outfolder: str
    :param resnames: list of resid names to search for in the pdb files. Default is ['LYS'], set in
        overall run through curation on aa_properties
    :type resnames: list
    :returns: A dictionary matching up all the residues of interest up with their corresponding
        PLDDT values for the given AF structure.
    :rtype: dict
    '''
    print(f'>> Finding plddt for AF structure {af_code_full}')

    dict_plddt = dict()
    os.makedirs(os.path.join(outfolder, 'curated'), exist_ok=True)
    af_plddt_json_path = os.path.join(outfolder, 'curated', f'{af_code_full}_PLDDT.json')

    try:
        with open(os.path.join(outfolder, "curated", af_code_full + ".pdb"), "r") as f:
            for line in f:
                try:
                    line = str(line)
                    if line[12:16].strip() == 'CA' and line[17:20].strip() in resnames:
                        chain = line[21]
                        resid = line[22:26].strip()
                        plddt = line[60:66].strip()
                        chain_res_key = f'{chain}_{str(resid)}'

                        dict_plddt.update({chain_res_key: plddt})

                        print(af_code_full + "; Chain: " + str(chain) +
                                "; Resid: " + str(resid) +
                                "; PLDDT: " + plddt)
                        chain, resid, plddt = '', '', ''

                except Exception as e:
                    chain = line[21]
                    resid = line[22:26].strip()
                    chain_res_key = f'{chain}_{str(resid)}'
                    plddt = f'Error {str(e)}'
                    dict_plddt.update({chain_res_key: plddt})
                    print(f'Error with AF PLDDT value for PDB code: {af_code_full}, chain: '
                          f'{chain}, resid: {resid}; Error: {str(e)}')

                    continue

    except Exception as e:
        print(f'Failed to obtain PLDDT data for {af_code_full}; error: {e}')
        dict_plddt = {('', ''): f'Error {str(e)}'}

    plddt_out = {af_code_full: dict_plddt}

    with open(af_plddt_json_path, 'w') as fout:
        json.dump(plddt_out, fout, indent=2, )

    print(f'>> Written file {af_plddt_json_path} for PLDDT values for {af_code_full}')

    return plddt_out


def report_af_plddt(outfolder='result'):
    '''
    Collate all PLDDT report json files from structures together into one overall csv file.

    .. rubric:: Method

    - If the PLDDT output file hasn't been created yet, create and add column headings
    - Open the PLDDT output file for appending data.
    - For each json PLDDT report available in the curated folder within the given outfolder

    :param outfolder: The output directory used to check that the PLDDT out file is in the correct
        place and allow appending to this
    :type outfolder: str
    '''
    cols = ['PDB_Code', 'Chain', 'Resid', 'PLDDT']
    os.makedirs(os.path.join(outfolder, 'curated'), exist_ok=True)
    plddt_file_path = os.path.join(outfolder, 'curated', 'AF_PLDDT_Output.csv')

    if not os.path.isfile(plddt_file_path):
        with open(plddt_file_path, 'w', newline='') as plddt_out_file:
            plddt_writer = csv.writer(plddt_out_file)
            plddt_writer.writerow(['PDB_Code', 'Chain', 'Resid', 'PLDDT'])

    plddt_json_files = [f for f in os.listdir(os.path.join(outfolder, 'curated'))
                        if '_PLDDT.json' in f]

    plddt_out_file = open(plddt_file_path, 'a', newline='')
    plddt_writer = csv.DictWriter(plddt_out_file, fieldnames=cols)
    for file in plddt_json_files:
        f_p = os.path.join(outfolder, 'curated', file)
        if os.path.isfile(f_p):
            with open(f_p) as fh:
                plddt_data = json.load(fh)

            af_code = list(plddt_data)[0]
            plddts = plddt_data[af_code]

            for k, plddt in plddts.items():
                if '_' in k:
                    chain = k.split('_')[0]
                    resid = k.split('_')[1]
                else:
                    chain = ''
                    resid = ''
                plddt_writer.writerow({'PDB_Code': af_code,
                                    'Chain': chain,
                                    'Resid': resid,
                                    'PLDDT': plddt})

    df_plddt = pd.read_csv(plddt_file_path)
    if len(df_plddt.drop_duplicates()) != len(df_plddt):
        print('>> Duplicates found in the PLDDT report entries, keeping the newest entries.')
        df_plddt = df_plddt.drop_duplicates(keep='last')

    plddt_out_file.close()



if __name__ == '__main__':
    download_AF_struc('AF-P0CG48-F1-model_v6', outfolder='result')
    #find_af_plddt('AF-P0CG48-F1-model_v6', outfolder='result')
