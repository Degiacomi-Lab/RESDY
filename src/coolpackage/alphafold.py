import csv
import os
import requests
from .helper import get_download_tool

def download_AF_struc(pdb, outfolder="result"):
    '''
    Download AlphaFold2 structures into the assembled folder.

    .. rubric:: Method

    - Check that curated directory exists in outdir, otherwise create one
    - Use the appropriate download tool to download the AF structure for the code given

    :param pdb: The AF code for the structure to extract the PLDDT values from
    :type pdb: str
    :param outfolder: The output directory used to know where the downloaded files should be
        written to
    :type outfolder: str

    .. rubric:: Example

    ::

        >>> download_AF_struc('AF-P0CG48-F1-model_v6', outfolder='test_plddt')
    '''

    download_path = os.path.join(outfolder, "curated")
    if not os.path.exists(download_path):
        os.makedirs(download_path)

    if os.path.exists(os.path.join(download_path, f'{pdb}.pdb')):
        print(f'>> AF structure for {pdb} is already downloaded, using previous version')
        return

    print(f"> downloading AlphaFold structure {pdb}")

    try:
        response = requests.get(url=f'https://alphafold.ebi.ac.uk/files/{pdb}.pdb', timeout=20)
        response.raise_for_status()
        open(os.path.join(download_path, f'{pdb}.pdb'), 'wb').write(response.content)

    except Exception as e:
        print(f'>> AF structure not found for {pdb}, error: {e}')
        raise Exception(f'AF structure not found for {pdb}: {e}') from e

    return

def find_af_plddt(af_code_full, outfolder="result", resnames=['LYS']):
    '''
    Obtain PLDDT (a measure of certainty where 100 is high and 70 low) value for each lysine in an
    alphafold structure.

    .. rubric:: Method

    - If the PLDDT output file hasn't been created yet, create and add column headings
    - Open the PLDDT output file for appending data.
    - Search the curated AF structure pdb file to extract the PLDDT values.
    - Append the PLDDT value to the output file and write to the console log.

    :param af_code_full: The AF code for the structure to extract the PLDDT values from
    :type af_code_full: str
    :param outfolder: The output directory used to check that the PLDDT out file is in the correct
        place and allow appending to this
    :type outfolder: str
    :param resnames: list of resid names to search for in the pdb files. Default is ['LYS'], set in
        overall run through curation on aa_properties
    :type resnames: list
    :returns: A dictionary matching up all the lysines with their corresponding PLDDT values for the
        given AF structure.
    :rtype: dict

    .. rubric:: Example

    ::

        >>>find_af_plddt('AF-P0CG48-F1-model_v4', outfolder='test_plddt')
        > Finding plddt
        AF-P0CG48-F1-model_v4; Resid No. A6; PLDDT: 93.79
        AF-P0CG48-F1-model_v4; Resid No. A11; PLDDT: 89.45
        AF-P0CG48-F1-model_v4; Resid No. A27; PLDDT: 94.28
        ...
    '''
    cols = ['PDB_Code', 'Chain', 'Resid', 'PLDDT']
    if not os.path.isfile(os.path.join(outfolder, "curated", "AF_PLDDT_Output.csv")):
        with open(os.path.join(outfolder, "curated", "AF_PLDDT_Output.csv"), 'w', newline='') as plddt_out_file:
            plddt_writer = csv.writer(plddt_out_file)
            plddt_writer.writerow(['PDB_Code', 'Chain', 'Resid', 'PLDDT'])

    plddt_out_file = open(os.path.join(outfolder, "curated", "AF_PLDDT_Output.csv"), 'a', newline='')
    plddt_writer = csv.DictWriter(plddt_out_file, fieldnames=cols)

    print(f'>> Finding plddt for AF structure {af_code_full}')

    dict_plddt = dict()

    try:
        with open(os.path.join(outfolder, "curated", af_code_full + ".pdb"), "r") as f:
            for line in f:
                try:
                    line = str(line)
                    if line[12:16].strip() == 'CA' and line[17:20].strip() in resnames:
                        chain = line[21]
                        resid = line[22:26].strip()
                        plddt = line[60:66].strip()
                        chain_resid = chain + resid

                        dict_plddt.update({chain_resid: plddt})

                        print(af_code_full + "; Chain: " + str(chain_resid[0]) + "; Resid: " + str(chain_resid[1:]) + "; PLDDT: " + plddt)
                        plddt_writer.writerow({'PDB_Code': af_code_full, 'Chain': chain_resid[0], 'Resid': chain_resid[1:], 'PLDDT': plddt})
                        chain, resid, plddt = '', '', ''

                except Exception as e:
                    print(f"Error {e}")
                    plddt_writer.writerow({'PDB_Code': af_code_full,
                                           'Chain': line[21:22],
                                           'Resid': line[22:26].strip(),
                                           'PLDDT': f'Error {str(e)}'})
                    continue

    except Exception as e:
        print(f'Failed to obtain PLDDT data for {af_code_full}; error: {e}')
        plddt_writer.writerow({'PDB_Code': af_code_full, 'Chain': '', 'Resid': '',
                               'PLDDT': f'Error {str(e)}'})

    finally:
        plddt_out_file.close()
    return dict_plddt


if __name__ == '__main__':
    download_AF_struc('AF-P0CG48-F1-model_v6', outfolder='result')
    #find_af_plddt('AF-P0CG48-F1-model_v6', outfolder='result')
