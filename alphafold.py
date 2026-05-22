import csv
import os
import re
import subprocess
from csv import writer
from helper import get_download_tool

def download_AF_struc(pdb, outfolder="result"):
    '''
    Download AlphaFold2 structures into the assembled folder.

    Method
    ------
    Check that curated directory exists in outdir, otherwise create one
    Use the appropriate download tool to download the AF structure for the code given


    Parameters
    ----------
    pdb : string
        The AF code for the structure to extract the PLDDT values from
    
    outfolder : string
        The outdirectory to used to know where the downloaded files should be written to

    Example
    -------
    >> download_AF_struc('AF-P0CG48-F1-model_v6', outfolder='test_plddt')
    '''

    oldcwd = os.getcwd()
    mypath = os.path.join(outfolder, "curated")
    if not os.path.exists(mypath):
        os.makedirs(mypath)

    os.chdir(mypath)

    print(f"> downloading AlphaFold structure {pdb}")

    tool = get_download_tool()
    try:
        if tool == "curl":
            line = f"curl -s -o {pdb}.pdb https://alphafold.ebi.ac.uk/files/{pdb}.pdb"
        elif tool == "wget":
            line = f"wget https://alphafold.ebi.ac.uk/files/{pdb}.pdb"
        else:
            raise RuntimeError("You don't have a commandline tool for downloading files")

        subprocess.check_call(line, shell=True)

    except Exception as e:
        os.chdir(oldcwd)
        raise Exception(f'AF structure not found for {pdb}: {e}') from e

    os.chdir(oldcwd)
    return

def find_af_plddt(af_code_full, outfolder="result"):
    '''
    Obtain PLDDT (a measure of certainty where 100 is high and 70 low) value
    for each lysine in an alphafold structure.

    Method
    ------
    If the PLDDT output file hasn't been created yet, create and add column headings
    Open the PLDDT output file for appending data.
    Search the curated AF structure pdb file to extract the PLDDT values.
    Append the PLDDT value to the output file and write to the console log.


    Parameters
    ----------
    af_code_full : string
        The AF code for the structure to extract the PLDDT values from
    
    outfolder : string
        The outdirectory to used to check that the PLDDT out file is in the correct place
        and allow appending to this

    Returns
    -------
    dict_plddt : dictionary
        A dictionary matching up all the lysines with their corresponding PLDDT values
        for the given AF structure.

    Example
    -------
    >>find_af_plddt('AF-P0CG48-F1-model_v4', outfolder='test_plddt')
    > Finding plddt
    AF-P0CG48-F1-model_v4; Resid No. A6; PLDDT: 93.79
    AF-P0CG48-F1-model_v4; Resid No. A11; PLDDT: 89.45
    AF-P0CG48-F1-model_v4; Resid No. A27; PLDDT: 94.28
    ...
    '''

    if not os.path.isfile(os.path.join(outfolder, "curated", "AF_PLDDT_Output.csv")):
        with open(os.path.join(outfolder, "curated", "AF_PLDDT_Output.csv"), 'w', newline='') as plddt_out_file:
            plddt_writer = csv.writer(plddt_out_file)
            plddt_writer.writerow(["Uniprot_Entry", "Resid", "PLDDT"])

    plddt_out_file = open(os.path.join(outfolder, "curated", "AF_PLDDT_Output.csv"), 'a', newline='')
    plddt_writer = csv.writer(plddt_out_file)


    print('> Finding plddt')

    # open .pdb file in assembled folder
    # columns = ['resid', 'chain', 'plddt']
    dict_plddt = dict()

    try:
        f = open(os.path.join(outfolder, "curated", af_code_full + ".pdb"), "r")

        # parse the file to find plddt value.
        for line in f:
            try:
                if re.search('CA  LYS', line):
                    strline = str(line)
                    data = strline.split()
                    if len(data[4]) > 1:
                        resid = (data[4])[1:]
                        plddt = data[9]
                        chain = (data[4])[0]
                    elif len(data[4]) == 1:
                        data = strline.split()
                        resid = data[5]
                        plddt = data[10]
                        chain = data[4]
                        chain_resid = chain + resid

                    # append to dictionary which is later merged into the main dataframe.
                    dict_plddt.update({chain_resid: plddt})

                    # print out what has been found for checking while the output is running
                    print(af_code_full + "; Resid No. " + str(chain_resid) + "; PLDDT: " + plddt)

                    # write the same output to the output file
                    plddt_writer.writerow([af_code_full, str(chain_resid), str(plddt)])

            except Exception as e:
                print(f"Error {e}")
                plddt_writer.writerow([af_code_full, str(chain_resid), f"Error {e}"])
                continue


    except Exception as e:
        print(f'Failed to obtain PLDDT data for {af_code_full}; error: {e}')
        plddt_writer.writerow([af_code_full, str(chain_resid), f"Error {e}"])
        f.close()
        return dict_plddt

    f.close()
    plddt_out_file.close()
    return dict_plddt
