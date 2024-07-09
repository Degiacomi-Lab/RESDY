import csv
import os, sys, re
import subprocess
from helper import get_download_tool
from csv import writer

def download_AF_struc(pdb, outfolder="result"):
    '''
    download AlphaFold2 structures into the assembled folder.
    '''
    
    oldcwd = os.getcwd()
    mypath = os.path.join(outfolder, "curated")
    if not os.path.exists(mypath):
        os.makedirs(mypath)
        
    os.chdir(mypath)

    print("> downloading AlphaFold structure %s"%pdb)

    tool = get_download_tool()
    try:
        if tool == "curl":
           line = "curl -s -o %s.pdb https://alphafold.ebi.ac.uk/files/%s.pdb"%(pdb, pdb)
        elif tool == "wget":
            line = "wget https://alphafold.ebi.ac.uk/files/" + pdb + ".pdb"
        else:
            raise RuntimeError("You don't have a commandline tool for downloading files")
        
        subprocess.check_call(line, shell=True)
        
    except Exception:
        os.chdir(oldcwd)
        raise Exception('AF structure not found for %s'%pdb)
        
    os.chdir(oldcwd)
    return

def find_AF_plddt(AF_code_full, outfolder="result"):
    '''
    Obtain PLDDT (a measure of certainty where 100 is high and 70 low) value
    for each lysine in an alphafold structure.
    '''

    if not os.path.isfile(os.path.join(outfolder, "curated", "AF_PLDDT_Output.csv")):
        with open(os.path.join(outfolder, "curated", "AF_PLDDT_Output.csv"), 'w', newline='') as PLDDT_out_file:
            writer = csv.writer(PLDDT_out_file)
            writer.writerow(["Uniprot_Entry", "Resid", "PLDDT"])
    else:
        PLDDT_out_file = open(os.path.join(outfolder, "curated", "AF_PLDDT_Output.csv"), 'a', newline='')
        writer = csv.writer(PLDDT_out_file)

    
    print('> Finding plddt')
    
    # open .pdb file in assembled folder
    # columns = ['resid', 'chain', 'plddt']
    dict_plddt = dict()

    try:
        f = open(os.path.join(outfolder, "curated", AF_code_full + ".pdb"), "r")
 
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
                    print(AF_code_full + "; Resid No. " + str(chain_resid) + "; PLDDT: " + plddt)

                    # write the same output to the output file
                    writer.writerow([AF_code_full, str(chain_resid), str(plddt)])

            except Exception as e:
                print("Error %s"%e)
                writer.writerow([AF_code_full, str(chain_resid), "Error %s"%e])
                continue
        
        

    except Exception as e:
        print("ERROR: %s"%e)
        print('Failed to obtain pLDDT data for ' + AF_code_full)
        writer.writerow([AF_code_full, str(chain_resid), "Error %s"%e])
        f.close()
        return dict_plddt

    f.close()
    PLDDT_out_file.close()
    return dict_plddt
