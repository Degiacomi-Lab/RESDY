import os, sys, re
import subprocess

def download_AF_struc(pdb):
    '''
    download AlphaFold2 structures into the assembled folder.
    '''
    
    oldcwd = os.getcwd()
    if not os.path.exists("assembled"):
        os.mkdir("assembled")
    os.chdir('assembled')
    
    print("> Downloading AlphaFold structure")
    try:
        if sys.platform == "win32":
           line = "curl -s -o %s.pdb https://alphafold.ebi.ac.uk/files/%s.pdb"%(pdb, pdb)
        else:
            line = "wget https://alphafold.ebi.ac.uk/files/" + pdb + ".pdb"
        
        subprocess.check_call(line, shell=True)
        
    except Exception:
        print('AF structure not found for %s'%pdb)
        pass
    
    os.chdir(oldcwd)
    #newcwd = os.getcwd

    return

def find_AF_plddt(AF_code_full):
    '''
    Obtain PLDDT (a measure of certainty where 100 is high and 70 low) value
    for each lysine in an alphafold structure.
    '''
    
    print('> Finding plddt')
    
    #Open .pdb file in assembled folder
    #columns = ['resid', 'chain', 'plddt']
    dict_plddt = dict()

    try:
        f = open(os.path.join("assembled", AF_code_full), "r")
 
        #Parse the file to find plddt value.
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
                        
                    #Append to dictionary which is later merged into the main dataframe.
                    dict_plddt.update({chain_resid: plddt})

            except Exception as e:
                print("Error %s"%e)
                continue

    except Exception as e:
        print("ERROR: %s"%e)
        print('Failed to obtain pLDDT data for ' + AF_code_full)
        f.close()
        return dict_plddt

    f.close()
    return dict_plddt
