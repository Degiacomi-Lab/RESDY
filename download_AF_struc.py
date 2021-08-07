import os
import subprocess

def download_AF_struc(pdb):
    oldcwd = os.getcwd()
    if not os.path.exists("assembled"):
        os.mkdir("assembled")
    os.chdir('assembled')
    subprocess.check_call("wget https://alphafold.ebi.ac.uk/files/" + pdb + ".pdb", shell=True)
    os.chdir(oldcwd)
    print('PDB Structure for ' + pdb + ' downloaded')
    return()