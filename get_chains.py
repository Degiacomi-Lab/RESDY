import biobox as bb
import subprocess
import os

def get_chains(PDBCODE_inpt):
    try:
        subprocess.check_call("wget https://files.rcsb.org/download/" + PDBCODE_inpt + ".pdb", shell=True)
        M = bb.Molecule()
        M.import_pdb(PDBCODE_inpt + '.pdb')
        df = M.data
        chain_column = df['chain'].tolist()
        unique_values_chain = list()
        for i in range(len(chain_column)):
            if chain_column[i] not in unique_values_chain:
                unique_values_chain.append(chain_column[i])

        os.remove(PDBCODE_inpt + '.pdb')
    except:
        print('Failed to obtain chain information for ' + PDBCODE_inpt)

    return(unique_values_chain)