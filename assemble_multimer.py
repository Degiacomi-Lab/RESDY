import biobox as bb
import pandas as pd
import os
import fileinput
#get_data.py and the autopatcher break the protein up into chains.
#This code 'reassembles' the protein into a multimer from the chains in the clean folder given the pdb code and chains the protein consists of.
def assemble_multimer(pdb_code, list_chains):
    #First makes sure there is an assembled folder.
    if not os.path.exists("assembled"):
        os.mkdir("assembled")

    #Next open defines the name of the assembly.
    try:
        name_of_assembly = pdb_code + '_assembled.pdb'
        Multi = bb.Multimer()

    #Next opens the pdb file for each chain in turn and appends to Multi
        for chain in list_chains:
            name_of_pdb_file = pdb_code + '_' + chain + '.pdb'
            patched_pdb_file = pdb_code + '_' + chain + '_' + 'patched.pdb'
            M = bb.Molecule()
            try:
                path = 'curate_PDB/clean/' + name_of_pdb_file
                M.import_pdb(path, include_hetatm=True)
                Multi.append(M)
            except:
                try:
                    path = 'curate_PDB/clean/' + patched_pdb_file
                    M.import_pdb(path, include_hetatm=True)
                    Multi.append(M)
                except Exception as e:
                    print("Error: %s"%e)
                    continue

    #Lastly writes out Multi as a .pdb file.

        path = 'assembled/' + name_of_assembly
        Multi.write_pdb(path)
        print('Success assembling ' + pdb_code)
        filename = path

        with fileinput.FileInput(filename, inplace = True) as f:
            for line in f:
                line = line.replace("TER","")
                print(line, end ='') 

    except Exception as e:
        print("Error: %s"%e)
        print('Failed to assemble ' + pdb_code)

    return()






    



    
    #A=bb.Molecule()
    #A.import_pdb('4XBJ_A.pdb')

    #B=bb.Molecule()
    #B.import_pdb('4XBJ_B.pdb')

    #C=bb.Molecule()
    #C.import_pdb('4XBJ_C_patched.pdb')

   #D=bb.Molecule()
    #D.import_pdb('4XBJ_D_patched.pdb')

    #Assem = bb.Assembly()
    #Assem.append(A, label='chain A')
    #Assem.append(B, label='chain B')
    #Assem.append(C, label='chain C')
    #Assem.append(D, label='chain D')
    #Assem.write_pdb('Assembly.pdb')

if __name__ == "__main__":
    try:    
        columns = ['Uniprot Entry', 'PDB Code', 'Method Structure Obtained by', 'Resolution', 'Chains']
        results_df = pd.DataFrame(columns=columns)
        chain = ['A', 'B']
        for i in range(len(chain)):
            data = ({'Uniprot Entry':'P0A434' , 'PDB Code': '1DPM', 'Method Structure Obtained by': 'X-ray', 'Resolution': '2.2', 'Chains': chain[i]})
            results_df = results_df.append(data, ignore_index=True)
        print(results_df)
        assemble_multimer(results_df)
    except Exception as e:
        print("ERROR: %s"%e)