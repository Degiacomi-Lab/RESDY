import biobox as bb
import pandas as pd
import os

def assemble_multimer(results_df):
    try:
        if not os.path.exists("assembled"):
            os.mkdir("assembled")

        list_of_pdb_codes = results_df['PDB Code']
        list_of_pdb_codes_no_dup = list()

        for entry in list_of_pdb_codes:
            if entry not in list_of_pdb_codes_no_dup:
                list_of_pdb_codes_no_dup.append(entry)
            else:
                continue
        
        for pdb_code in list_of_pdb_codes_no_dup:
            df_one_pdb_code = results_df.where(results_df['PDB Code'] == pdb_code)
            name_of_assembly = pdb_code + '_assembled.pdb'
            for i in range(len(df_one_pdb_code)):
                Multi = bb.Multimer()
                list_of_chains = df_one_pdb_code['Chains']
                for chain in list_of_chains:
                        name_of_pdb_file = pdb_code + '_' + chain + '.pdb'
                        patched_pdb_file = pdb_code + '_' + chain + '_' + 'patched.pdb'
                        label = 'Chain ' + chain
                        M = bb.Molecule()
                        try:
                            path = 'curate_PDB/clean/' + name_of_pdb_file
                            M.import_pdb(path)
                            Multi.append(M)
                        except:
                            try:
                                path = 'curate_PDB/clean/' + patched_pdb_file
                                M.import_pdb(path)
                                Multi.append(M)
                            except Exception as e:
                                print("Error: %s"%e)
                                continue
            path = 'assembled/' + name_of_assembly
            Multi.write_pdb(path)
            print('Success assembling ' + pdb_code)
    except Exception as e:
        print("Error: %s"%e)
        print('Error Assembling Multimer')
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
        chain = ['A', 'B', 'C', 'D']
        for i in range(len(chain)):
            data = ({'Uniprot Entry': 'P0A6B4', 'PDB Code': '4XBJ', 'Method Structure Obtained by': 'X-ray', 'Resolution': '2.2', 'Chains': chain[i]})
            results_df = results_df.append(data, ignore_index=True)
        print(results_df)
        assemble_multimer(results_df)
    except Exception as e:
        print("ERROR: %s"%e)