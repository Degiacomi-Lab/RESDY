import biobox as bb
import pandas as pd
import math
import os
from biobox.measures.calculators import sasa
def break_up_and_calculate_sasa(pdb_code):
    try:
        print('BREAKING UP MOLECULE')
        print(pdb_code)
        list_of_index = list()
        list_of_sasa = list()
        list_of_resid = list()
        list_of_chains = list()
        chain_resid_list = list()

        M = bb.Molecule()
        path = 'assembled/' + pdb_code
        print(path)
        M.import_pdb(path)
        df = M.data

        lys_coords, lys_idx = M.atomselect('*','LYS', 'NZ', use_resname=True, get_index=True)
        df = M.data

        for entry in lys_idx:
            chain = df.at[entry, 'chain']
            list_of_chains.append(chain)

        for entry in lys_idx:
            resid = df.at[entry, 'resid']
            list_of_resid.append(resid)
        print(list_of_chains)
        all_coords, idx = M.atomselect('*','*','*', get_index=True)



        for j in range(len(lys_coords)):
            try:
                list_close_points = list()
                
                for i in range(len(all_coords)):
                    x_dist = (((lys_coords[j])[0] - (all_coords[i])[0])**2)
                    y_dist = (((lys_coords[j])[1] - (all_coords[i])[1])**2)
                    z_dist = (((lys_coords[j])[2] - (all_coords[i])[2])**2)
                    distance = math.sqrt(x_dist + y_dist + z_dist)
                    if distance < 15:
                        list_close_points.append(idx[i])
                #print(list_close_points)
                M.write_pdb('teststruc.pdb', index=list_close_points, split_struc=False)
                S = bb.Molecule()
                S.import_pdb('temp_struc.pdb')
                chain = list_of_chains[j]
                resid = list_of_resid[j]
                pts_2, indx_2 = S.atomselect(chain, [resid], ["CB", "CG", "CD", "CE", "NZ"],  use_resname=False, get_index=True)
                x = sasa(S, targets=indx_2, probe=1.4, n_sphere_point=960, threshold=0.05)
                print(x[0])
                chain_resid_list.append(str(chain + str(resid)))
                list_of_sasa.append(x[0])
            except:
                print('Error obtaining SASA for ' + pdb_code + ' index value ' + str(j))
                continue

        df = pd.DataFrame({'Assembled Index': lys_idx, 'chain': list_of_chains, 'resid': list_of_resid, 'sasa': list_of_sasa, 'Chain_Resid':chain_resid_list})
        print(df)
        os.remove('temp_struc.pdb')
    except Exception as e:
        print("Error %s"%e)
        print('Error obtaining sasa data for ' + pdb_code)
        

    return(df)
if __name__ == "__main__":
    try:    
        print(break_up_and_calculate_sasa('4XBJ_assembled.pdb'))
    except Exception as e:
        print("ERROR: %s"%e)