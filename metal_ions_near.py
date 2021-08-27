import biobox as bb
import pandas as pd
import numpy as np
import glob
import math




def identity_local_metal_ions(pdb_code):
    list_of_chains = list()
    list_of_resid = list()
    columns = ['PDB Code', 'Lysine Index', 'Local Amino Acid', 'Distance']
    local_aa_df = pd.DataFrame(columns=columns)
    print('BREAKING UP MOLECULE')
    print(pdb_code)

    M = bb.Molecule()
    path = 'assembled/' + pdb_code
    print(path)
    M.import_pdb(path, include_hetatm=True)
    df = M.data
    df.to_csv('1dpm_info.csv')
    print(df)
    lys_coords, lys_idx = M.atomselect('*','LYS', 'NZ', use_resname=True, get_index=True)


    for entry in lys_idx:
        chain = df.at[entry, 'chain']
        list_of_chains.append(chain)

    for entry in lys_idx:
        resid = df.at[entry, 'resid']
        list_of_resid.append(resid)

    all_coords, idx = M.atomselect('*','*','*', get_index=True)


    for j in range(len(lys_coords)):
        try:
            list_close_points = list()
            
            for i in range(len(all_coords)):
                x_dist = (((lys_coords[j])[0] - (all_coords[i])[0])**2)
                y_dist = (((lys_coords[j])[1] - (all_coords[i])[1])**2)
                z_dist = (((lys_coords[j])[2] - (all_coords[i])[2])**2)
                distance = math.sqrt(x_dist + y_dist + z_dist)
                if distance < 10:
                    index_of_aa = idx[i]
                    aa_resid = df.at[index_of_aa, 'resid']
                    lys_resid = list_of_resid[j]

                    if aa_resid == lys_resid:
                        continue
                    else:
                        identity = df.at[index_of_aa, 'resname']
                        d = ({'PDB Code': pdb_code, 'Lysine Index': lys_idx[j], 'Local Amino Acid': identity, 'Distance':distance})
                        local_aa_df = local_aa_df.append(d, ignore_index=True)
                    #print(local_aa_df)
        except Exception as e:
            print('Error %s'%e)
    local_aa_df.to_csv('local_aa_near_lysine.csv')
    return(local_aa_df)
