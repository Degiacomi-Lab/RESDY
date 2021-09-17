import re
import os, shutil
import subprocess

import pandas as pd
import numpy as np

import biobox as bb

#call propka to calculate the pKa of all groups in the assembled/*.pdb file
#and then parse the .pka file for the pka of lysine side chains in the protein.
def calculate_pKa(code, outdir="Output"):

    print('> Obtaining pKa')
    code_for_df = code[:4]
    path = os.path.join("assembled", code)
    no_pdb = code[:-4]
    
    pkaoutdir = os.path.join(outdir, "propkaoutput")
    if not os.path.exists(pkaoutdir):
        os.mkdir(pkaoutdir)

    error_file_name = os.path.join(pkaoutdir, "%s_propka_errors.txt"%code.split(".")[0])

    try:   
        f = open(error_file_name, 'w')
        process = subprocess.Popen(['python', '-m', 'propka', path],
                            stdout=f, stderr=f)
        stdout, stderr = process.communicate()
        f.close()
                
    except Exception as e:
        f.close()
        try:
            shutil.move(pkafile, os.path.join(pkaoutdir, pkafile))
        except:
            pass

        raise Exception('Failed to obtain pKa data. %s.'%e)

    
    try:
        propka_lys_fails = parse_propka_errors(error_file_name)
    except:
        raise Exception("Failed extracting propka errors. %s"%e)
        
    try:
        pkafile = code_for_df + '_assembled.pka'
        propres = open(pkafile)
        AF_struc = False
    except:
        try:
            code_pka = (code[:-4]) + '.pka'
            pkafile = code_pka
            propres = open(pkafile)
            AF_struc = True
        except Exception as e:
            raise Exception('Failed to find ' + code + '.pka')
            
    lys_number = list()
    pkas = list()
    chain = list()
    chain_resid_list = list()
    buried_percentage = list()
    
    try:
        for line in propres:
            if re.search('^   LYS' , line):
                try:

                    line = line[6:]
                    line = line.split()
                    
                    chain_resid_list.append(line[1] + line[0])
                    lys_number.append(line[0])
                    chain.append(line[1])
                    pkas.append(line[2])
                    
                except Exception as e:
                    print("Error %s"%e)
                    continue

        propres.close()
            
    except Exception as e:
        propres.close()
        shutil.move(pkafile, os.path.join(pkaoutdir, pkafile))
        raise Exception('Failure parsing ' + code + '.pka')

    if AF_struc == True:
        code_for_df = code[:-4]
    elif AF_struc == False:
        code_for_df = code_for_df

    try:
        df = pd.DataFrame({'resid':lys_number, 'chain':chain, 'pKa':pkas, 'PDB Code':code_for_df, 'Chain_Resid': chain_resid_list})
        df.sort_values(by=['pKa'], inplace=True)
        convert_dict = {'pKa': float}
        df = df.astype(convert_dict)
        convert_dict = {'chain': float}
        df = df.where(df['pKa'] < 20)
        df = df.dropna()
        df = df.drop_duplicates(subset=None, keep='first', inplace=False, ignore_index=True)
        
    except Exception as e:
        raise Exception('Failed to construct pKa dataframe. %s'%e)

    finally:
        try:
            shutil.move(pkafile, os.path.join(pkaoutdir, pkafile))
        except:
            print("Warning: could not move file %s to destination directory %s"%(pkafile, pkaoutdir))
            pass
        
    return df, propka_lys_fails


#Forms small structures which include just the atoms surrounding the lysine of interest
#and computes SASA from that.
def break_up_and_calculate_sasa(pdb_code):

    try:
        print('> Breaking up molecule (for SASA calculation)')
        list_of_index = list()
        list_of_sasa = list()
        list_of_resid = list()
        list_of_chains = list()
        chain_resid_list = list()

        #Accesses .pdb file from assembled/
        M = bb.Molecule()
        path = os.path.join("assembled", pdb_code)
        #path = "assembled%s%s"%(os.sep, pdb_code)
        M.import_pdb(path, include_hetatm=True)
        df = M.data

        #Finds the coordinates and index of all lysine residues in the protein.
        lys_coords, lys_idx = M.atomselect('*', ['LYS'], 'NZ', use_resname=True, get_index=True)
        df = M.data

        #Use get_subset...
        #Finds the chain and resid of each lysine.
        for entry in lys_idx:
            chain = df.at[entry, 'chain']
            list_of_chains.append(chain)

        for entry in lys_idx:
            resid = df.at[entry, 'resid']
            list_of_resid.append(resid)
            
        #Finds the coordinates and index of every atom in the molecule.
        all_coords, idx = M.atomselect('*','*','*', get_index=True)

    except Exception as e:
        raise Exception("%s"%e)
    
    print('> Obtaining SASA')
    #For each lysine it works out the distance between the lys NZ and the each atom in the protein.
    for j in range(len(lys_coords)):
        list_close_points = list()
        
        for i in range(len(all_coords)):
            try:
                x_dist = (((lys_coords[j])[0] - (all_coords[i])[0])**2)
                y_dist = (((lys_coords[j])[1] - (all_coords[i])[1])**2)
                z_dist = (((lys_coords[j])[2] - (all_coords[i])[2])**2)
                distance = np.sqrt(x_dist + y_dist + z_dist)
                if distance < 15:
                    list_close_points.append(idx[i])
            except:
                continue

        #If the atoms are close to the lys NZ they are included in a small .pdb structure.
        try:
            M.write_pdb('temp_struc.pdb', index=list_close_points, split_struc=False)

            S = bb.Molecule()
            S.import_pdb('temp_struc.pdb', include_hetatm=True)
            chain = list_of_chains[j]
            resid = list_of_resid[j]

            #SASA is calculated for that lysine in the small molecule.
            pts_2, indx_2 = S.atomselect(chain, [resid], ["CB", "CG", "CD", "CE", "NZ"],  use_resname=False, get_index=True)
            x = bb.sasa(S, targets=indx_2, probe=1.4, n_sphere_point=960, threshold=0)
            chain_resid_list.append(str(chain + str(resid)))
            list_of_sasa.append(x[0])

        except:
            print('Error obtaining SASA for ' + pdb_code + ' index value ' + str(j))
            list_of_sasa.append(None)
            continue

    #The results are appended to a df which is given as outpit
    try:
        df = pd.DataFrame({'Assembled Index': lys_idx, 'chain': list_of_chains, 'resid': list_of_resid, 'sasa': list_of_sasa, 'Chain_Resid':chain_resid_list})

    except Exception as e:
        raise Exception('Error obtaining SASA data. %s'%e)

    try:
        os.remove('temp_struc.pdb')
        
    except Exception as e:
        print("Error %s"%e)
        print('Failed to remove temporary pdb structure for ' + pdb_code)
        pass
        
    return df


#This function parses the propka output file and appends the chain_resid of any lysines mentioned into list_remove.
#This list is returned to main and later the residues in it are removed from the df.
def parse_propka_errors(path):
    
    print('> Checking for errors and warnings in PROPKA')
    f = open(path, 'r')
    list_remove = list()
    cnt = 0
    for line in f:
        cnt += 1
        lys_raw = re.findall('LYS [\d]*[\s][\w]*', line)

        for line in lys_raw:
            words = line.split(' ')
            resid = words[1]
            chain = words[2]
            chain_resid = chain + resid
            list_remove.append(chain_resid)

        lys_raw_2 = re.findall('[\d]*-LYS \(\w\)', line)

        for line in lys_raw_2:
            words = line.split(' ')
            chain = (words[1])[1:-1]
            words_2 = line.split('-')
            resid = words_2[0]
            chain_resid = chain + resid
            list_remove.append(chain_resid)
    
    
    f.close()
    
    # if the file is completely empty, it is pointless to keep it. Let's wipe it!
    if cnt == 0:
        os.remove(path)
    
    
    return list_remove


#if __name__ == "__main__":
