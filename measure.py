import re
import os, shutil
import subprocess
import glob

import pandas as pd
import numpy as np

import biobox as bb


class Measure(object):
    
    def __init__(self, df_input, outfolder="result"):
        
        self.df_input = df_input
        self.folder = os.path.join(outfolder, "curated")
               
        # Check that all files in DataFrame appear at least once in folder
        files1=[os.path.basename(c).split(".")[0] for c in glob.glob(os.path.join(self.folder, "*pdb"))] #to find AlphaFold entries
        files2=[os.path.basename(c).split("-")[0] for c in glob.glob(os.path.join(self.folder, "*pdb"))] #to find PDB entries
        for f in df_input["PDB Code"].values:
            if f not in files1 and f not in files2:
                print("WARNING: %s not found in folder %s"%(f, self.folder))

        self.pkaoutdir = os.path.join(outfolder, "propkaoutput")
        if not os.path.exists(self.pkaoutdir):
            os.makedirs(self.pkaoutdir)

        columns = ['Uniprot Entry', 'PDB Code', 'Method', 'Resolution', 'Chain', 'Resid']
        self.df = pd.DataFrame(columns=columns)

        # measures to carry out [label for DataFrame column, and function evaluating a file]
        # functions must return a dataframe [chain, resid, measure]
        self.measures = [["pKa", self.calculate_pka], ["sasa", self.calculate_sasa]]
        
    
    def measure_dataframe(self):

        files = glob.glob(os.path.join(self.folder, "*pdb"))

        for index, row in self.df_input.iterrows():

            PDBCODE = row["PDB Code"]
            uniprot_code = row["Uniprot Entry"]
            chains = row["Chains"].split("/")
            
            print("\n# UNIPROT: %s PDB: %s, chain(s): %s"%(uniprot_code, PDBCODE, " ".join(chains)))
            
            # calculate features values from all PDB files associated with specific DataFrame entry
            for f in files:
                
                if PDBCODE not in f:
                    continue

                print("\n> File: %s"%f)
  
                # create temporary DataFrame for data of current file,
                # to be then appended to main DataFrame self.df
                columns = ['Uniprot Entry', 'PDB Code', 'Method', 'Resolution', 'Chain', 'Resid']
                df = pd.DataFrame(columns=columns)

                # append to temporary DataFrame all lysines in the file of interest             
                M = bb.Molecule(f)
                _, idxs = M.atomselect("*", ["LYS"], ["CA"], get_index=True, use_resname=True)
                for i in idxs:

                    #save only lysine entries from chain of interest
                    if M.data["chain"].values[i] not in chains:
                        continue
                
                    data = ({'Uniprot Entry': uniprot_code,
                    'PDB Code': f.split(".")[0],
                    'Method': row["Method"],
                    'Resolution': row["Resolution"],
                    'Chain': M.data["chain"].values[i],
                    'Resid': M.data["resid"].values[i]})
   
                    df = pd.concat([df, pd.DataFrame.from_records(data, index=[0])], ignore_index=True)
                   
                print(">> %s lysines of interest found"%len(df))

                # iterate over measures to carry out (according to self.measures)
                for meas in self.measures:
                    print(">> evaluating %s..."%meas[0])
                    try:
                        df[meas[0]] = np.nan # create new column for measure
                        result = meas[1](f) # run measurement
                        df = self._combine_dataframes(df, result, meas[0]) #insert measures into temporary DataFrame
                    except Exception as e:
                        print("ERROR: %s"%e)
                        continue
                
                #append temporary DataFrame with all measures on a single file to main DataFrame
                self.df = pd.concat([self.df, df], ignore_index=True)

                
    def _combine_dataframes(self, target, to_merge, col_name):
        '''
        target is a DataFrame to be filled with data, to_merge contains the data.
        Values to insert are indexed in both array by two columns: Chain and Resid.
        '''
        
        for i in range(len(target)):
            
            chain_value = target.loc[i, "Chain"]
            resid_value = target.loc[i, "Resid"]
            
            idx = np.where((to_merge["Chain"] == chain_value) & (to_merge["Resid"].astype(int) == resid_value))
            if len(idx[0]) == 0:
                continue
            
            target.at[i, col_name] = to_merge.loc[idx[0][0], col_name]
          
        return target

            
    def calculate_pka(self, path):
        '''
        Call PROPKA to calculate the pKa of a file, parse the .pka file to extract lysine data
        parse errors, and return a dataframe containing all measurements not yielding an error.
        '''
        
        code_for_df = os.path.basename(path).split(".")[0]
        error_file_name = os.path.join(self.pkaoutdir, "%s_propka_errors.txt"%code_for_df)
    
        try:   
            f = open(error_file_name, 'w')
            process = subprocess.Popen(['python', '-m', 'propka', path],
                                stdout=f, stderr=f)
            stdout, stderr = process.communicate()
            f.close()
                    
        except Exception as e:
            f.close()
            
            try:
                shutil.move(code_for_df, os.path.join(self.pkaoutdir, code_for_df))
            except:
                pass
            
            raise Exception('Failed to obtain pKa data. %s.'%e)
    
        try:
            propka_lys_fails = self.parse_propka_errors(error_file_name)
        except Exception as e:
            raise Exception("Failed extracting propka errors. %s"%e)
           
        try:
            pkafile = code_for_df + '.pka'
            propres = open(pkafile)
        except Exception:
            raise Exception('Failed to find %s'%pkafile)
                
        lys_number = list()
        pkas = list()
        chain = list() 
        try:
            for line in propres:
                if re.search('^   LYS' , line):
                    try:
    
                        line = line[6:]
                        line = line.split()

                        # reject adding entries associated with errors in structure (as per logfile)
                        if len(propka_lys_fails)>0:
                            idx = np.where((propka_lys_fails["Chain"] == line[1]) & (propka_lys_fails["Resid"].astype(int) == int(line[0])))
                            if len(idx[0])>0:
                                continue    
            
                        lys_number.append(line[0])
                        chain.append(line[1])
                        pkas.append(line[2])
                        
                    except Exception as e:
                        print("> Error %s"%e)
                        continue
    
            propres.close()
            shutil.move(pkafile, os.path.join(self.pkaoutdir, pkafile))
        
        except Exception as e:
            propres.close()
            shutil.move(pkafile, os.path.join(self.pkaoutdir, pkafile))
            raise Exception('Failure parsing %s.pka. %s'%(code_for_df, e))
            
        try:
            df = pd.DataFrame({'Resid':lys_number,
                               'Chain': chain,
                               'pKa':pkas})
            
            df.sort_values(by=['pKa'], inplace=True)
            df = df.dropna()
            df = df.drop_duplicates(subset=None, keep='first', inplace=False, ignore_index=True)
            
        except Exception as e:
            raise Exception('Failed to construct pKa dataframe. %s'%e)
    
        return df
    

    def parse_propka_errors(self, path):
        '''
        parse the PROPKA output file and appends unique chain and resid of any lysines mentioned a DataFrame.
        This list is returned to main and later the residues in it are removed from the df.
        '''
        
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
                if resid != "" and chain != "":
                    list_remove.append([chain, resid])  
    
            lys_raw_2 = re.findall('[\d]*-LYS \(\w\)', line)
    
            for line in lys_raw_2:
                words = line.split()
                chain = (words[1])[1:-1]
                words_2 = line.split('-')
                resid = words_2[0]
                if resid != "" and chain != "":
                    list_remove.append([chain, resid])  
        
        f.close()
        
        # if the file is completely empty, let's just wipe it!
        if cnt == 0:
            os.remove(path)
           
        if len(list_remove) == 0:
            return []
        else:
            return pd.DataFrame(np.array(list_remove), columns=["Chain", "Resid"]).drop_duplicates()
    
    
    def calculate_sasa(self, path):
        '''
        Form small structures which include just the atoms surrounding the lysine of interest,
        and compute the SASA from that.
        '''
        
        try:
            ###print('>>> Breaking up molecule (for SASA calculation)')
            list_of_sasa = list()
            list_of_resid = list()
            list_of_chains = list()
    
            #read PDB file
            M = bb.Molecule()
            M.import_pdb(path, include_hetatm=True)
            df = M.data
    
            #Find the coordinates and index of all lysine residues in the protein.
            lys_coords, lys_idx = M.atomselect('*', ['LYS'], 'NZ', use_resname=True, get_index=True)
            df = M.data
    
            #Use get_subset...
            #Find the chain and resid of each lysine.
            for entry in lys_idx:
                chain = df.at[entry, 'chain']
                list_of_chains.append(chain)
    
            for entry in lys_idx:
                resid = df.at[entry, 'resid']
                list_of_resid.append(resid)
                
            #Find the coordinates and index of every atom in the molecule.
            all_coords, idx = M.atomselect('*','*','*', get_index=True)
    
        except Exception as e:
            raise Exception("%s"%e)
        
        ##print('>>> Obtaining SASA')
        #For each lysine it works out the distance between the lys NZ,
        #and the each atom in the protein.
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
    
            #if the atoms are close to the lys NZ they are included in a small .pdb structure.
            try:
                M.write_pdb('temp_struc.pdb', index=list_close_points, split_struc=False)
    
                S = bb.Molecule()
                S.import_pdb('temp_struc.pdb', include_hetatm=True)
                chain = list_of_chains[j]
                resid = list_of_resid[j]
    
                #SASA is calculated for that lysine in the small molecule.
                pts_2, indx_2 = S.atomselect(chain, [resid], ["CB", "CG", "CD", "CE", "NZ"],  use_resname=False, get_index=True)
                x = bb.sasa(S, targets=indx_2, probe=1.4, n_sphere_point=960, threshold=0)
                list_of_sasa.append(x[0])
    
            except:
                print('Error obtaining SASA at index value ' + str(j))
                list_of_sasa.append(None)
                continue
    
        #append results to a df which is given as output
        try:
            df = pd.DataFrame({'Chain': list_of_chains,
                               'Resid': list_of_resid,
                               'sasa': list_of_sasa})
    
        except Exception as e:
            raise Exception('Error obtaining SASA data. %s'%e)
    
        try:
            os.remove('temp_struc.pdb')
            
        except Exception as e:
            print("Error %s"%e)
            print('Failed to remove temporary pdb structure.')
            pass
        
        return df
        

if __name__ == "__main__":
    pass
    #M = Measure("")