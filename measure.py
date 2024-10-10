import re
import os
import logging
import datetime
import shutil
import subprocess
import glob
import time
import pandas as pd
import numpy as np
import biobox as bb

# AEV packages
try:
    from ase import Atoms
    import torch
    import torchani
except Exception as e:
    print(f'Packages required for AEV calculation are not available, will not be able to calculate AEVs. Error: {e}')

try:
    from Bio.PDB import PDBParser
    from Bio.PDB.ResidueDepth import min_dist, get_surface, residue_depth
except Exception as e:
    print(f"biopython and msms unavailable. Unable be able to calculate residue depth. Error: {e}")


class Measure(object):

    def __init__(self, df_input, outdir="result", activate_log=False, log_path='measure_log.txt',
                 features=['propka', 'pkaANI', 'sasa', 'depth', 'aev', 'das']):

        self.activate_log = False
        if activate_log:
            self.activate_log = activate_log
            self.log_path = os.path.join(outdir, log_path)

            # define a logger
            self.logger = logging.getLogger('MeasureLog')
            self.logger.setLevel(level = logging.DEBUG)

            formatter = logging.Formatter('%(message)s') # as simple as possible
            handler = logging.FileHandler(self.log_path, encoding = 'UTF-8')
            handler.setLevel(logging.INFO)
            handler.setFormatter(formatter)

            self.logger.addHandler(handler)

        self._setup_measures(features)

        # for restarting
        self.current_index = 0

        # document failed pdb files
        self.wrong_pdb_file = []

        self.outdir = outdir
        self.df_input = df_input
        self.folder = os.path.join(outdir, "curated")

        # Check that all files in DataFrame appear at least once in folder
        # find all AlphaFold entries
        files_af=[os.path.basename(c).split(".")[0] for c in glob.glob(os.path.join(self.folder, "*pdb"))]
        # find all PDB entries
        files_pdb=[os.path.basename(c).split("-")[0] for c in glob.glob(os.path.join(self.folder, "*pdb"))]
        for f in df_input["PDB Code"].values:
            if f not in files_af and f not in files_pdb:
                print(f'WARNING: {f} not found in folder {self.folder}')

        self.pkaoutdir = os.path.join(outdir, "propkaoutput")
        if not os.path.exists(self.pkaoutdir):
            os.makedirs(self.pkaoutdir)

        self.PDB_only = False

        if 'Uniprot_Entry' in self.df_input.columns: 
            columns = ['Uniprot_Entry', 'PDB Code', 'Method', 'Resolution', 'Chain', 'Resid']
            self.df = pd.DataFrame(columns=columns)
        else:
            self.PDB_only = True
            columns = ['PDB Code', 'Chain', 'Resid', 'propka', 'pkaANI', 'sasa', 'depth']
            self.df = pd.DataFrame(columns = columns)

    def _setup_measures(self, features):
        '''
        convert a list of features into a measuring protocol
        '''

        # measures to carry out [label for DataFrame column, and function evaluating a file]
        # functions must return a dataframe [chain, resid, measure]
        self.measures = []
        for m in features:
            if m == "propka":
                self.measures.append([m, self.calculate_pka_propka])
            elif m == "pkaANI":
                self.measures.append([m, self.calculate_pkaANI])
            elif m == "sasa":
                self.measures.append([m, self.calculate_sasa])
            elif m == "depth":
                self.measures.append([m, self.calculate_depth])
            elif m == 'aev':
                self.measures.append([m, self.calculate_aevs])
            elif m == 'das':
                self.measures.append([m, self.calculate_das])
            else:
                raise Exception(f"measure {m} unknown")

    def save_state(self, outname="measures.csv"):
        '''
        Function saves a csv file of all of the measurements calculated through measure_dataframe()
        File automatically saved in the output directory that has been set
         previously when setting up the measures class
        Option to customise the name of the output file through outname parameter

        Parameters
        ----------
        outname : string
            the name of the csv file that the output is written to

        Example
        -------
        M.save_state(outname='measures.csv')
        '''
        self.df.to_csv(os.path.join(self.outdir, outname), index_label=False, index=False)


    def measure_dataframe(self):
        '''
        DOCSTRING WRITING NEEDED

        Method
        ------

        Parameters
        ----------

        Returns
        -------

        Example
        -------

        '''
        # use a different method if handling pdb codes only
        if self.PDB_only:
            return 'Call PDB_only method instead'

        files = glob.glob(os.path.join(self.folder, "*pdb"))
        # remove the files which have pkaani in the name as these are output files from pkaani
        files = [file for file in files if 'pkaani' not in file]

        first_index = self.df_input.index[0]
        total_structures = len(files)
        current_structure = 0
        overall_st = time.time()
        print(f'Total number of structures to analyse: {total_structures}')

        for index, row in self.df_input.iterrows():
            pdb_code = row["PDB Code"]
            chains = row["Chains"].split("/")
            self.current_index = index

            uniprot_code = row["Uniprot_Entry"]
            print(f'\n# UNIPROT: {uniprot_code} PDB: {pdb_code}, chain(s): {" ".join(chains)}')

            # calculate features values from all PDB files associated with specific DataFrame entry
            for f in files:

                if pdb_code not in f:
                    continue
                print(f'Analysing structure {current_structure}/{total_structures}')
                tstart = time.time()
                print("\n> File: %s"%f)

                # create temporary DataFrame for data of current file,
                # to be then appended to main DataFrame self.df

                columns = ['Uniprot_Entry', 'PDB Code', 'Method', 'Resolution', 'Chain', 'Resid']
                df_currentfile = pd.DataFrame(columns=columns)

                # append to temporary DataFrame all lysines in the file of interest
                try:
                    M = bb.Molecule(f) # sometimes bb does not work with a pdb file
                except:
                    self.wrong_pdb_file.append(f)
                    continue
                
                _, idxs = M.atomselect("*", ["LYS"], ["CA"], get_index=True, use_resname=True)
                for i in idxs:

                    #save only lysine entries from chain of interest
                    if M.data["chain"].values[i] not in chains:
                        continue
                    
                    
                    data = ({'Uniprot_Entry': uniprot_code,
                        'PDB Code': f.split(".")[0],
                        'Method': row["Method"],
                        'Resolution': row["Resolution"],
                        'Chain': M.data["chain"].values[i],
                        'Resid': M.data["resid"].values[i]})

                    df_currentfile = pd.concat([df_currentfile, pd.DataFrame.from_records(data, index=[0])], ignore_index=True)

                print(f">> {len(df_currentfile)} lysines of interest found")

                # iterate over measures to carry out (according to self.measures)
                for meas in self.measures:
                    print(f">> evaluating {meas[0]}...")
                    try:
                        df_currentfile[meas[0]] = np.nan # create new column for measure
                        result = meas[1](f) # run measurement
                        df_currentfile = self._combine_dataframes(df_currentfile, result, meas[0]) #insert measures into temporary DataFrame

                    except Exception as e:
                        print(f"ERROR: {e}")
                        continue

                current_structure += 1
                print(">> file processed in %4.2f seconds."%(time.time()-tstart))
                average_time_per_file = round(((time.time()- overall_st) / current_structure), 2)
                print(f'>> Time average per file: {average_time_per_file} seconds.')
                sec_remaining = average_time_per_file * (total_structures + 1 - current_structure)
                time_remaining_str = str(datetime.timedelta(seconds=sec_remaining))
                print(f'Predicted time remaining: {time_remaining_str}')

                # document the data to a log file
                if self.activate_log:
                    if df_currentfile.empty is False:
                        try:
                            self.logger.info(df_currentfile.to_string().strip('    Uniprot_Entry                    PDB Code Method Resolution Chain Resid    pKa       sasa'))
                            self.logger.info('--------------------------------------------------------------------------')
                        except Exception as e:
                            print(f'Error in logging: {e}')

                #append temporary DataFrame with all measures on a single file to main DataFrame
                self.df = pd.concat([self.df, df_currentfile], ignore_index=True)

            # remove possible duplicated rows (if restarted)
            if index == first_index:
                try:
                    self.df = self.df.drop_duplicates(subset=None, keep='first', inplace=False, ignore_index=True)
                except Exception as e_one:
                    try:
                        self.df = self.df.loc[self.df.astype(str).drop_duplicates(subset=None, keep='first', inplace=False, ignore_index=True)]
                    except Exception as e_two:
                        print(f'Failed to remove duplicates from measurement dataframe: {e_two}')
                        pass


    def recover_from_log(self, log_path):
        '''
        DOCSTRING WRITING NEEDED

        Method
        ------

        Parameters
        ----------

        Returns
        -------

        Example
        -------

        '''
        if self.PDB_only:
            return 'Function not callable.'

        columns = ['Uniprot_Entry', 'PDB Code', 'Method', 'Resolution', 'Chain', 'Resid', 'pKa', 'sasa']
        log_to_df = pd.DataFrame(columns=columns)

        with open(log_path) as inf:
            for line in inf:
                line = line.replace('--------------------------------------------------------------------------',' ')
                parts = line.split()
                if len(parts) == 0:
                    continue
                parts = parts[1:]
                data = dict(zip(columns, parts))
                log_to_df = log_to_df.append(data, ignore_index=True)

        return log_to_df


    def restart_measure(self, log_path='measure_log.txt'):
        '''
        A function to restart the measurements calculations
        Useful if the initial run of the measurements crashes or gets stuck
        Works out how far along the simulation was by running an analysis of the measures log file

        Method
        ------
        Read over the measures log file and collate a list of files that have been analysed
        Remove the final value from the list as this may not have been done properly
        Remove completed files from files to do
        Restart measure_dataframe() with the new list

        Parameters
        ----------
        log_path : string
            The name of the measures log file
            By default takes the name 'measures_log.txt'

        Example
        -------
        >> M.restart_measure()
        '''
        # TODO is there a way to restart the PDB_only measurements? Might need to produce on if not
        if self.PDB_only:
            return 'restart_measure() function not callable when using PDB_only'

        print(len(self.df_input))

        # 1. Analyse the measures log file to create a list of files that were analysed
        log_path = os.path.join(self.outdir, log_path)
        proteins_completed = []
        words_to_ignore = ['sasa', 'depth', 'pKa', 'PDB', 'Code', 'Method', 'Chain', 'Resid', 'Resolution', 'Uniprot', 'Entry']
        with open(file=log_path, mode='r') as lpf:
            for line in lpf:
                line = line.replace('--------------------------------------------------------------------------',' ')
                parts = line.split()
                if len(parts) == 0:
                    continue
                protein_code = parts[1].split('/')[-1]
                if protein_code not in words_to_ignore:
                    proteins_completed.append(protein_code)

        # remove the last protein from list incase it wasn't completed fully
        final_protein = proteins_completed[-1]
        proteins_completed = [c for c in proteins_completed if c != final_protein]
        proteins_completed = list(set(proteins_completed))

        # 2. Update df_input to only have the files which haven't been analysed yet
        idx_to_remove = []
        for i, r in self.df_input.iterrows():
            if r['Uniprot_Entry'] in proteins_completed:
                idx_to_remove.append(i)
        self.df_input = self.df_input.drop(idx_to_remove)
        print(len(self.df_input))
        # 3. Restart the measure_dataframe() with the new file list
        self.measure_dataframe()


    def _combine_dataframes(self, target, to_merge, col_name):
        '''
        target is a DataFrame to be filled with data, to_merge contains the data.
        Values to insert are indexed in both array by two columns: Chain and Resid.
        '''
        # e.g. self._combine_dataframes(df, result, meas[0])

        for i in range(len(target)):

            chain_value = target.loc[i, "Chain"]
            resid_value = target.loc[i, "Resid"]

            idx = np.where((to_merge["Chain"] == chain_value) & (to_merge["Resid"].astype(int) == resid_value))
            if len(idx[0]) == 0:
                continue

            # this if statement allows you to add the lists of the aevs into the overall dataframe
            if col_name == 'aev':
                target['aev'] = target['aev'].astype('object')

            target.at[i, col_name] = to_merge.loc[idx[0][0], col_name]

        return target


    def measure_PDB_only(self):

        if not self.PDB_only:
            return 'Calling the wrong method.'

        files = glob.glob(os.path.join(self.folder, "*pdb"))

        for _, row in self.df_input.iterrows():
            PDBCODE = row['PDB Code']


            # calculate features values from all PDB files associated with specific DataFrame entry
            for f in files:

                if PDBCODE not in f:
                    continue

                tstart = time.time()
                print(f"\n> File: {f}")

                # create temporary DataFrame for data of current file,
                # to be then appended to main DataFrame self.df

                columns = ['PDB Code', 'Chain', 'Resid']
                df = pd.DataFrame(columns=columns)


                # append to temporary DataFrame all lysines in the file of interest
                try:
                    M = bb.Molecule(f) # sometimes bb does not work with a pdb file
                except Exception as e:
                    self.wrong_pdb_file.append(f)
                    continue

                _, idxs = M.atomselect("*", ["LYS"], ["CA"], get_index=True, use_resname=True)
                for i in idxs:

                    #save only lysine entries from chain of interest
                    #if M.data["chain"].values[i] not in chains:
                    #    continue


                    data = ({'PDB Code': f.split(".")[0],
                        'Chain': M.data["chain"].values[i],
                        'Resid': M.data["resid"].values[i]})

                    df = pd.concat([df, pd.DataFrame.from_records(data, index=[0])], ignore_index=True)

                print(f">> {len(df)} lysines of interest found")

                # iterate over measures to carry out (according to self.measures)
                for meas in self.measures:
                    print(f">> evaluating {meas[0]}...")
                    try:
                        df[meas[0]] = np.nan # create new column for measure
                        result = meas[1](f) # run measurement
                        df = self._combine_dataframes(df, result, meas[0]) #insert measures into temporary DataFrame

                    except Exception as e:
                        print(f"ERROR: {e}")
                        continue

                processing_time = round((time.time()-tstart), 2)
                print(f">> file processed in {processing_time} sec.")

                # document the data to a log file
                if self.activate_log:
                    if df.empty is False:
                        try:
                            self.logger.info(df.to_string().strip('PDB Code        Chain Resid    pKa       sasa'))
                            self.logger.info('--------------------------------------------------------------------------')
                        except Exception as e:
                            print(f'Error in logging measurements: {e}')

                #append temporary DataFrame with all measures on a single file to main DataFrame
                self.df = pd.concat([self.df, df], ignore_index=True)



            '''
            for f in files:
                if PDBCODE not in f:
                    continue
                    
                tstart = time.time()
                print("\n> File: %s"%f)
                
                result_list = list()
                for meas in self.measures:
                    print(">> evaluating %s..."%meas[0])
                    try:
                       
                        result = meas[1](f) # run measurement
                        result_list.append(result)
                        
                    except Exception as e:
                        print("ERROR: %s"%e)
                        continue
                try:
                    # add sasa data to pka data
                    to_merge = result_list[1]
                    target = result_list[0]
                    
                    for i in range(len(to_merge)):      
    
                        chain_value = to_merge.loc[i, "Chain"]
                        resid_value = to_merge.loc[i, "Resid"]
    
                        idx = np.where((target["Chain"] == chain_value) & (target["Resid"] == resid_value))
    
                        if len(idx[0]) == 0:
                            row = to_merge.loc[i,:]
                            target = target.append(row, ignore_index = True)
                            continue
    
                        target.at[idx[0][0], 'sasa'] = to_merge.loc[i, 'sasa']
                   
                    # final thing to do: make sure the order of the columns is correct
                    target = target[['Chain', 'Resid', 'pKa', 'sasa']]
                    print(">> %s lysines of interest found"%len(target))
                    
                    # todo1: add a col i.e. PDB Code
                    target.insert(0, 'PDB Code', [PDBCODE] * len(target))
                    
                    # todo2: append to the main df
                    self.df = pd.concat([self.df, target], ignore_index=True)
                    
                except Exception as e:
                    print("ERROR: %s"%e)
                print(">> file processed in %4.2f sec."%(time.time()-tstart))   
                '''


    def calculate_pka_propka(self, path):
        '''
        Call PROPKA to calculate the pKa of a file, parse the .pka file to extract lysine data
        parse errors, and return a dataframe containing all measurements not yielding an error.
        '''

        code_for_df = os.path.basename(path).split(".")[0]
        error_file_name = os.path.join(self.pkaoutdir, f"{code_for_df}_propka_errors.txt")

        test_path = code_for_df + '.pka'
        if not os.path.isfile(test_path):
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

                raise Exception(f'Failed to obtain pKa data (PROPKA): {e}') from e

        try:
            propka_lys_fails = self.parse_propka_errors(error_file_name)
        except Exception as e:
            raise Exception(f"Failed extracting PROPKA errors from output file. {e}")

        try:
            pkafile = code_for_df + '.pka'
            propres = open(pkafile)
        except Exception as e:
            raise Exception(f'Failed to find {pkafile} output file to read') from e

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

                        lys_number.append(int(line[0]))
                        chain.append(line[1])
                        pkas.append(float(line[2]))

                    except Exception as e:
                        print(f"> Error {e}")
                        continue

            propres.close()
            shutil.move(pkafile, os.path.join(self.pkaoutdir, pkafile))

        except Exception as e:
            propres.close()
            shutil.move(pkafile, os.path.join(self.pkaoutdir, pkafile))
            raise Exception(f'Failure parsing {code_for_df}.pka, error: {e}') from e

        try:
            df_propka = pd.DataFrame({'Resid':lys_number,
                               'Chain': chain,
                               'propka':pkas})

            df_propka.sort_values(by=['propka'], inplace=True)
            df_propka = df_propka.dropna()
            df_propka = df_propka.drop_duplicates(subset=None, keep='first', inplace=False, ignore_index=True)

        except Exception as e:
            raise Exception(f'Failed to construct pKa (PROPKA) dataframe. {e}') from e

        return df_propka


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


    def calculate_pkaANI(self, path):

        code_for_df = os.path.basename(path).split(".")[0]
        pdb_path = path.split(".")[0]
        test_path = pdb_path + '_pka.log'
        if not os.path.isfile(test_path):
            try:
                _ = subprocess.run(['pkaani', '-i', path])
            except Exception as e:
                print(e)
                if "[Errno 2] No such file or directory: 'pkaani'" == str(e):
                    raise Exception(f'Error: Failed to obtain pkaANI data: {e}. Is pkaANI installed correctly? If so, reload environment and try again.')
                else:
                    raise Exception(f"Failed to obtain pkaANI data through running pkaANI, error: {e}")

        try:
            log_file = pdb_path + '_pka.log'
            propres = open(log_file)
        except Exception as e:
            raise Exception(f'Failed to find pkaANI log file: {e}') from e

        lys_number = []
        pkas = []
        chains = []
        try:
            for line in propres:
                if re.search('^LYS', line):
                    line = line[4:]
                    info = line.split()

                    lys_number.append(int(info[0]))
                    chains.append(info[1])
                    pkas.append(float(info[2]))

            propres.close()
        except Exception as e:
            propres.close()
            raise Exception(f'Failure parsing pkaANI log file for: {code_for_df}_pka.log. {e}') from e

        try:
            df_pkaani = pd.DataFrame({'Resid':lys_number,
                         'Chain': chains,
                         'pkaANI':pkas})

            df_pkaani.sort_values(by=['pkaANI'], inplace=True)
            df_pkaani = df_pkaani.dropna()
            df_pkaani = df_pkaani.drop_duplicates(subset=None, keep='first', inplace=False, ignore_index=True)
        except Exception as e:
            raise Exception(f'Failed to construct pkaANI dataframe, error: {e}') from e

        return df_pkaani


    def calculate_sasa(self, path):
        '''
        Calculate the solvent accessible surface area of the NZ atom within the lysine structure

        Method
        ------
        Form small structures which include just the atoms surrounding the lysine of interest.
        Small structures are classified as any atoms within 15 angstroms of the NZ of the lysines.
        A new biobox moleucle is created for the substructure and SASA is calculated from that.
        The SASA calculation uses the bb.sasa() function.

        Parameters
        ----------
        path : string
            The path of the pdb file that SASA is being calculated for.

        Returns
        -------

        '''

        try:
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

            #Find the chain and resid number of each lysine.
            list_of_resid = list(M.data['resid'][lys_idx])
            list_of_chains = list(M.data['chain'][lys_idx])

            #Find the coordinates and index of every atom in the molecule.
            all_coords, idx = M.atomselect('*','*','*', get_index=True)

        except Exception as e:
            raise Exception(f'{e}') from e

        #For each lysine it works out the distance between the lys NZ,
        #and the each atom in the protein.
        for j, lys_coord in enumerate(lys_coords):
            list_close_points = list()

            for i, coord in enumerate(all_coords):
                try:
                    x_dist = (lys_coord[0] - coord[0])**2
                    y_dist = (lys_coord[1] - coord[1])**2
                    z_dist = (lys_coord[2] - coord[2])**2
                    distance = np.sqrt(x_dist + y_dist + z_dist)
                    if distance < 15:
                        list_close_points.append(idx[i])
                except:
                    continue

            #if the atoms are close to the lys NZ they are included in a small .pdb structure.
            try:
                S = M.get_subset(idxs=list_close_points) 
                chain = list_of_chains[j]
                resid = list_of_resid[j]
                #print([chain, resid])

                #SASA is calculated for that lysine in the small molecule.
                pts_2, indx_2 = S.atomselect(chain, [resid], ["CB", "CG", "CD", "CE", "NZ"], 
                                            use_resname=False, get_index=True)
                #print([pts_2, indx_2, S.data['radius']])
                x = bb.sasa(S, targets=indx_2, probe=1.4, n_sphere_point=960, threshold=0)
                list_of_sasa.append(x[0])

            except:
                print('Error obtaining SASA at index value ' + str(j))
                list_of_sasa.append(None)
                continue

        #append results to a df which is given as output
        try:
            df_sasa = pd.DataFrame({'Chain': list_of_chains,
                                'Resid': list_of_resid,
                                'sasa': list_of_sasa})

        except Exception as e:
            raise Exception(f'Error obtaining SASA data. {e}')

        return df_sasa


    def calculate_depth(self, path):

        try:
            M = bb.Molecule(path)
            pos, idx = M.atomselect("*", "*", "NZ", get_index=True)
        except Exception as e:
            raise Exception(f">> could not find NZ atoms within atomic structure - {e}")

        try:
            parser = PDBParser()
            structure = parser.get_structure('structure', path)
            surface = get_surface(structure[0])
        except Exception as e:
            raise Exception(f">> could not get biopython structure - {e}")


        results = []
        for i in range(len(pos)):
            chain = M.data.loc[idx[i], ["chain"]].values[0]
            resid = M.data.loc[idx[i], ["resid"]].values[0]

            mychain = structure[0][chain]
            myres = mychain[int(resid)]

            try:
                #dist = min_dist(pos[i], surface)
                rd = residue_depth(myres, surface)
            except Exception as e:
                raise Exception(f">> failed getting min_dist - {e}")

            results.append([chain, resid, rd])

        df_depth = pd.DataFrame(results, columns=["Chain", "Resid", "depth"])

        return df_depth


    def calculate_aevs(self, path):
        '''
        Calculate the Atomic Environment Vectors (AEVs) of the NZ atom within the lysine structure

        Method
        ------
        Uses the ANI-2x AEV calculator to calculate the AEVs
        Option available to use cuaev accelerated AEV calculation, can also just be run with a cpu
        For each NZ atom withing the lysines of the protein, a substructure is created 
            including all atoms within a cutoff distance
        The cutoff distance is set at 6A currently as this was the minimum distance needed
            for all information and agrees with pkaANI cutoff set
        The AEV is a vector with length 1008 representing the environment for the lysine

        Parameters
        ----------
        path : string
            The path of the pdb file that SASA is being calculated for.

        Returns
        -------
        df_aevs : dataframe
            Dataframe with information on chain, residue number and AEV output. Outline:
            Chain   Resid   aev
            x       x       [x]

        Example
        -------
        >> print(calculate_aevs(1ubq.pdb))
        Chain Resid                                                aev
        0     A     6  [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ...
        1     A    11  [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ...
        2     A    27  [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ...
        3     A    29  [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ...
        4     A    33  [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ...
        5     A    48  [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ...
        6     A    63  [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, ...
        '''

        # define output dataframe
        df_aevs = pd.DataFrame(columns=["Chain", "Resid", "aev"])

        # 1: prepare the biobox structure, take the species and coordinates and convert to Atoms structure, find the locations of the NZ atoms within the lysines in the strucure
        try:
            M = bb.Molecule(path)
            coords_nz, idx_nz = M.atomselect("*", "LYS", "NZ", use_resname=True, get_index=True)
            all_coords, idx = M.atomselect('*','*','*', get_index=True)
            list_resids = list(M.data['resid'][idx_nz])
            list_chains = list(M.data['chain'][idx_nz])
        except Exception as e:
            print(f'AEV Calculations: 1 - could not create atomic structure representation: {e}')
            return

        # 2: Iterate over the protein structure to cut out substructures and calculate an AEV at each of these.
        try:
            for j, lys_coord in enumerate(coords_nz):
                # 2.1: for the NZ atom of the lysine, find all the atoms within the cutoff distance and create a substructure
                list_close_points = []
                distance_cut_off = 6  # current cutoff for substructure from analysis done on different cutoffs and matching pkaANI
                for i, coord in enumerate(all_coords):
                    try:
                        x_dist = (lys_coord[0] - coord[0])**2
                        y_dist = (lys_coord[1] - coord[1])**2
                        z_dist = (lys_coord[2] - coord[2])**2
                        distance = np.sqrt(x_dist + y_dist + z_dist)
                        if distance < distance_cut_off:
                            list_close_points.append(idx[i])
                    except Exception:
                        continue

                S = M.get_subset(idxs=list_close_points)
                chain = list_chains[j]
                resid = list_resids[j]
                temp_atom_species = S.data['atomtype']
                temp_coords = S.coordinates[0]  # take the coords from the molecule read in through biobox
                temp_structure = Atoms(temp_atom_species, temp_coords)
                temp_idx_nz = S.atomselect("*", "LYS", "NZ", use_resname=True, get_index=True)[1]
                aevs = None

                # 2.2: calculate the AEV for the subset of the protein and add this to the output dataframe
                try:
                    device_specs = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
                    ANI = torchani.models.ANI2x(periodic_table_index=True).to(device=device_specs)
                    species = ANI.species_to_tensor(temp_structure.get_chemical_symbols()).unsqueeze(0)
                    ani_coords = torch.tensor(temp_structure.get_positions(), dtype=torch.float32).unsqueeze(0)
                    species = species.to(device_specs)
                    ani_coords = ani_coords.to(device_specs)
                    aevs = ANI.aev_computer((species, ani_coords)).aevs
                    lys_nz_location = list_close_points.index(idx_nz[j])
                    aevs = aevs[0,lys_nz_location,:]
                    aevs = aevs.tolist()
                except Exception as e:
                    print(f'AEV Calculations: could not create AEV for resid {idx_nz[j]} of protein {path}, error: {e}')

                # 2.3: Append the new AEV to the output dataframe
                aev_to_append = {'Chain': chain, 'Resid': resid, 'aev': aevs}
                df_aevs = pd.concat([df_aevs, pd.DataFrame([aev_to_append])], ignore_index=True)

        except torch.cuda.OutOfMemoryError:
            # potential that calculating the AEVs could overload the gpu, if too much memory, catch this and skip the file
            print(f'CUDA memory error with file: {path}, skipping')
            return df_aevs
        except MemoryError:
            # potential that calculating the AEVs could overload the cpu, if too much memory, catch this and skip the file
            print(f'CPU memory error with file: {path}, skipping')
            return df_aevs
        except Exception as e:
            print(f'AEV Calculations: 2 - could not create the AEVs for the protein for protein {path}, error: {e}')
            return df_aevs

        # 4: if everything has worked, return the dataframe with the AEVs for the protein
        #print(df_aevs)
        return df_aevs


    def calculate_das(self, path):
        '''
        Calculate the Dynamically Accessible Surface (DAS) of the NZ atom in the lysine structure
        This is effectively the number of positions that the NZ atom can take within the structure of the protein

        Method
        ------
        Uses biobox functionality to calculate the value
        Create a molecule for the protein structure from the bb.Molecule class
        Use the bb.Xlink class to setup the linking module
        Use the hidden method .__get_half_sphere() to work out the das value
        As the density of points in the sphere of the NZ atom of the lysine is constant,
            the das value is the number of points that are accessible


        Parameters
        ----------
        path : string
            The path of the pdb file that DAS is being calculated for.

        Returns
        -------
        df_das : dataframe
            Dataframe with information on chain, residue number and DAS output. Outline:
            Chain   Resid   das
            x       x       [x]

        Example
        -------
        >> print(calculate_das(1ubq.pdb))
        Chain  Resid  das
        0     A      6   36
        1     A     11   47
        2     A     27   22
        3     A     29   37
        4     A     33   46
        5     A     48   34
        6     A     63   30
        '''

        # 1: Load in the structure and locate all the NZ atoms within the lysines, calculate the list of chains and list of resids to go with this
        try:
            M = bb.Molecule(path)
            idx_nz = M.atomselect('*', 'LYS', 'NZ', use_resname=True, get_index=True)[1]
            lys_res_nums = list(M.data['resid'][idx_nz])
            list_chains = list(M.data['chain'][idx_nz])
        except Exception as e:
            print(f'DAS Calculation: 1 - could not load and identify the NZ atoms within the lysines of the structure: {e}')

        # 2: Setup the Xlink module and create the half spheres
        try:
            XL = bb.Xlink(M)
            das_output = []
            for lys_nz_idx in idx_nz:
                # the parameteres (pts_surf, thresh, radii) for the _get_half_sphere are already set for lysine residues
                # therefore the only parameter that needs to be set is i: this is the index of the atom of interest within the lysine
                half_sphere_coords = XL._get_half_sphere(i=lys_nz_idx)
                # as the density of points created by the get half sphere is constant for any setup,
                # therefore can just count the number of coordinates that are returned for a measure for SASA Path
                das_output.append(len(half_sphere_coords))
        except Exception as e:
            print(f'DAS Calculation: 2 - Failed to calculate the half spheres for the NZ atoms within the lysines: {e}')

        # 3: Create dataframe to return
        df_das = pd.DataFrame(columns=["Chain", "Resid", "das"])
        try:
            df_das['Chain'] = list_chains
            df_das['Resid'] = lys_res_nums
            df_das['das'] = das_output
        except Exception as e:
            print(f'DAS Calcualtion: 3 - Failed to create datafame to append to the overall dataframe: {e}')

        return df_das



if __name__ == "__main__":


    f1 = "Demo{os.sep}curated{os.sep}1M2E-alt-1.pdb"

    from uniprot import Uniprot
    from protein import PDB

    print("Scanning UNIPROT...")
    UP = Uniprot()
    UP.get_protein_data("P0CG47")

    df = UP.df.iloc[7:9]

    print("Gathering proteins")
    PDB = PDB()
    PDB.gather_proteins(df)

    print("Measuring...")
    M = Measure(PDB.df)
    M.measure_dataframe()
