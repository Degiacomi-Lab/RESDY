import re
import os
import io
import logging
import datetime
import shutil
import subprocess
import glob
import time
from datetime import date
from multiprocessing import cpu_count
from multiprocessing import Manager
from multiprocessing.pool import Pool
from contextlib import redirect_stdout
import pandas as pd
import numpy as np
import biobox as bb
import matplotlib.pyplot as plt
#import dill


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

# Frustration packages
try:
    import frustratometer
except Exception as e:
    print(f"frustratometer unavailable. Unable to calculate frustration. Error: {e}")

# Melodia packages
try:
    import melodia_py as mel
except Exception as e:
    print(f"melodia unavailable. Unable to calculate melodia. Error: {e}")


class Measure(object):

    def __init__(self, df_input, outdir="result", activate_log=False, log_path='measure_log.txt',
                 features=['propka', 'pkaANI', 'sasa', 'depth', 'aev', 'das', 'seqcharge'],
                 parallel=False, include_modified=False, report_errors= True):
        '''
        Initialisation of the Measure class. This class provides all the resources to measure
        specific quantities for the protein structures given as input

        Parameters
        ----------
        df_input -> dataframe
            The input dataframe containing information on the structures over which the measurements
            will be done. This is usually the output given from the curation steps (pdb.df).
            This contains the filepath, Uniprot_Entry, PDB_Code, etc
            TODO: finish this description with full set of required column names here.
        outdir -> string
            The name of the directory where the measurement output will be written to.
        activate_log -> bool
            By default a log is produced for the measurements, the option here enables a more
            detailed log of the measurements work for debugging. TODO Check this
        log_path -> string
            The name of the output file which contains the log of the measurements.
            This file can be used to create the measurement csv file through using the
            recover_from_log() function.
        features -> list
            The list of measurements that you wish to use on the given structures. Select which
            of the following options to use: 'propka', 'pkaANI', 'sasa', 'depth', 'aev',
            'das', 'seqcharge', 'melodia', 'frustration', 'density', 'legolas', 'writhing',
            'curvature', 'torsion', 'arc_length', 'phi', 'psi'
        parallel -> bool
            Option to run the measurements in parallel.
        include_modified -> bool
            Option to include lysines that have been seen to be modified in the measurements
            analysis. If False, only lysines of type 'LYS' will be included in the measurements.
            If True, lysines of types 'LYN' will be included in the measurements as well as all
            'LYS' residues. In either case, a column will be included stating if the measured
            residue is a modified one.
        report_errors -> bool
            Option to record any of the protein files which are giving errors when measures
            calculations are being performed. This will write the file and the error to a separate
            text document labelled "measures_errors_{date}.txt".

        '''
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

        self.outdir = outdir
        self.df_input = df_input
        self.folder = os.path.join(outdir, "curated")

        self.features = features
        self._setup_measures(features)
        pd.set_option("display.max_columns", None)
        pd.reset_option('display.max_rows')

        # for restarting
        self.current_index = 0
        self.progress_index = 0
        self.pdb_only_files_to_ignore = []

        # document failed pdb files
        self.wrong_pdb_file = []
        self.report_errors = report_errors
        if report_errors:
            new_file_name = f'meaures_errors_{date.today()}.txt'
            while os.path.exists(new_file_name):
                if '_no' in new_file_name:
                    error_file_num = int(new_file_name.split('_no')[-1].split('.')[0])
                    new_file_name = f'measure_errors_{date.today()}_no{(error_file_num + 1)}.txt'
                else:
                    new_file_name = f'measure_errors_{date.today()}_no{1}.txt'
            with open(new_file_name, 'w') as error_f1:
                error_f1.write(f'New measures errors file created at {datetime.datetime.now()}\n')
        self.error_filename = new_file_name
        print(self.error_filename)

        # modified lysine management
        self.include_modified = include_modified

        # for parallel measurements
        self.parallel = parallel
        self.files_to_analyse = []
        self.parallel_items = {}

        self.legolas_loc = '/home/gweston/Documents/extra_packages/legolas-main/test'

        # Check that all files in DataFrame appear at least once in folder
        # find all AlphaFold entries
        files_af=[os.path.basename(c).split(".")[0] for c in glob.glob(os.path.join(self.folder, "*pdb"))]
        # find all PDB entries
        files_pdb=[os.path.basename(c).split("-")[0] for c in glob.glob(os.path.join(self.folder, "*pdb"))]
        for f in df_input["PDB_Code"].values:
            if f not in files_af and f not in files_pdb:
                print(f'WARNING: {f} not found in folder {self.folder}')

        self.pkaoutdir = os.path.join(outdir, "propkaoutput")
        if not os.path.exists(self.pkaoutdir):
            os.makedirs(self.pkaoutdir)

        self.PDB_only = False

        if 'Uniprot_Entry' in self.df_input.columns:
            columns = ['Uniprot_Entry', 'PDB_Code', 'Method', 'Resolution', 'Chain', 'Resid']
            self.df = pd.DataFrame(columns=columns)
        else:
            self.PDB_only = True
            columns = ['PDB_Code', 'Chain', 'Resid']
            self.df = pd.DataFrame(columns = columns)

    def _setup_measures(self, features):
        '''
        convert a list of features into a measuring protocol
        '''

        # measures to carry out [label for DataFrame column, and function evaluating a file]
        # functions must return a dataframe [chain, resid, measure]
        self.measures = []
        melodia_added = False
        frustration_added = False
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
            elif m == 'seqcharge':
                self.measures.append([m, self.calculate_seqcharge])
            elif m == 'legolas':
                self.measures.append([m, self.calculate_legolas])
                self.legolas_output_path = os.path.join(self.outdir, 'legolas')
                if not os.path.exists(self.legolas_output_path):
                    os.mkdir(self.legolas_output_path)
            elif m in ['frustration', 'density']:
                if not frustration_added:
                    self.measures.append(['frustration', self.calculate_frustration])
                    frustration_added = True
            elif m == 'melodia':
                self.measures.append([m, self.calculate_melodia])
                melodia_added = True
                self.features += ['curvature', 'writhing', 'torsion', 'arc_length', 'phi', 'psi']
                self.features.remove('melodia')
            elif m in ['curvature', 'writhing', 'torsion', 'arc_length', 'phi', 'psi']:
                if not melodia_added:
                    self.measures.append(['melodia', self.calculate_melodia])
                    melodia_added = True
            else:
                raise Exception(f"measure {m} unknown")


    def measure_data(self):
        '''
        Determine the appropriate measures function to call based on the combination of
        running PDB_only and in parallel, reducing the number individual functions that
        the user will have to call themselves.

        Example
        -------
        M.measure_data()
        '''
        match (self.PDB_only, self.parallel):
            case (False, False):
                # not PDB only and not parallel
                self.measure_dataframe()
            case (False, True):
                # not PDB only and parallel:
                self.measure_dataframe_parallel()
            case (True, False):
                # PDB only and not parallel
                self.measure_PDB_only()
            case (True, True):
                # PDB only and parallel
                print('This setup does not currently have a method, please change the setup')
        
        self._cleanup_calculation_files()


    def _report_error_to_file(self, measurement_stage, path, error):
        '''
        Helper function to remove redundant code writing errors in the measurements to the
        measurement error log file.
        
        Parameters
        ----------
        measurement_stage -> string
            The stage of measurements that has caused the error with the file, eg propka 1
        path -> string
            The path of the pdb file that the measurement has been attempted on
        error -> string
            The error that has been produced at that step of the measurement when it has been
            attempted to extract features from the pdb file
        
        Example
        -------
        self._report_error_to_file('propka 1', path, e)
        '''
        with open(self.error_filename, 'a') as e_f:
            e_f.writelines('--------------------------------------------------------------------------\n')
            e_f.writelines(f'{measurement_stage} calc error\n')
            e_f.writelines(path + '\n')
            e_f.writelines(error + '\n')


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
        # sort by uniprot code to give order to output after parallel run
        if not self.PDB_only:
            self.df = self.df.sort_values(by='Uniprot_Entry')
        self.df.to_csv(os.path.join(self.outdir, outname), index_label=False, index=False)


    def measure_dataframe_parallel(self):
        '''
        TODO FINISH THIS
        Function to measure specified features for all the structure files curated earlier in the programme.
        Will take a list of the required proteins, finds associated curated structures and runs the reequired measurement functions.
        Results are saved to memory and a log file produced at the same time. (M.save_state() can be used to save the data to a csv)

        Method
        ------
        Create list of files that have been curated into the self.outdir directory.
        Iterate over the list of the files, check if structure file is

        Parameters
        ----------

        Example
        -------
        >> M.measure_dataframe()
        '''
        # use a different method if handling pdb codes only
        if self.PDB_only:
            return 'Call PDB_only method instead'

        files = glob.glob(os.path.join(self.folder, "*pdb"))
        # remove the files which have pkaani in the name as these are output files from pkaani
        files = [file for file in files if 'pkaani' not in file]
        self.files_to_analyse = files

        first_index = self.df_input.index[0]
        total_structures = len(files)
        current_structure = 0
        overall_st = time.time()
        print(f'Total number of structures to analyse: {total_structures}')

        # 2 options: either linear or parallel measurements
        match self.parallel:
            case True:
                # determine how to parallelise
                n_cores_to_use = cpu_count() - 2
                print('>> Measurements running in parallel')
            case False:
                # use a singular core for step by step processing
                n_cores_to_use = 1
                print('>> Measurements running step by step')

        with Manager() as manager:
            # create lock to avoid multiple parts writing to output files at the same time
            lock = manager.Lock()
            ns_measures = manager.Namespace()
            ns_measures.df = self.df
            # prepare the inputs for the parallelisation
            items = []
            for i, r in self.df_input.iterrows():
                pdb_code = r["PDB_Code"]
                chains = r["Chains"].split("/")
                method = r['Method']
                res = r['Resolution']
                uniprot_code = r["Uniprot_Entry"]
                file_details = [uniprot_code, pdb_code, method, res, chains]
                items.append([file_details, lock, ns_measures])
            with Pool(n_cores_to_use, maxtasksperchild=10) as pool:
                result = pool.starmap_async(self._measure_file, items)
                result.wait()
                print(result)
                self.df = ns_measures.df

        # remove possible duplicated rows (if restarted)
        try:
            self.df = self.df.drop_duplicates(subset=None, keep='first', inplace=False, ignore_index=True)
        except Exception as e_one:
            try:
                self.df = self.df.loc[self.df.astype(str).drop_duplicates(subset=None, keep='first', inplace=False, ignore_index=True)]
            except Exception as e_two:
                print(f'Failed to remove duplicates from measurement dataframe: {e_two}')
                pass


    def measure_dataframe(self):
        '''
        TODO FINISH THIS
        Function to measure specified features for all the structure files curated earlier in the programme.
        Will take a list of the required proteins, finds associated curated structures and runs the reequired measurement functions.
        Results are saved to memory and a log file produced at the same time. (M.save_state() can be used to save the data to a csv)

        Method
        ------
        Create list of files that have been curated into the self.outdir directory.
        Iterate over the list of the files, check if structure file is

        Parameters
        ----------

        Example
        -------
        >> M.measure_dataframe()
        '''
        # use a different method if handling pdb codes only
        if self.PDB_only:
            return 'Call PDB_only method instead'

        if self.parallel:
            return 'Call measure_dataframe_parallel instead'

        files = glob.glob(os.path.join(self.folder, "*pdb"))
        # remove the files which have pkaani in the name as these are output files from pkaani
        files = [file for file in files if 'pkaani' not in file]
        self.files_to_analyse = files

        first_index = self.df_input.index[0]
        total_structures = len(files)
        current_structure = 0
        overall_st = time.time()
        print(f'Total number of structures to analyse: {total_structures}')

        for i, r in self.df_input.iterrows():
            pdb_code = r["PDB_Code"]
            chains = r["Chains"].split("/")
            method = r['Method']
            res = r['Resolution']
            uniprot_code = r["Uniprot_Entry"]
            self.current_index = i
            files_list = self.files_to_analyse

            # calculate features values from all PDB files associated with specific DataFrame entry
            for f in files_list:
                # check the file for the required pbd code, if not there, skip
                if pdb_code not in f:
                    continue

                print(f'Analysing structure {current_structure}/{total_structures}')
                tstart = time.time()
                print(f"\n> File: {f}")

                # create temporary DataFrame for data of current file,
                # to be then appended to main DataFrame self.df

                columns = ['Uniprot_Entry', 'PDB_Code', 'Method', 'Resolution', 'Chain', 'Resid']
                df_currentfile = pd.DataFrame(columns=columns)

                # append to temporary DataFrame all lysines in the file of interest
                try:
                    M = bb.Molecule(f) # sometimes bb does not work with a pdb file
                except Exception as e:
                    self.wrong_pdb_file.append(f)
                    print(f'Failed to produce bb for pdb file with error: {e}')
                    continue

                _, idxs = M.atomselect("*", ["LYS"], ["CA"], get_index=True, use_resname=True)
                for i in idxs:

                    #save only lysine entries from chain of interest
                    if M.data["chain"].values[i] not in chains:
                        continue

                    data = ({'Uniprot_Entry': uniprot_code,
                        'PDB_Code': f.split(".")[0],
                        'Method': method,
                        'Resolution': res,
                        'Chain': M.data["chain"].values[i],
                        'Resid': M.data["resid"].values[i]})

                    df_currentfile = pd.concat([df_currentfile, pd.DataFrame.from_records(data, index=[0])], ignore_index=True)

                print(f">> {len(df_currentfile)} lysines of interest found")

                # iterate over measures to carry out (according to self.measures)
                for meas in self.measures:
                    print(f">> evaluating {meas[0]}...")
                    try:
                        #df_currentfile[meas[0]] = np.nan # create new column for measure
                        result = meas[1](f)
                        df_currentfile = self._combine_dataframes(df_currentfile, result, meas[0])

                    except Exception as e:
                        print(f"ERROR: {e}")
                        continue

                processing_time = round((time.time()-tstart), 2)
                print(f">> file processed in {processing_time} seconds.")
                #average_time_per_file = round(((time.time()- overall_st) / current_structure), 2)
                #print(f'>> Time average per file: {average_time_per_file} seconds.')
                #sec_remaining = average_time_per_file * (total_structures + 1 - current_structure)
                #time_remaining_str = str(datetime.timedelta(seconds=sec_remaining))
                #print(f'Predicted time remaining: {time_remaining_str}')

                # document the data to a log file
                if self.activate_log:
                    if df_currentfile.empty is False:
                        pd.set_option('display.max_colwidth', None,
                                    'display.width', None,
                                    'max_seq_items', None,
                                    "display.max_rows", None)
                        try:
                            self.logger.info(df_currentfile)
                            self.logger.info('--------------------------------------------------------------------------')
                        except Exception as e:
                            print(f'Error in logging: {e}')

                        # reset the pandas display options back to default for regular displaying
                        pd.reset_option('display.max_colwidth')
                        pd.reset_option('display.width')
                        pd.reset_option('max_seq_items')
                        pd.reset_option('display.max_rows')

                #append temporary DataFrame with all measures on a single file to main DataFrame
                if not df_currentfile.empty:
                    self.df = pd.concat([self.df, df_currentfile], ignore_index=True)

        # remove possible duplicated rows (if restarted)
        try:
            self.df = self.df.drop_duplicates(subset=None, keep='first', inplace=False, ignore_index=True)
        except Exception as e_one:
            try:
                print(f'Failed to remove duplicates from measurement dataframe, trying new method: {e_one}')
                self.df = self.df.loc[self.df.astype(str).drop_duplicates(subset=None, keep='first', inplace=False, ignore_index=True)]
            except Exception as e_two:
                print(f'Failed to remove duplicates from measurement dataframe: {e_two}')
                pass



    def _measure_file(self, file_details, lock, ns):
        '''
        Take a file and calculate the required measurements for this.
        Return the dataframe of the calculations to the overall self.df

        Method
        ------
        Take the given information about the file and through the given file_list,
        find the files to analyse. Through the structure, identify all lysines and
        create a temporary dataframe for the results. Iterate over the required measurements
        (self.measures) and insert results into the temporary dataframe. Append temporary
        dataframe to main dataframe.

        Parameters
        ----------
        file_details : list
            list of details for the file that has been selected to be calculated
            takes the form of [uniprot_code, pdb_code, method, res, chains]
        
        lock : multiprocessing manager lock
            lock used to stop processes writing to output files and dataframes at the same time

        ns : multiprocessing manager namespace
            allows appending the dataframe results to a shared out dataframe,
            this is then transferred to self.df

        Example
        -------
        self._measure_file(file_details, files_list)
        '''

        # decompose the file_details into variables
        uniprot_code, pdb_code, method, res, chains = file_details
        files_list = self.files_to_analyse

        # calculate features values from all PDB files associated with specific DataFrame entry
        for f in files_list:
            # check the file for the required pbd code, if not there, skip
            if pdb_code not in f:
                continue

            terminal_out_statements = []
            #print(f'Analysing structure {current_structure}/{total_structures}')
            tstart = time.time()
            #print(f"\n> File: {f}")
            terminal_out_statements.append(f"\n> File: {f}")

            # create temporary DataFrame for data of current file,
            # to be then appended to main DataFrame self.df

            columns = ['Uniprot_Entry', 'PDB_Code', 'Method', 'Resolution', 'Chain', 'Resid']
            df_currentfile = pd.DataFrame(columns=columns)

            # append to temporary DataFrame all lysines in the file of interest
            try:
                M = bb.Molecule(f) # sometimes bb does not work with a pdb file
            except Exception as e:
                self.wrong_pdb_file.append(f)
                terminal_out_statements.append(f'Failed to produce bb for pdb file with error: {e}')
                continue

            _, idxs = M.atomselect("*", ["LYS"], ["CA"], get_index=True, use_resname=True)
            for i in idxs:

                #save only lysine entries from chain of interest
                if M.data["chain"].values[i] not in chains:
                    continue


                data = ({'Uniprot_Entry': uniprot_code,
                    'PDB_Code': f.split(".")[0],
                    'Method': method,
                    'Resolution': res,
                    'Chain': M.data["chain"].values[i],
                    'Resid': M.data["resid"].values[i]})

                df_currentfile = pd.concat([df_currentfile, pd.DataFrame.from_records(data, index=[0])], ignore_index=True)

            #print(f">> {len(df_currentfile)} lysines of interest found")
            terminal_out_statements.append(f">> {len(df_currentfile)} lysines of interest found")

            # iterate over measures to carry out (according to self.measures)
            for meas in self.measures:
                #print(f">> evaluating {meas[0]}...")
                terminal_out_statements.append(f">> evaluating {meas[0]}...")
                try:
                    df_currentfile[meas[0]] = np.nan # create new column for measure
                    out_print_trap = io.StringIO()
                    with redirect_stdout(out_print_trap):
                        result = meas[1](f) # run measurement
                    terminal_out_statements.append(out_print_trap.getvalue())
                    df_currentfile = self._combine_dataframes(df_currentfile, result, meas[0]) #insert measures into temporary DataFrame

                except Exception as e:
                    #print(f"ERROR: {e}")
                    terminal_out_statements.append(f"ERROR: {e}")
                    continue

            processing_time = round((time.time()-tstart), 2)
            #print(f">> file processed in {processing_time} seconds.")
            terminal_out_statements.append(f">> file processed in {processing_time} seconds.")
            #average_time_per_file = round(((time.time()- overall_st) / current_structure), 2)
            #print(f'>> Time average per file: {average_time_per_file} seconds.')
            #sec_remaining = average_time_per_file * (total_structures + 1 - current_structure)
            #time_remaining_str = str(datetime.timedelta(seconds=sec_remaining))
            #print(f'Predicted time remaining: {time_remaining_str}')

            with lock:
                # write all terminal outputs for file
                for statement in terminal_out_statements:
                    print(statement)

                # document the data to a log file
                if self.activate_log:
                    if df_currentfile.empty is False:
                        pd.set_option('display.max_colwidth', None,
                                      'display.width', None,
                                      'max_seq_items', None,
                                      "display.max_rows", None)
                        try:
                            self.logger.info(df_currentfile)
                            self.logger.info('--------------------------------------------------------------------------')
                        except Exception as e:
                            print(f'Error in logging: {e}')

                        # reset the pandas display options back to default for regular displaying
                        pd.reset_option('display.max_colwidth')
                        pd.reset_option('display.width')
                        pd.reset_option('max_seq_items')
                        pd.reset_option('display.max_rows')

                #append temporary DataFrame with all measures on a single file to main DataFrame
                if not df_currentfile.empty:
                    ns.df = pd.concat([ns.df, df_currentfile], ignore_index=True)


    def recover_from_log(self, log_path):
        '''
        Take the log file produced through running measure_dataframe() and convert this to a csv

        Method
        ------
        Read in the log file (measure_log.txt)
        Work out the columns from the header
        If the headers can't be worked out, ask for input to match up columns
        Read in data
        Sets self.df to be the data output recovered from the log file.

        Parameters
        ----------
        log_path : string
            the file name for the log file to convert

        Returns
        -------
        log_to_df : dataframe
            Dataframe containing all the measurements that were in the given log file

        Example
        -------
        M.recover_from_log()
        '''
        if self.PDB_only:
            return 'Function not callable.'

        base_columns = ['Uniprot_Entry', 'PDB_Code', 'Method', 'Resolution', 'Chain', 'Resid']
        log_to_df = pd.DataFrame(columns=base_columns)
        log_path = os.path.join(self.outdir, log_path)
        print(f'Recovering measured data from file: {log_path}')
        print('WARNING: could take up to a few minutes depending on the number of measurements completed.')

        test_lines = 0
        columns_all_set = False
        potential_col_names = {'1': 'propka', '2': 'pkaANI', '3': 'sasa',
                               '4': 'depth', '5': 'aev', '6': 'das'}
        with open(log_path, "rb") as f:
            num_lines = sum(1 for _ in f)
        curr_line = 0

        with open(log_path) as inf:
            for line in inf:
                curr_line += 1
                # check if it is a header line, check if doesn't start with number or -
                if line[0].isalpha() or line[0] == ' ':
                    # found a header line
                    parts = line.split()
                    # check that columns have been written to the log file correctly
                    if columns_all_set:
                        continue
                    elif len(parts) <= 6 and not columns_all_set:
                        print('Columns were not set correctly in the log file.')
                        print(f'The first 6 columns are assumed to be: {base_columns}')
                        continue
                    elif len(parts) >= 6 and not columns_all_set:
                        # if all seems correct with the writing
                        # check that all the columns can be found in the current columns, if not, add in
                        for part in parts:
                            if part not in base_columns:
                                base_columns.append(part)
                        continue

                line_splitter_bool = all(a == '-' for a in line.strip())
                if line_splitter_bool:
                    line = ''

                # split information into parts keeping the AEV as one unit
                parts = re.split(r'([\w.,\/-]+)|(\[.+?\])', line)
                if len(parts) == 0:
                    continue
                # remove None elemnts from matching and reduce all space values to ''
                parts = [elmnt.strip() for elmnt in parts if elmnt is not None]
                # remove '' elements from the list
                parts = [elmnt for elmnt in parts if elmnt != '']
                parts = parts[1:]

                if len(parts) == 0:
                    continue
                num_parts = len(parts)
                if num_parts != len(base_columns):
                    while num_parts != len(base_columns):
                        print('Need to set a column header')
                        print(f'Options for columns are: {potential_col_names}')
                        print(f'Please enter the number corresponding to the header required for the column which contains the following value: {parts[len(base_columns)]}')
                        new_header_val = input('Enter the number for the new column header: ')
                        while True:
                            if not new_header_val.isnumeric():
                                new_header_val = input('Enter the number for the new column header: ')
                            elif 1 <= int(new_header_val) <= len(potential_col_names):
                                break
                            else:
                                new_header_val = input('Enter the number for the new column header: ')
                        base_columns.append(potential_col_names[new_header_val])
                    columns_all_set = True
                data = dict(zip(base_columns, parts))
                log_to_df = pd.concat([log_to_df, pd.DataFrame.from_records(data, index=[0])], ignore_index=True)
                test_lines += 1
                print(f'Progress analysing log file: {round((curr_line/num_lines)*100, 2)} %\r', end='', flush=True)

        self.df = log_to_df
        print('Data recovered from log file')
        print(f'Numer of measurements read: {len(log_to_df)}')
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

        if self.PDB_only:
            return 'restart_measure() function not callable when using PDB_only'

        print('Restarting measurements')
        # 1. Analyse the measures log file to create a list of files that were analysed
        log_path = os.path.join(self.outdir, log_path)
        print(f'Finding measured proteins from log file: {log_path}')
        proteins_completed = []
        with open(log_path, "rb") as f:
            num_lines = sum(1 for _ in f)
        curr_line = 0
        with open(file=log_path, mode='r') as lpf:
            for line in lpf:
                curr_line += 1
                if line[0].isalpha() or line[0] == ' ' or line[0] == '-':
                    continue
                parts = line.split()
                protein_code = parts[1]
                proteins_completed.append(protein_code)
                print(f'Progress analysing log file: {round((curr_line/num_lines)*100, 2)} %\r', end='', flush=True)
        print('Measured proteins recovered from log file')

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
        # 3. Restart the measure_dataframe() with the new file list
        print(f'Continuing measurements. {len(self.df_input)} proteins to measure.')
        self.measure_dataframe_parallel()


    def _combine_dataframes(self, target, to_merge, col_name):
        '''
        target is a DataFrame to be filled with data, to_merge contains the data.
        Values to insert are indexed in both array by two columns: Chain and Resid.
        '''
        # e.g. self._combine_dataframes(df, result, meas[0])

        for i, r in target.iterrows():

            chain_value = r["Chain"]
            resid_value = r["Resid"]

            idx = np.where((to_merge["Chain"] == chain_value) & (to_merge["Resid"].astype(int) == resid_value))
            if len(idx[0]) == 0:
                continue

            # account for measurements that have special cases
            # aevs - add the list of aevs in one column to the overall dataframe
            if col_name == 'aev':
                target['aev'] = target['aev'].astype('object')
            # melodia - check over all the required features to add and add these back in to the overall dataframe
            if col_name == 'melodia':
                melodia_features = ['curvature', 'writhing', 'torsion', 'arc_length', 'phi', 'psi']
                for feature in self.features:
                    if feature in melodia_features:
                        target.at[i, feature] = to_merge.loc[idx[0][0], feature]
            elif col_name == 'frustration':
                frust_features = ['frustration', 'density']
                for feature in self.features:
                    if feature in frust_features:
                        target.at[i, feature] = to_merge.loc[idx[0][0], feature]
            else:
                target.at[i, col_name] = to_merge.loc[idx[0][0], col_name]
        
        return target


    def measure_PDB_only(self):
        '''
        TODO FINISH THIS
        Function to measure specified features for a set of pdb files. Takes a list of pdb files,
        finds associated curated structures and runs the required measurement functions.
        Results are saved to memory and a log file produced at the same time.
        M.save_state() can be used to save the data to a csv.

        Method
        ------

        Example
        -------
        >> M.measure_PDB_only()
        '''
        if not self.PDB_only:
            print('Called measure_PDB_only() when running not on PDB_only. Call measure_dataframe() instead or change to run PDB_only.')
            return 'Calling the wrong method.'

        files = glob.glob(os.path.join(self.folder, "*pdb"))
        tstart_overall = time.time()
        num_pdb_files = len(self.df_input) + self.progress_index

        # remove any potential duplicates from the input dataframe
        self.df_input = self.df_input.drop_duplicates()
        if 'completed' not in self.df_input.columns:
            self.df_input['completed'] = False

        for pdb_idx, row in self.df_input.iterrows():
            pdb_code = row['PDB_Code']

            if row['completed']:
                continue

            # calculate features values from all PDB files associated with specific DataFrame entry
            for f in files:

                if pdb_code not in f:
                    continue

                if f in self.pdb_only_files_to_ignore:
                    continue

                tstart = time.time()
                print(f"\n> File: {f}")

                # create temporary DataFrame for data of current file,
                # to be then appended to main DataFrame self.df

                columns = ['PDB_Code', 'Chain', 'Resid']
                df_currentfile = pd.DataFrame(columns=columns)


                # append to temporary DataFrame all lysines in the file of interest
                try:
                    M = bb.Molecule(f) # sometimes bb does not work with a pdb file
                except Exception as e:
                    print(f'Failed to create biobox molecule for file {f} with error: {e}')
                    self.wrong_pdb_file.append(f)
                    continue

                df_idx, idxs = M.atomselect("*", ["LYS"], ["CA"], get_index=True, use_resname=True)
                for i in idxs:

                    #save only lysine entries from chain of interest
                    #if M.data["chain"].values[i] not in chains:
                    #    continue


                    data = ({'PDB_Code': f.split(".")[0],
                        'Chain': M.data["chain"].values[i],
                        'Resid': M.data["resid"].values[i]})

                    df_currentfile = pd.concat([df_currentfile, pd.DataFrame.from_records(data, index=[0])], ignore_index=True)

                print(f">> {len(df_currentfile)} lysines of interest found")

                # iterate over measures to carry out (according to self.measures)
                for meas in self.measures:
                    print(f">> evaluating {meas[0]}...")
                    try:
                        #df_currentfile[meas[0]] = np.nan # create new column for measure  # change GW 11.03.25 - dont need this, new column created anyway, leaving in incase removing creates problems later
                        result = meas[1](f) # run measurement
                        df_currentfile = self._combine_dataframes(df_currentfile, result, meas[0]) #insert measures into temporary DataFrame

                    except Exception as e:
                        print(f"ERROR adding the measurements for file {f} to the dataframe: {e}")
                        continue

                processing_time = round((time.time()-tstart), 2)
                print(f">> file processed in {processing_time} sec.")

                # document the data to a log file
                if self.activate_log:
                    if df_currentfile.empty is False:
                        try:
                            pd.set_option('display.max_colwidth', None,
                                        'display.width', None,
                                        'max_seq_items', None,
                                        "display.max_rows", None)
                            try:
                                self.logger.info(df_currentfile)
                                self.logger.info('--------------------------------------------------------------------------')
                            except Exception as e:
                                print(f'Error in logging: {e}')

                            # reset the display options back to default for regular displaying
                            pd.reset_option('display.max_colwidth')
                            pd.reset_option('display.width')
                            pd.reset_option('max_seq_items')
                            pd.reset_option('display.max_rows')
                        except Exception as e:
                            print(f'Error in logging measurements: {e}')

                #append temporary DataFrame with all measures on a single file to main DataFrame
                if not df_currentfile.empty:
                    self.df = pd.concat([self.df, df_currentfile], ignore_index=True)

                if f.replace('.pdb', '') == pdb_code:
                    self.df_input.at[pdb_idx, 'completed'] = True
                    break

            try:
                avg_time_per_file = (time.time() - tstart_overall) / (pdb_idx + 1)
                time_remaining = datetime.timedelta(seconds=int(round((len(self.df_input) - (pdb_idx + 1)) * avg_time_per_file, 0)))
                perc_prog_measure = round(((pdb_idx + self.progress_index + 1)/num_pdb_files)*100, 2)
                print(f'>> Progress calculating measurements: {perc_prog_measure}%. Predicted time remaining: {time_remaining}s \r', end='', flush=True)
            except Exception as e:
                print(f'>> Broken progress updater: {e} \r', end='', flush=True)



    def restart_measure_pdb_only(self, log_path='measure_log.txt'):
        '''
        A function to restart the measurements calculations for the PDB only function.
        Useful if the initial run of the measurements crashes or gets stuck.
        Works out how far along the simulation was by running an analysis of the measures log file.

        Method
        ------
        Read over the measures log file and collate a list of files that have been analysed.
        Remove the final value from the list as this may not have been done properly.
        Remove completed files from files to do.
        Restart measure_dataframe() with the new list.

        Parameters
        ----------
        log_path : string
            The name of the measures log file
            By default takes the name 'measures_log.txt'

        Example
        -------
        >> M.restart_measure_pdb_only()
        '''

        if not self.PDB_only:
            print('>> restart_measure_pdb_only() function not callable when not running PDB only')
            return 'restart_measure_pdb_only() function not callable when not running PDB only'

        print('>> Preparing to restart measurements on PDB only')
        # 1. Analyse the measures log file to create a list of files that were analysed
        log_path = os.path.join(self.outdir, log_path)
        print(f'>> Finding measured proteins from log file: {log_path}')
        proteins_completed = []
        with open(log_path, "rb") as f:
            num_lines = sum(1 for _ in f)

        if num_lines == 0:
            print('>> No previous measures data is found in the specified log file. Make sure the log file stated is correct or run measure_pdb_only() from start.')
            return 'No previous measures data is found in the specified log file. Make sure the log file stated is correct or run measure_pdb_only() from start.'

        curr_line = 0
        with open(file=log_path, mode='r') as lpf:
            for line in lpf:
                curr_line += 1
                if line[0].isalpha() or line[0] == ' ' or line[0] == '-':
                    continue
                parts = line.split()
                protein_code = parts[1].split('/')[-1]
                proteins_completed.append(protein_code)
                percent_prog = round((curr_line/num_lines)*100, 2)
                print(f'Progress analysing log file: {percent_prog} %\r', end='', flush=True)
        print('>> Measured proteins recovered from log file')

        # 2. remove the last protein from list incase it wasn't completed fully
        final_protein = proteins_completed[-1]
        proteins_completed = [c for c in proteins_completed if c != final_protein]
        proteins_completed = list(set(proteins_completed))
        self.progress_index = len(proteins_completed)
        # TODO GW 14.01.25 - need to remove this measurement from the measures_log.txt file
        #                    eventually currently doesn't matter too much as will just be
        #                    removed with remove duplicates later

        # 3. Update df_input to only have the files which haven't been analysed yet
        if 'completed' not in self.df_input.columns:
            self.df_input['completed'] = False
        old_len_df_input = len(self.df_input)
        idx_to_remove = []
        files = [a.split('/')[-1] for a in glob.glob(os.path.join(self.folder, "*pdb"))]
        for i, r in self.df_input.iterrows():
            matched_pdb_files = [a.replace('.pdb', '') for a in files if r['PDB_Code'] in a]
            for recover_file in proteins_completed:
                # case 1: exact match code and file - for measuring data from simulations mainly
                if r['PDB_Code'] == recover_file:
                    idx_to_remove.append(i)
                    self.df_input.at[i, 'completed'] = True
                    break
                # case 2: PDB files renamed by curation that are not the PDB code alone
                if recover_file in matched_pdb_files:
                    if recover_file not in self.pdb_only_files_to_ignore:
                        full_recover_file = self.folder + os.sep + recover_file + '.pdb'
                        self.pdb_only_files_to_ignore.append(full_recover_file)

            perc_prog_remove = round((i/len(self.df_input))*100, 2)
            print(f'Progress removing measured files: {perc_prog_remove} %\r', end='', flush=True)

        self.df_input = self.df_input.drop(idx_to_remove)
        new_len_df_input = len(self.df_input)
        lines_df_input_removed = old_len_df_input - new_len_df_input
        # TODO GW 16.01.25 - updated verison for this will use column of completed for everything here, change over to this rather than removing it from the df_input
        # 4. Restart the measure_dataframe() with the new file list
        files_left_to_calc = len(self.df_input) - len(self.pdb_only_files_to_ignore)
        print(f'>> {lines_df_input_removed} exact matches in PDB codes removed from the input list that have already been calculated.')
        print(f'>> {len(self.pdb_only_files_to_ignore)} files to ignore in measurements that have already been calculated.')
        print(f'>> Continuing measurements. {files_left_to_calc} proteins to measure.')
        self.measure_PDB_only()


    def recover_from_log_PDB_only(self, log_path):
        '''
        Take the log file produced through running measure_PDB_only() and convert this to a csv

        Method
        ------
        Read in the log file (measure_log.txt) or other given name.
        Work out the columns from the header.
        If the headers can't be worked out, ask for input to match up columns.
        Read in data.
        Sets self.df to be the data output recovered from the log file.

        Parameters
        ----------
        log_path : string
            the file name for the log file to convert

        Returns
        -------
        log_to_df : dataframe
            Dataframe containing all the measurements that were in the given log file

        Example
        -------
        M.recover_from_log_PDB_only()
        '''
        if not self.PDB_only:
            return 'Function not callable.'

        base_columns = ['PDB_Code', 'Chain', 'Resid']
        log_to_df = pd.DataFrame(columns=base_columns)
        log_path = os.path.join(self.outdir, log_path)
        print(f'Recovering measured data from file: {log_path}')
        print('WARNING: could take up to a few minutes depending on the number of measurements completed.')

        test_lines = 0
        columns_all_set = False
        potential_col_names = {'1': 'propka', '2': 'pkaANI', '3': 'sasa',
                               '4': 'depth', '5': 'aev', '6': 'das',
                               '7': 'Other'}
        with open(log_path, "rb") as f:
            num_lines = sum(1 for _ in f)
        curr_line = 0

        with open(log_path) as inf:
            set_header_line = ''
            dataframe_columns = []
            for line in inf:
                curr_line += 1
                # check if it is a header line, check if doesn't start with number or -
                if line[0].isalpha() or line[0] == ' ':
                    # found a header line
                    if line != set_header_line:
                        parts = line.split()
                        # check that columns have been written to the log file correctly
                        if len(parts) <= 3 and not columns_all_set:
                            print('Columns were not set correctly in the log file.')
                            print(f'The first 6 columns are assumed to be: {base_columns}')
                            continue
                        elif len(parts) >= 3 and not columns_all_set:
                            # if all seems correct with the writing
                            # check that all the columns can be found in the current columns, if not, add in
                            for part in parts:
                                if part not in base_columns:
                                    base_columns.append(part)
                                if part not in dataframe_columns:
                                    dataframe_columns.append(part)
                                    if len(log_to_df) != 0:
                                        log_to_df[part] = None
                            continue

                line_splitter_bool = all(a == '-' for a in line.strip())
                if line_splitter_bool:
                    line = ''

                # split information into parts keeping the AEV as one unit
                parts = re.split(r'([\w.,\/-]+)|(\[.+?\])', line)
                if len(parts) == 0:
                    continue
                # remove None elemnts from matching and reduce all space values to ''
                parts = [elmnt.strip() for elmnt in parts if elmnt is not None]
                # remove '' elements from the list
                parts = [elmnt for elmnt in parts if elmnt != '']
                parts = parts[1:]

                if len(parts) == 0:
                    continue
                num_parts = len(parts)

                # Case 1: Setting the columns when the columns have been messed up and aren't the
                #         same as the data in the log file
                if num_parts > len(base_columns):
                    while num_parts != len(base_columns):
                        print('Need to set a column header')
                        if len(potential_col_names) != 0:
                            print(f'Options for columns are: {potential_col_names}')
                            print(f'Please enter the number corresponding to the header required for the column which contains the following value: {parts[len(base_columns)]}')
                            new_header_val = input('Enter the number for the new column header: ')
                            while True:
                                if not new_header_val.isnumeric():
                                    new_header_val = input('Enter the number for the new column header: ')
                                elif 1 <= int(new_header_val) <= len(potential_col_names):
                                    break
                                else:
                                    new_header_val = input('Enter the number for the new column header: ')
                            if new_header_val == '7':
                                new_header_name = input('Other selected, please enter a unique name for the column: ')
                                base_columns.append(new_header_name)
                            else:
                                base_columns.append(potential_col_names[new_header_val])
                                del potential_col_names[new_header_val]
                        else:
                            new_header = input(f'No more suggested columns available, please enter your column name for the column containing this value:  {parts[len(base_columns)]}')
                            base_columns.append(new_header)
                    columns_all_set = True
                data = dict(zip(base_columns, parts))
                for col in dataframe_columns:
                    if col not in data:
                        data[col] = None

                log_to_df = pd.concat([log_to_df, pd.DataFrame.from_records(data, index=[0])], ignore_index=True)
                test_lines += 1
                print(f'Progress analysing log file: {round((curr_line/num_lines)*100, 2)} %\r', end='', flush=True)

        self.df = log_to_df
        print('Data recovered from log file')
        print(f'Numer of measurements read: {len(log_to_df)}')
        return log_to_df


    def _cleanup_calculation_files(self):
        '''
        Function to remove any temporary or result files created through the calculation
        of the measurements within this class. While all are meant to have been moved
        at the time of calculation, occasionally this fails and leaves some behind.
        Note: please add specific subprocesses if need to add extra cleanup items into
        this function.

        Method
        ------
        Call subprocess calls to move specific sets of files to a specific directory.

        Example
        -------
        >> M._cleanup_calculation_files()
        '''
        dir_files = [f for f in os.listdir() if os.path.isfile(os.path.join(os.getcwd(),f))]
        nmr_cs_file = False
        nmr_parquet_file = False
        propka_pka_file = False
        propka_error_file = False
        for f in dir_files:
            if f.endswith('_cs.csv'):
                nmr_cs_file = True
            elif f.endswith('_cs.parquet'):
                nmr_parquet_file = True
            elif f.endswith('.pka'):
                propka_pka_file = True
            elif f.endswith('_propka_errors.txt'):
                propka_error_file = True
        if nmr_cs_file:
            subprocess.run(f'mv *_cs.csv {self.outdir}{os.sep}legolas{os.sep}', shell=True, check=False)
        if nmr_parquet_file:
            subprocess.run(f'mv *_cs.parquet {self.outdir}{os.sep}legolas{os.sep}', shell=True, check=False)
        if propka_pka_file:
            subprocess.run(f'mv *.pka {self.outdir}{os.sep}propkaoutput{os.sep}', shell=True, check=False)
        if propka_error_file:
            subprocess.run(f'mv *_propka_errors.txt {self.outdir}{os.sep}propkaoutput{os.sep}', shell=True, check=False)
        if self.report_errors:
            with open(self.error_filename, 'r') as f:
                for count, line in enumerate(f):
                    pass
            if count <= 2:
                os.remove(self.error_filename)



    def calculate_pka_propka(self, path):
        '''
        Call PROPKA to calculate the pKa of a file, parse the .pka file to extract lysine data
        parse errors, and return a dataframe containing all measurements not yielding an error.
        
        Method
        ------
        Check if propka has been run before on this protein, otherwise run PROPKA3 on the given
        pdb file. Use the function _parse_propka_errors() to identify any lysines within the structure
        that did not caclulate correctly before searching the output file, extracting the pka
        values produced and writing them to df_propka to output.

        Parameters
        ----------
        path : string
            The path of the pdb file that pKa is being calculated for with PROPKA3.

        Returns
        -------
        df_propka : dataframe
            Dataframe with information on chain, residue number and propka output. Outline:
            Chain   Resid   propka
            x       x       x

        Example
        -------
        >> print(calculate_propka(1ubq.pdb))
        Chain  Resid  propka
        0     A      6   x
        1     A     11   x
        2     A     27   x
        3     A     29   x
        4     A     33   x
        5     A     48   x
        6     A     63   x
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

                if self.report_errors: self._report_error_to_file('PROPKA 1', path, str(e))
                raise Exception(f'Failed to obtain pKa data (PROPKA): {e}') from e

        try:
            propka_lys_fails = self._parse_propka_errors(error_file_name)
        except Exception as e:
            if self.report_errors: self._report_error_to_file('PROPKA 2', path, str(e))
            raise Exception(f"Failed extracting PROPKA errors from output file. {e}")

        try:
            pkafile = code_for_df + '.pka'
            propres = open(pkafile)
        except Exception as e:
            if self.report_errors: self._report_error_to_file('PROPKA 3', path, str(e))
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
                                print(f'Error associated with lysine {int(line[0])} on chain {line[1]} found in the log file at index {idx[0]} for pKa calculation. Value not taken through to measurements.')
                                continue

                        lys_number.append(int(line[0]))
                        chain.append(line[1])
                        pkas.append(float(line[2]))

                    except Exception as e:
                        if self.report_errors: self._report_error_to_file('PROPKA 4', path, str(e))
                        print(f"> Error {e}")
                        continue

            propres.close()
            shutil.move(pkafile, os.path.join(self.pkaoutdir, pkafile))

        except Exception as e:
            propres.close()
            shutil.move(pkafile, os.path.join(self.pkaoutdir, pkafile))
            if self.report_errors: self._report_error_to_file('PROPKA 5', path, str(e))
            raise Exception(f'Failure parsing {code_for_df}.pka, error: {e}') from e

        try:
            df_propka = pd.DataFrame({'Resid':lys_number,
                                      'Chain': chain,
                                      'propka':pkas})

            df_propka.sort_values(by=['propka'], inplace=True)
            df_propka = df_propka.dropna()
            df_propka = df_propka.drop_duplicates(subset=None, keep='first', inplace=False, ignore_index=True)

        except Exception as e:
            if self.report_errors: self._report_error_to_file('PROPKA 6', path, str(e))
            raise Exception(f'Failed to construct pKa (PROPKA) dataframe. {e}') from e

        return df_propka


    def _parse_propka_errors(self, path):
        '''
        parse the PROPKA output file and appends unique chain and resid of any lysines mentioned a DataFrame.
        This list is returned to main and later the residues in it are removed from the df.
        
        Method
        ------
        Go over the propka errors output file that is produced when running. Identify the lines which
        contain information about the lysines within the protein structure analysed that have errors
        associated with them. Extract the chain and resid number from this and append to a dataframe
        to return which contains a set of data on the lysines to remove from the read output pka values.

        Parameters
        ----------
        path : string
            The path of the errors output file from the PROPKA analysis of the pdb file of interest.

        Returns
        -------
        dataframe
            Dataframe with information on chain, residue number for lysines with calculation errors.
            Outline:
            Chain   Resid
            x       x
        '''

        f = open(path, 'r')
        list_remove = list()
        cnt = 0
        for line in f:
            cnt += 1
            lys_raw = re.findall(r'LYS [\d]*[\s][\w]*', line)

            for line in lys_raw:
                words = line.split(' ')
                resid = words[1]
                chain = words[2]
                if resid != "" and chain != "":
                    list_remove.append([chain, resid])

            lys_raw_2 = re.findall(r'[\d]*-LYS \(\w\)', line)

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
        '''
        A second method for calculating the pKa of the NZ atom of the lysines within the protein
        structure, this time using the external program, pKaANI.

        Method
        ------
        Call a subprocess to open the pkaani program with the desired pdb file given through path.
        Find the log file produced from running this and search this to find the values produced
        for LYS. Translate the data found within the log file to the dataframe to be returned.

        Parameters
        ----------
        path : string
            The path of the pdb file that pKa is being calculated for with pKaANI.

        Returns
        -------
        df_pkaani : dataframe
            Dataframe with information on chain, residue number and pkaani output. Outline:
            Chain   Resid   pkaani
            x       x       x

        Example
        -------
        >> print(calculate_pkaani(1ubq.pdb))
        Chain  Resid  sasa
        0     A      6   x
        1     A     11   x
        2     A     27   x
        3     A     29   x
        4     A     33   x
        5     A     48   x
        6     A     63   x
        '''
        code_for_df = os.path.basename(path).split(".")[0]
        pdb_path = path.split(".")[0]
        test_path = pdb_path + '_pka.log'
        if not os.path.isfile(test_path):
            try:
                _ = subprocess.run(['pkaani', '-i', path])
            except Exception as e:
                print(e)
                if "[Errno 2] No such file or directory: 'pkaani'" == str(e):
                    if self.report_errors: self._report_error_to_file('pkaANI 1', path, str(e))
                    raise Exception(f'Error: Failed to obtain pkaANI data: {e}. Is pkaANI installed correctly? If so, reload environment and try again.')
                else:
                    if self.report_errors: self._report_error_to_file('pkaANI 1', path, str(e))
                    raise Exception(f"Failed to obtain pkaANI data through running pkaANI, error: {e}")

        try:
            log_file = pdb_path + '_pka.log'
            propres = open(log_file)
        except Exception as e:
            if self.report_errors: self._report_error_to_file('pkaANI 2', path, str(e))
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
            if self.report_errors: self._report_error_to_file('pkaANI 3', path, str(e))
            raise Exception(f'Failure parsing pkaANI log file for: {code_for_df}_pka.log. {e}') from e

        try:
            df_pkaani = pd.DataFrame({'Resid':lys_number,
                                      'Chain': chains,
                                      'pkaANI':pkas})

            df_pkaani.sort_values(by=['pkaANI'], inplace=True)
            df_pkaani = df_pkaani.dropna()
            df_pkaani = df_pkaani.drop_duplicates(subset=None, keep='first', inplace=False, ignore_index=True)
        except Exception as e:
            if self.report_errors: self._report_error_to_file('pkaANI 4', path, str(e))
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
        df_sasa : dataframe
            Dataframe with information on chain, residue number and sasa output. Outline:
            Chain   Resid   sasa
            x       x       x

        Example
        -------
        >> print(calculate_sasa(1ubq.pdb))
        Chain  Resid  sasa
        0     A      6   x
        1     A     11   x
        2     A     27   x
        3     A     29   x
        4     A     33   x
        5     A     48   x
        6     A     63   x
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
            if self.report_errors: self._report_error_to_file('SASA 1', path, str(e))
            raise Exception(f'SASA calc error: {e}') from e

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
                print(f'Error obtaining SASA at index value {str(j)}')
                list_of_sasa.append(None)
                if self.report_errors: self._report_error_to_file('Depth 1', path, f'Error obtaining SASA at index value {str(j)}')
                continue

        #append results to a df which is given as output
        try:
            df_sasa = pd.DataFrame({'Chain': list_of_chains,
                                'Resid': list_of_resid,
                                'sasa': list_of_sasa})

        except Exception as e:
            if self.report_errors: self._report_error_to_file('SASA 3', path, str(e))
            raise Exception(f'Error obtaining SASA data. {e}')

        return df_sasa


    def calculate_depth(self, path):
        '''
        Calculate the depth of the NZ atom from the surface of the protein within
        the overall protein structure.

        Method
        ------
        Identify all NZ atoms within the protein structure through biobox. Use the PDBparser
        from biopython to calculate the surface of the protein. Loop over the identified positions
        of all of the NZ atoms and use the residue_depth function contained in biopython
        to extract the minimum distances from the surface for each of the lysines.


        Parameters
        ----------
        path : string
            The path of the pdb file that the depth of the lysines are being calculated for.

        Returns
        -------
        df_depth : dataframe
            Dataframe with information on chain, residue number and depth output. Outline:
            Chain   Resid   depth
            x       x       x

        Example
        -------
        >> print(calculate_depth(1ubq.pdb))
        Chain  Resid  depth
        0     A      6   x
        1     A     11   x
        2     A     27   x
        3     A     29   x
        4     A     33   x
        5     A     48   x
        6     A     63   x
        '''
        try:
            M = bb.Molecule(path)
            pos, idx = M.atomselect("*", "*", "NZ", get_index=True)
        except Exception as e:
            if self.report_errors: self._report_error_to_file('Depth 1', path, str(e))
            raise Exception(f">> DEPTH error: could not find NZ atoms within atomic structure - {e}")

        try:
            parser = PDBParser()
            structure = parser.get_structure('structure', path)
            surface = get_surface(structure[0])
        except Exception as e:
            if self.report_errors: self._report_error_to_file('Depth 2', path, str(e))
            raise Exception(f">> DEPTH error: could not get biopython structure - {e}")


        results = []
        for i, coord in enumerate(pos):
            chain = M.data.loc[idx[i], ["chain"]].values[0]
            resid = M.data.loc[idx[i], ["resid"]].values[0]
            mychain = structure[0][chain]
            myres = mychain[int(resid)]
            try:
                #dist = min_dist(pos[i], surface)
                rd = residue_depth(myres, surface)
            except Exception as e:
                if self.report_errors: self._report_error_to_file('Depth 3', path, str(e))
                raise Exception(f">> DEPTH error: failed getting min_dist - {e}")

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
            if self.report_errors: self._report_error_to_file('AEV 1', path, str(e))
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
                    if self.report_errors: self._report_error_to_file('AEV 1.1', path, str(e))
                    print(f'AEV Calculations: could not create AEV for resid {idx_nz[j]} of protein {path}, error: {e}')

                # 2.3: Append the new AEV to the output dataframe
                aev_to_append = {'Chain': chain, 'Resid': resid, 'aev': aevs}
                df_aevs = pd.concat([df_aevs, pd.DataFrame([aev_to_append])], ignore_index=True)

        except torch.cuda.OutOfMemoryError:
            # potential that calculating the AEVs could overload the gpu, if too much memory, catch this and skip the file
            print(f'AEV calc error: CUDA memory error with file: {path}, skipping')
            if self.report_errors: self._report_error_to_file('AEV 2', path, 'CUDA memory error with file')
            return df_aevs
        except MemoryError:
            # potential that calculating the AEVs could overload the cpu, if too much memory, catch this and skip the file
            print(f'AEV calc error: CPU memory error with file: {path}, skipping')
            if self.report_errors: self._report_error_to_file('AEV 2', path, 'CPU memory error with file')
            return df_aevs
        except Exception as e:
            print(f'AEV Calculations: 2 - could not create the AEVs for the protein for protein {path}, error: {e}')
            if self.report_errors: self._report_error_to_file('AEV 2', path, str(e))
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
            if self.report_errors: self._report_error_to_file('DAS 1', path, e)
            print(f'DAS Calculation: 1 - could not load and identify the NZ atoms within the lysines of the structure: {e}')

        # 2: Setup the Xlink module and create the half spheres
        try:
            XL = bb.Xlink(M)
            das_output = []
        except Exception as e:
            if self.report_errors: self._report_error_to_file('DAS 2', path, str(e))
            print(f'DAS Calculation: 2 - Failed to setup the Xlink biobox class: {e}')

        for i, lys_nz_idx in enumerate(idx_nz):
            try:
                # the parameteres (pts_surf, thresh, radii) for the _get_half_sphere are already set
                # for lysine residues therefore the only parameter that needs to be set is i: this
                # is the index of the atom of interest within the lysine
                half_sphere_coords = XL._get_half_sphere(i=lys_nz_idx)
                # as the density of points created by the get half sphere is constant for any setup,
                # therefore can just count the number of coordinates that are returned for a measure for SASA Path
                das_output.append(len(half_sphere_coords))
            except Exception as e:
                if self.report_errors: self._report_error_to_file('DAS 2', path, str(e))
                print(f'DAS Calculation: 2 - Failed to calculate the half spheres for the NZ atoms on lysine no {lys_res_nums[i]}: {e}')
                das_output.append(None)

        # 3: Create dataframe to return
        df_das = pd.DataFrame(columns=["Chain", "Resid", "das"])
        try:
            df_das['Chain'] = list_chains
            df_das['Resid'] = lys_res_nums
            df_das['das'] = das_output
        except Exception as e:
            if self.report_errors: self._report_error_to_file('DAS 3', path, str(e))
            print(f'DAS Calculation: 3 - Failed to create datafame to append to the overall dataframe: {e}')

        return df_das
    

    def calculate_seqcharge(self, path, num_add_aa=10):
        '''
        Calculate the Sequence Charge of the local sequence around a LYS of interest.
        This is a single value number representing the summation of the charges of the amino acids
        over the specified number of amino acids either side of the lysine.

        Method
        ------
        Take the PDB file and extract the overall sequence using BioBox. Identify all lysines
        within the structure and any shift that has taken place in the PDB file compared to
        the Uniprot sequence. Extract sequences for the lysine of interest and calculate a
        value for the charge based on the summation of charged residues within the sequence.

        Parameters
        ----------
        path : string
            The path of the pdb file that DAS is being calculated for.

        num_add_aa : int
            The number of amino acids to include either side of the lysine of interest.
            An optional parameter which is set to 10 by default.

        Returns
        -------
        df_seqcharge : dataframe
            Dataframe with information on chain, residue number and seqcharge output. Outline:
            Chain   Resid   seqcharge
            x           x           x

        Example
        -------
        >> print(self.calculate_seqcharge(1M2F-alt-1.pdb))
        Chain   Resid  seqcharge
        0     A      95         -3
        '''

        # 1: Extract the overall sequence for the protein given
        try:
            M = bb.Molecule(path)
            idx_nz = M.atomselect('*', 'LYS', 'NZ', use_resname=True, get_index=True)[1]
            lys_res_nums = list(M.data['resid'][idx_nz])
            list_chains = list(M.data['chain'][idx_nz])

            c_alpha_idxs = M.atomselect('*', '*', 'CA', use_resname=True, get_index=True)[1]

            # TODO GW 16.04.25 - eventually will need to add in ability to use  letter codes and charges
            #                    for the 3 lettter cases rather than the 1 letter cases which when using
            #                    modified residues may run into problems

            protein_letters_dict = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
                                    'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G',
                                    'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
                                    'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
                                    'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
                                    'HIE': 'H', 'HID': 'H', 'HIP': 'H', 'LYN': 'K',
                                    'ASX': 'B', 'GLX': 'Z', 'SEC': 'U', 'PYL': 'O',
                                    'XAA': 'X', 'XLE': 'J', 'PSER': 'p', 'PTHR': 't',
                                    'PTYR': 'y', 'MELYS': 'k', 'MEARG': 'r', 'ACLYS': 'k',
                                    'KCX': 'X', 'LYE': 'X'}  
            # KCX and LYE down as X so that they are not treated as positive K

            def _catch(func, *args, handle=lambda e : e, **kwargs):
                try:
                    return func(*args, **kwargs)
                except Exception as e:
                    print(f"Could not convert {e} to a 1 letter code: using 'X' instead")
                    return 'X'

            sequence = ''.join([_catch(lambda : protein_letters_dict[a.upper()]) for a in list(M.data['resname'][c_alpha_idxs])])


        except Exception as e:
            if self.report_errors: self._report_error_to_file('Seqcharge 1', path, str(e))
            print(f'SeqCharge Calculation: 1 - could not extract the sequence from the protein file given: {e}')
            return pd.DataFrame(columns=["Chain", "Resid", "seqcharge"])

        # 2: Extract local sequences based on the overall chain, calculate charge score and add to output
        seqcharge_output = []
        for i, lys_res_idx in enumerate(lys_res_nums):
            try:
                # get the index of the start position of the residues to work out the shift
                shift_val = int(M.data['resid'].iloc[0]) - 1
                seq_lys_index = lys_res_idx - shift_val - 1

                start_idx = seq_lys_index - num_add_aa
                end_idx = seq_lys_index + num_add_aa + 1
                start_null = 0
                end_null = 0

                if start_idx < 0:
                    start_null = - start_idx
                    start_idx = 0

                if end_idx >= len(sequence):
                    end_null = - end_idx
                    end_idx = len(sequence)

                seq = ('-' * start_null) + sequence[start_idx:end_idx] + ('-' * end_null)
                seq_split = list(seq)
                if seq_split[10] != 'K':
                    print(f'A lysine was not found at the desired position {lys_res_idx+1} read in for PDB file {path}; sequence -> {seq}')
                    continue

                pos_aa = ['K', 'H', 'R']
                neg_aa = ['D', 'E']
                count = 0
                for aa in seq:
                    if aa in pos_aa:
                        count += 1
                    elif aa in neg_aa:
                        count -= 1
                seqcharge_output.append(count)

            except Exception as e:
                if self.report_errors: self._report_error_to_file('Seqcharge 2', path, str(e))
                print(f'SeqCharge Calculation 2: Could not calculate a charge for lysine at position {lys_res_idx}, error: {e}')
                seqcharge_output.append(None)

        # 3: Create dataframe to return
        df_seqcharge = pd.DataFrame(columns=["Chain", "Resid", "seqcharge"])
        try:
            df_seqcharge['Chain'] = list_chains
            df_seqcharge['Resid'] = lys_res_nums
            df_seqcharge['seqcharge'] = seqcharge_output
        except Exception as e:
            if self.report_errors: self._report_error_to_file('Seqcharge 3', path, str(e))
            print(f'SeqCharge Calculation: 3 - Failed to create datafame to append to the overall dataframe: {e}')

        #print(df_seqcharge)
        return df_seqcharge


    def calculate_melodia(self, path):
        '''
        Call the Melodia package to calculate data for the following structural features of
        the lysines of interest within the structure: curvature, arc-length, phi, psi

        Method
        ------
        Call melodia on the path of the pdb file that has been passed to the function and
        create dataframe from the results. Process the dataframe to remove any of the 
        calculated features that were not asked for.

        Parameters
        ----------
        path : string
            The path of the pdb file that DAS is being calculated for.

        Returns
        -------
        df_melodia : dataframe
            Dataframe with information on chain, residue number and desired melodia output.
            Outline for all features:
            Chain   Resid   curvature   writhing    torsion   arc-length  phi psi
            x           x           x          x          x            x    x   x

        Example
        -------
        >> print(self.calculate_melodia(1ubq.pdb))
        Chain   Resid    curvature  writhing    torsion  arc-length  phi psi
        0
        '''
        # Melodia 1 - Calculating geometry using melodia-py
        try:
            melodia_results = mel.geometry_from_structure_file(path)
            if isinstance(melodia_results, pd.Series):
                melodia_results = melodia_results.to_frame().T

        except Exception as e:
            if self.report_errors: self._report_error_to_file('Melodia 1', path, str(e))
            print(f'Melodia 1: Error processing input file - {path} with error: {e}')

        # Melodia 2 - Formatting and filtering
        try:
            #melodia_results.rename({"chain": "Chain", "order": "Resid", "curvature": "melodia"}, axis="columns", inplace = True)
            melodia_results.rename({"chain": "Chain", "order": "Resid"}, axis="columns", inplace = True)
            lys_results = melodia_results['name'] == 'LYS'
            df_melodia = melodia_results[lys_results].copy()
            df_melodia.reset_index(inplace=True, drop=True)
            cols_to_drop = ['code', 'id', 'model', 'curvature', 'writhing', 'torsion', 'phi', 'psi', 'name', 'arc_length']
            cols_to_drop = [col for col in cols_to_drop if col not in self.features]
            df_melodia.drop(labels=cols_to_drop, axis = 'columns', inplace=True)

        except Exception as e:
            if self.report_errors: self._report_error_to_file('Melodia 2', path, str(e))
            print(f'Melodia 2: Unable to reformat melodia output correctly for input {path} with error: {e}')

        #print(df_melodia)
        return df_melodia


    def calculate_frustration(self, path):
        '''
        Use the Frustratometer package to identify the frustration metric for the lysines
        of interest.

        Method
        ------
        Load in protein structure into frustratometer before using AWSEM to create a model for
        this with desired parameters. Use this model to calculate the 

        Parameters
        ----------
        path : string
            The path of the pdb file that DAS is being calculated for.

        Returns
        -------
        df_frustration : dataframe
            Dataframe with information on chain, residue number and desired output from Frustratometer.
            Outline for all features:
            Chain   Resid   frustration   density
            x           x             x         x

        Example
        -------
        >> print(self.calculate_frustration(1ubq.pdb))
            PDB_Code    Chain   Resid   frustration     density
        0       1ubq        A       6     -1.203796    4.006602
        '''
        # temp - modules required to add into readme - openmm, pdbfixer
        # Frustratometer 1 - create structure and AWSEM model
        df_frustration = pd.DataFrame()
        try:
            M = bb.Molecule(path)
            idx_nz = M.atomselect('*', 'LYS', 'NZ', use_resname=True, get_index=True)[1]
            lys_res_nums = list(M.data['resid'][idx_nz])
            list_chains = list(M.data['chain'][idx_nz])
            df_frustration['Chain'] = list_chains
            df_frustration['Resid'] = lys_res_nums

            out_print_trap = io.StringIO()
            with redirect_stdout(out_print_trap):
                frust_struc = frustratometer.Structure(path)
                model_single_resids = frustratometer.AWSEM(frust_struc, min_sequence_separation_contact=2)
            del out_print_trap
        except Exception as e:
            print(f'Frustratometer calculation 1 - failed to create frustratometer structure or AWSEM model with error: {e}')
            df_frustration['frustration'] = None
            df_frustration['density'] = None
            if self.report_errors: self._report_error_to_file('Frustratometer 1', path, str(e))
            return df_frustration

        # Frustratometer 2 - use model to calculate outputs
        try:
            single_residue_awsem_frustration = model_single_resids.frustration(kind='singleresidue')
            resid_densities = model_single_resids.rho_r
        except Exception as e:
            if self.report_errors: self._report_error_to_file('Frustratometer 2', path, str(e))
            print(f'Frustratometer calculation 2 - failed to create frustratometer outputs: {e}')

        # Frustration 3 - extract lysine values from the outputs and append to output dataframe
        try:
            lys_frustration = []
            lys_density = []
            for lys in lys_res_nums:
                lys_frustration.append(single_residue_awsem_frustration[lys-1])
                lys_density.append(resid_densities[lys-1])
            df_frustration['frustration'] = lys_frustration
            df_frustration['density'] = lys_density

            try:
                pdb_code = path.split('/')[-1]
                cleaned_code_to_remove = pdb_code.split('.')[0] + '_cleaned.pdb'
                os.remove(cleaned_code_to_remove)
            except Exception as ef:
                print(f'Failed to remove cleaned pdb for frustratometer calculation with error {ef}')
        except Exception as e:
            if self.report_errors: self._report_error_to_file('Frustratometer 3', path, str(e))
            print(f'Frustratometer calculation 3 - failed to append data to return dataframe: {e}')

        #print(df_frustration)
        return df_frustration


    def calculate_legolas(self, path):
        '''
        Calculate 15N nmr data using legolas
        
        Method
        ------
        Use biobox to extract the positions of the lysines within the the protein
        structure given in the path. Then change directory to the path of legolas
        and run legolas.py on the desired protein structure.

        Parameters
        ----------
        path : string
            The path of the pdb file that legolas is being calculated for.
        
        Example
        -------
        >> print(self.calculate_legolas(1ubq.pdb))
                             PDB_Code   Chain  Resid   legolas
        0     data/curated/1UBQ-alt-1       A      6   121.614
        1     data/curated/1UBQ-alt-1       A     11   121.192
        2     data/curated/1UBQ-alt-1       A     27   118.507
        3     data/curated/1UBQ-alt-1       A     29   119.557
        4     data/curated/1UBQ-alt-1       A     33   117.238
        5     data/curated/1UBQ-alt-1       A     48   119.989
        6     data/curated/1UBQ-alt-1       A     63   121.946
        '''

        # 1: Load in the structure and locate all the NZ atoms within the lysines, calculate the list of chains and list of resids to go with this
        try:
            M = bb.Molecule(path)
            idx_nz = M.atomselect('*', 'LYS', 'NZ', use_resname=True, get_index=True)[1]
            lys_res_nums = list(M.data['resid'][idx_nz])
            list_chains = list(M.data['chain'][idx_nz])
            result_filename = os.path.join(self.legolas_output_path, path.split(f'{os.sep}')[-1].split('.')[0] + '_cs.csv')
            if os.path.exists(result_filename):
                already_exists = True
                print(f'Legolas: The file {result_filename} already exists, using this file')
            else:
                already_exists = False
            modified_struc = False
            if not already_exists:
                if 'HIE' in list(M.data['resname']):
                    print('modifying the structure to change HIE to HIS')
                    M.data.loc[M.data['resname'] == 'HIE', 'resname'] = 'HIS'
                    M.write_pdb('temp_legolas.pdb')
                    modified_struc = True
        except Exception as e:
            print(f'Legolas: 1 - could not load and identify the NZ atoms within the lysines of the structure: {e}')
            if self.report_errors: self._report_error_to_file('LEGOLAS 1', path, str(e))
            return pd.DataFrame(columns=['Chain', 'Resid', 'legolas'])

        # 2: change location to legolas directory and run the legolas program on the specified pdb before changing back to working directory
        try:
            if modified_struc:
                pdb_absolute_path = 'temp_legolas.pdb'
                result_filename = 'temp_legolas_cs.csv'
            else:
                pdb_absolute_path = os.path.join(os.getcwd(), path)
                result_filename = path.split(f'{os.sep}')[-1].split('.')[0] + '_cs.csv'
            if not already_exists:
                legolas_prog = '/home/gweston/Documents/extra_packages/legolas-main/test/legolas.py'
                subprocess.run(['python', legolas_prog, pdb_absolute_path, '-atype', 'N'])
            else:
                result_filename = os.path.join(self.legolas_output_path, result_filename)
            df_nmr = pd.read_csv(result_filename)
            nmr_results_list = list(df_nmr['CHEMICAL_SHIFT'])
            lys_nmr_vals = []
            for idx in lys_res_nums:
                lys_nmr_vals.append(nmr_results_list[(idx - 1)])
            #os.remove(result_filename)
            if not already_exists:
                os.remove(result_filename.split('.')[0] + '.parquet')
                os.rename(result_filename, path.split(f'{os.sep}')[-1].split('.')[0] + '_cs.csv')
                result_filename = path.split(f'{os.sep}')[-1].split('.')[0] + '_cs.csv'
                shutil.move(result_filename, self.legolas_output_path)
                if modified_struc:
                    os.remove('temp_legolas.pdb')

        except Exception as e:
            print(f'Legolas 2: Failed to run the legolas program and extract the 15N nmr shifts for the protein: {e}')
            if self.report_errors: self._report_error_to_file('LEGOLAS 2', path, str(e))
            return pd.DataFrame(columns=['Chain', 'Resid', 'legolas'])

        # 3: Create dataframe to return
        df_legolas = pd.DataFrame(columns=["Chain", "Resid", "legolas"])
        try:
            df_legolas['Chain'] = list_chains
            df_legolas['Resid'] = lys_res_nums
            df_legolas['legolas'] = lys_nmr_vals
        except Exception as e:
            if self.report_errors: self._report_error_to_file('LEGOLAS 3', path, str(e))
            print(f'Legolas: 3 - Failed to create datafame to append to the overall dataframe: {e}')

        return df_legolas

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
    M.measure_dataframe_parallel()
