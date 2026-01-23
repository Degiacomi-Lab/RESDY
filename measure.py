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
from ast import literal_eval
import pandas as pd
import numpy as np
import biobox as bb
import matplotlib.pyplot as plt
#import dill
#import features
from features.aev import AEV
from features.charge import Charge
from features.das import DAS
from features.depth import Depth
from features.frustration import Frustration
from features.nmr import NMR
from features.pka import PKA
from features.sasa import SASA
from features.structure import Structure
from features.flexibility import Flexibility


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
    '''
    Class to handle functions used in calling feature functions and managing how these are
    called and return a dataframe which contains the results after.
    '''

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
        self.legolas_aevs = True
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
        if self.report_errors: self._setup_report_errors_file()

        # modified lysine management
        self.include_modified = include_modified

        # for parallel measurements
        self.parallel = parallel
        self.files_to_analyse = []
        self.parallel_items = {}

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
        Convert a list of features into a measuring protocol. If ['all'] given as input for
        the features, this will convert the features list to a list containing all current
        possible features.
        
        Parameters
        ----------
        features : list
            The list of features that are required to measure over the set of proteins
        '''

        # measures to carry out [label for DataFrame column, and function evaluating a file]
        # functions must return a dataframe [chain, resid, measure]
        if 'all' in features:
            features = ['propka', 'pkaANI', 'sasa', 'depth', 'aev', 'seqcharge', 'legolas',
                        'melodia', 'aev_legolas', 'frustration', 'density', 'das']
            self.features = features
        if 'melodia' in self.features: self.features.append(self.features.pop(self.features.index('melodia')))
        self.measures = []
        melodia_features = []
        melodia_added = False
        frustration_added = False
        legolas_added = False
        #TODO can this handle calculations where you actually want multiple methods for same feature calculating
        for m in features:
            if m in ['propka', 'pkaANI']:
                pka = PKA(outdir=self.outdir, calc_method=m)
                self.measures.append([m, pka.calculate_pka])
            elif m == 'pka':
                print('Please enter which pKa calculation method you would like to use: ' \
                      'propka or pkaANI')
                while not input('propka or pkaANI:') in ['propka', 'pkaANI']:
                    print('Please enter either propka or pkaANI')  #TODO this needs testing but should work
                pka = PKA(outdir=self.outdir, calc_method='propka')
                self.measures.append([m, pka.calculate_propka])
            elif m == 'sasa':
                sasa = SASA()
                self.measures.append([m, sasa.calculate_sasa])
            elif m == "depth":
                depth = Depth()
                self.measures.append([m, depth.calculate_depth])
            elif m == 'aev':
                aev = AEV()
                self.measures.append([m, aev.calculate_aevs])
            elif m == 'das':
                das = DAS()
                self.measures.append([m, das.calculate_das])
            elif m == 'seqcharge':
                charge = Charge()
                self.measures.append([m, charge.calculate_seqcharge])
            elif m == 'flexibility':
                flex = Flexibility()
                self.measures.append([m, flex.calculate_flexibility])
            elif m == 'legolas':
                if self.legolas_aevs:
                    if 'aev_legolas' not in self.features:
                        self.features.append('aev_legolas')
                    nmr = NMR(outdir=self.outdir, legolas_aevs=True)
                else:
                    nmr = NMR(outdir=self.outdir, legolas_aevs=False)
                self.measures.append([m, nmr.calculate_legolas])
                legolas_added = True
            elif m == 'aev_legolas':
                if not legolas_added:
                    nmr = NMR(outdir=self.outdir, legolas_aevs=True)
                    self.measures.append(['legolas', nmr.calculate_legolas])
                    legolas_added = True
            elif m in ['frustration', 'density']:
                if not frustration_added:
                    frustration = Frustration()
                    self.measures.append(['frustration', frustration.calculate_frustration])
                    frustration_added = True
            elif m == 'melodia':
                structure = Structure(melodia_features=['all'])
                self.measures.append([m, structure.calculate_melodia])
                melodia_added = True
                self.features += ['curvature', 'writhing', 'torsion', 'arc_length', 'phi', 'psi']
                self.features.remove('melodia')
            elif m in ['curvature', 'writhing', 'torsion', 'arc_length', 'phi', 'psi']:
                if not melodia_added:
                    melodia_features += [m]
                    structure = Structure(melodia_features=melodia_features)
                    self.measures.append(['melodia', structure.calculate_melodia])
                    melodia_added = True
            else:
                raise Exception(f"measure {m} unknown")

    def _setup_report_errors_file(self):
        '''
        Function to set up the file where errors produced through running the Measure
        class will be written to such that they are easier to look over after running,
        rather than trawling through output.
        '''
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


    def restart_measure_data(self):
        '''
        Determine the appropriate measures function to call based on the combination of
        running PDB_only and in parallel, reducing the number individual functions that
        the user will have to call themselves. Different to measure_data() as this will
        restart the measurements from final previous point rather than starting again.

        Example
        -------
        M.restart_measure_data()
        '''
        match (self.PDB_only, self.parallel):
            case (False, False) | (False, True):
                # (not PDB only and not parallel) or (not PDB only and parallel):
                self.restart_measure()
            case (True, False):
                # PDB only and not parallel
                self.restart_measure_pdb_only()
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
        Will take a list of the required proteins, finds associated curated structures and runs the required measurement functions.
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

        # remove possible duplicated rows in dataframe
        try:
            self.df.drop_duplicates(subset=None, keep='first', inplace=True, ignore_index=True)
            print('\n>> Removed duplicates from measurement dataframe.')
        except Exception as e_one:
            try:
                print(f'\n>> Failed to remove duplicates from measurement dataframe, trying new method: {e_one}')
                self.df = self.df.astype(str).drop_duplicates(subset=None, keep='first', inplace=False, ignore_index=True)
                print('>> Removed duplicates from measurement dataframe using new method.')
            except Exception as e_two:
                print(f'>> Failed to remove duplicates from measurement dataframe: {e_two}')
                pass


    def measure_dataframe(self):
        '''
        TODO FINISH THIS
        Function to measure specified features for all the structure files curated earlier in the programme.
        Will take a list of the required proteins, finds associated curated structures and runs the required measurement functions.
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

            if i != 0:
                avg_time_per_pdb = (time.time() - overall_st) / i
                pred_time_remaining = avg_time_per_pdb * (len(self.df_input) - i)
            else:
                pred_time_remaining = 'undefined'
            print(f'Analysing PDB code ({pdb_code}) {i}/{len(self.df_input)}. Predicted time remaining: {pred_time_remaining}')

            # calculate features values from all PDB files associated with specific DataFrame entry
            for f in files_list:
                # check the file for the required pbd code, if not there, skip
                if pdb_code not in f:
                    continue

                tstart = time.time()
                print(f"\n> Calculating for measurements for file: {f}")

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
                        print(f"Error iterating measures, potential dataframe combination problem: {e}")
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

        # remove possible duplicated rows in dataframe
        try:
            self.df.drop_duplicates(subset=None, keep='first', inplace=True, ignore_index=True)
            print('\n>> Removed duplicates from measurement dataframe.')
        except Exception as e_one:
            try:
                print(f'\n>> Failed to remove duplicates from measurement dataframe, trying new method: {e_one}')
                self.df = self.df.astype(str).drop_duplicates(subset=None, keep='first', inplace=False, ignore_index=True)
                print('>> Removed duplicates from measurement dataframe using new method.')
            except Exception as e_two:
                print(f'>> Failed to remove duplicates from measurement dataframe: {e_two}')
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
        if self.parallel:
            self.measure_dataframe_parallel()
        else:
            self.measure_dataframe()


    def _combine_dataframes(self, target, to_merge, col_name):
        '''
        Function to combine the dataframe produced by a measurement function into the
        main dataframe containing all the measurements.
        target is a DataFrame to be filled with data, to_merge contains the data.
        Values to insert are indexed in both array by two columns: Chain and Resid.

        Parameters
        ----------
        target : DataFrame
            DataFrame to be filled with data.
        to_merge : DataFrame
            to_merge contains the new data to merge.
        col_name : string
            Name of the column which the new data is from.
        
        Example
        -------
        self._combine_dataframes(df, result, meas[0])
        '''
        for i, r in target.iterrows():

            chain_value = r["Chain"]
            resid_value = r["Resid"]

            idx = np.where((to_merge["Chain"] == chain_value) & (to_merge["Resid"].astype(int) == resid_value))
            if len(idx[0]) == 0:
                continue

            # account for measurements that have special cases
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
            elif col_name == 'legolas':
                legolas_features = ['legolas', 'aev_legolas']
                for feature in self.features:
                    if feature in legolas_features:
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
        print('\n>> Cleaning up leftover files from measures calculations...')
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
        print('>> Unused file cleanup complete.')


if __name__ == "__main__":


    file_one = "Demo{os.sep}curated{os.sep}1M2E-alt-1.pdb"

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
