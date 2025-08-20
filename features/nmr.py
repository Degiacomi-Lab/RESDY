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
from features.error_reporting import report_error_to_file


# AEV packages
try:
    from ase import Atoms
    import torch
    import torchani
except Exception as e:
    print(f'Packages required for AEV calculation are not available, will not be able to calculate AEVs. Error: {e}')


class NMR():
    '''
    Class to house the different methods for calculating 15N nmr values for structures
    '''

    def __init__(self, outdir, legolas_aevs=False):
        '''
        Initialise the NMR class

        Parameters
        ----------
        outdir : string
            The output directory that measurements will be saved to.

        legolas_aevs : bool
            Toggle setting to indicate if you want the legolas programme to dump the AEVs from
            the calculation of the 15N nmr values. Default is False.
        '''

        self.outdir = outdir
        self.legolas_aevs = legolas_aevs

        self.legolas_output_path = os.path.join(self.outdir, 'legolas')
        if not os.path.exists(self.legolas_output_path):
            os.mkdir(self.legolas_output_path)

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
            report_error_to_file('LEGOLAS 1', path, str(e))
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

            #  AEV section from this -if include_aev:
            if self.legolas_aevs:
                try:
                    new_aev_filename = f'{path.split(f"{os.sep}")[-1].split(".")[0]}_legolasaev.txt'
                    if already_exists:
                        if os.path.exists(os.path.join(self.legolas_output_path, new_aev_filename)):
                            with open(os.path.join(self.legolas_output_path, new_aev_filename)) as f:
                                aevs = f.readlines()
                        else:
                            aevs = []
                    elif os.path.exists('tmp_aevs_protein.txt'):
                        with open('tmp_aevs_protein.txt', 'r') as f:
                            aevs = f.readlines()
                        os.rename('tmp_aevs_protein.txt', new_aev_filename)
                        shutil.move(new_aev_filename, self.legolas_output_path)
                    else:
                        print('Legolas AEVs: tmp_aevs_protein.txt file not found')
                        aevs = []


                    lys_aevs = []
                    for idx in lys_res_nums:
                        #lys_aevs.append(literal_eval(aevs[(idx - 1)]))
                        lys_aevs.append(aevs[(idx - 1)])


                    if os.path.exists('tmp_aevs_protein.txt'):
                        os.remove('tmp_aevs_protein.txt')
                    #print(lys_aevs)

                except Exception as e:
                    print(f'Legolas AEVs: failed to extract aev data from legolas: {e}')
                    report_error_to_file('LEGOLAS AEV 1', path, str(e))


        except Exception as e:
            print(f'Legolas 2: Failed to run the legolas program and extract the 15N nmr shifts for the protein: {e}')
            report_error_to_file('LEGOLAS 2', path, str(e))
            return pd.DataFrame(columns=['Chain', 'Resid', 'legolas'])

        # 3: Create dataframe to return
        df_legolas = pd.DataFrame(columns=["Chain", "Resid", "legolas"])
        try:
            df_legolas['Chain'] = list_chains
            df_legolas['Resid'] = lys_res_nums
            df_legolas['legolas'] = lys_nmr_vals
            if self.legolas_aevs:
                df_legolas['aev_legolas'] = str(lys_aevs)
        except Exception as e:
            report_error_to_file('LEGOLAS 3', path, str(e))
            print(f'Legolas: 3 - Failed to create datafame to append to the overall dataframe: {e}')

        return df_legolas
    