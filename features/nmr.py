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

    def __init__(self, outdir, legolas_aevs=False, include_modified=False,
                 aa_properties = {'non_modified_codes': ['LYS', 'LYSN'],
                                'modified_codes': ['LYE', 'KCX'],
                                'atom_select_names_nonmod': ['NZ'],
                                'atom_select_names_modified': ['NZ', 'N07']},
                 error_filename = 'measure_errors.txt'):
        '''
        Initialise the NMR class

        Parameters
        ----------
        outdir : string
            The output directory that measurements will be saved to.
        legolas_aevs : bool
            Toggle setting to indicate if you want the legolas programme to dump the AEVs from
            the calculation of the 15N nmr values. Default is False.
        include_modified : bool
            Toggle to include residues which have been modified within the featurisation
        aa_properties -> dict
            Properties of the amino acid of interest to investigate modification sites for.
            Defaults to lysine for carbamylation. Properties are the 3 letter codes for
            non modified ('non_modified_codes') and modified ('modified_codes') and the atom
            names for non modified ('atom_select_names_nonmod') and modified ('atom_select_names_modified')
        error_filename : str
            Name of the text file passed through from overall measures to write any errors from
            calculating features out to.
        '''

        self.outdir = outdir
        self.legolas_aevs = legolas_aevs
        self.include_modified = include_modified
        self.aa_properties = aa_properties
        self.error_filename = error_filename

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

            M_n = M.get_subset(M.atomselect('*', '*', 'N', use_resname=True, get_index=True)[1])
            if self.include_modified: idx_n_res_interest = M_n.atomselect('*', (self.aa_properties['non_modified_codes'] + self.aa_properties['modified_codes']), 'N', use_resname=True, get_index=True)[1]
            else: idx_n_res_interest = M_n.atomselect('*', self.aa_properties['non_modified_codes'], 'N', use_resname=True, get_index=True)[1]
            list_modified = list(a in self.aa_properties['modified_codes'] for a in list(M_n.data['resname']))
            df_legolas = M_n.data[['resid', 'chain', 'resname']]
            if self.include_modified: df_legolas = df_legolas.assign(**{'Modified': list_modified})

        except Exception as e:
            print(f'Legolas: 1 - could not load and identify the NZ atoms within the lysines of the structure: {e}')
            report_error_to_file('LEGOLAS 1', path, str(e), self.error_filename)
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
                # Note: if you are not GW and running this, you will need to change this path to your own installation path!!
                legolas_prog = '/home/gweston/Documents/extra_packages/legolas-main/test/legolas.py'
                subprocess.run(['python', legolas_prog, pdb_absolute_path, '-atype', 'N'])
            else:
                result_filename = os.path.join(self.legolas_output_path, result_filename)
            df_nmr = pd.read_csv(result_filename)
            df_legolas = pd.concat([df_legolas, df_nmr.loc[:, ~df_nmr.columns.str.contains('^Unnamed')]], axis=1)
            #print(df_legolas_new)

            if not already_exists:
                os.remove(result_filename.split('.')[0] + '.parquet')
                os.rename(result_filename, path.split(f'{os.sep}')[-1].split('.')[0] + '_cs.csv')
                result_filename = path.split(f'{os.sep}')[-1].split('.')[0] + '_cs.csv'
                shutil.move(result_filename, self.legolas_output_path)
                if modified_struc:
                    os.remove('temp_legolas.pdb')

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

                    if aevs != []:
                        df_legolas = df_legolas.assign(**{'legolas_aev': [str(aevs[i]) for i in list(M.atomselect('*', '*', 'N', use_resname=True, get_index=True)[1])]})

                    if os.path.exists('tmp_aevs_protein.txt'):
                        os.remove('tmp_aevs_protein.txt')

                except Exception as e:
                    print(f'Legolas AEVs: failed to extract aev data from legolas: {e}')
                    report_error_to_file('LEGOLAS AEV 1', path, str(e))
            
            df_legolas = df_legolas.iloc[idx_n_res_interest].drop(columns=['resname', 'ATOM_TYPE', 'RESIDUE_ID', 'CHEMICAL_SHIFT_STD']).rename(columns={'CHEMICAL_SHIFT': 'legolas', 'chain': 'Chain', 'resid': 'Resid'})

        except Exception as e:
            print(f'Legolas 2: Failed to run the legolas program and extract the 15N nmr shifts for the protein: {e}')
            report_error_to_file('LEGOLAS 2', path, str(e), self.error_filename)
            return pd.DataFrame(columns=['Chain', 'Resid', 'legolas'])

        return df_legolas


if __name__ == '__main__':
    nmr = NMR(outdir='result', legolas_aevs=True, include_modified=False)
    #print(nmr.calculate_legolas(path=f'result{os.sep}curated{os.sep}1UBQ-alt-1.pdb'))
    print(nmr.calculate_legolas(path=f'2I1V-alt-1.pdb'))
    #print(nmr.calculate_legolas(path=f'1nsk_AmberMod0000.pdb'))
