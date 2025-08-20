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

# Melodia packages
try:
    import melodia_py as mel
except Exception as e:
    print(f"melodia unavailable. Unable to calculate melodia. Error: {e}")

class Structure():
    '''
    Class to house the different methods for calculating structural measurement values
    for structures.

    Parameters
    ----------

    '''

    def __init__(self, melodia_features=['all']):
        '''
        Initialise the Structure class
        
        Parameters
        ----------
        melodia_features : list
            List of features which are calculated through melodia which has been requested
            when the Measure class is initialised. Default is set to ['all'].
        '''

        self.melodia_features = melodia_features
        if self.melodia_features == ['all']:
            self.melodia_features = ['curvature', 'writhing', 'torsion', 'arc_length', 'phi', 'psi']

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
            report_error_to_file('Melodia 1', path, str(e))
            print(f'Melodia 1: Error processing input file - {path} with error: {e}')

        # Melodia 2 - Formatting and filtering
        try:
            melodia_results.rename({"chain": "Chain", "order": "Resid"}, axis="columns", inplace = True)
            lys_results = melodia_results['name'] == 'LYS'
            df_melodia = melodia_results[lys_results].copy()
            df_melodia.reset_index(inplace=True, drop=True)
            cols_to_drop = ['code', 'id', 'model', 'curvature', 'writhing',
                            'torsion', 'phi', 'psi', 'name', 'arc_length']
            cols_to_drop = [col for col in cols_to_drop if col not in self.melodia_features]
            df_melodia.drop(labels=cols_to_drop, axis = 'columns', inplace=True)

        except Exception as e:
            report_error_to_file('Melodia 2', path, str(e))
            print(f'Melodia 2: Unable to reformat melodia output correctly for input {path} with error: {e}')

        return df_melodia
