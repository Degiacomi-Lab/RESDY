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

class DAS():
    '''
    Class to house the different methods for calculating dynamically accessible surface
    area (das) values for structures
    '''
    
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