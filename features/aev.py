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
from error_reporting import _report_error_to_file


# AEV packages
try:
    from ase import Atoms
    import torch
    import torchani
except Exception as e:
    print(f'Packages required for AEV calculation are not available, will not be able to calculate AEVs. Error: {e}')


class AEV():
    '''
    Class to house the different methods for calculating representations for the local structure
    of lysines within the protein structure termed atomic environment vectors (AEVs). Potential AEV
    representations are currently:
    1. ANI-2x AEVs
    2. LEGOLAS ANI-2x AEVs
    3. Coarse-grain representation AEVs 
    '''
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
            _report_error_to_file('AEV 1', path, str(e))
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
                    _report_error_to_file('AEV 1.1', path, str(e))
                    print(f'AEV Calculations: could not create AEV for resid {idx_nz[j]} of protein {path}, error: {e}')

                # 2.3: Append the new AEV to the output dataframe
                aev_to_append = {'Chain': chain, 'Resid': resid, 'aev': aevs}
                df_aevs = pd.concat([df_aevs, pd.DataFrame([aev_to_append])], ignore_index=True)

        except torch.cuda.OutOfMemoryError:
            # potential that calculating the AEVs could overload the gpu, if too much memory, catch this and skip the file
            print(f'AEV calc error: CUDA memory error with file: {path}, skipping')
            _report_error_to_file('AEV 2', path, 'CUDA memory error with file')
            return df_aevs
        except MemoryError:
            # potential that calculating the AEVs could overload the cpu, if too much memory, catch this and skip the file
            print(f'AEV calc error: CPU memory error with file: {path}, skipping')
            _report_error_to_file('AEV 2', path, 'CPU memory error with file')
            return df_aevs
        except Exception as e:
            print(f'AEV Calculations: 2 - could not create the AEVs for the protein for protein {path}, error: {e}')
            _report_error_to_file('AEV 2', path, str(e))
            return df_aevs

        # 4: if everything has worked, return the dataframe with the AEVs for the protein
        #print(df_aevs)
        return df_aevs
