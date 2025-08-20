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


class SASA():
    '''
    Class to house the different methods for calculating solvent accessible surface area
    (SASA) values for structures
    '''

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
            report_error_to_file('SASA 1', path, str(e))
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
                report_error_to_file('Depth 1', path, f'Error obtaining SASA at index value {str(j)}')
                continue

        #append results to a df which is given as output
        try:
            df_sasa = pd.DataFrame({'Chain': list_of_chains,
                                'Resid': list_of_resid,
                                'sasa': list_of_sasa})

        except Exception as e:
            report_error_to_file('SASA 3', path, str(e))
            raise Exception(f'Error obtaining SASA data. {e}')

        return df_sasa
