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

    def __init__(self, include_modified = False, aa_properties = {'non_modified_codes': ['LYS', 'LYSN'],
                                                                'modified_codes': ['LYE', 'KCX'],
                                                                'atom_select_names_nonmod': ['NZ'],
                                                                'atom_select_names_modified': ['NZ', 'N07']},
                  error_filename = 'measure_errors.txt'):
        '''
        Initialise the SASA class, include any global variables that are required from
        measures in here.

        Parameters
        ----------
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
        self.include_modified = include_modified
        self.aa_properties = aa_properties
        self.error_filename = error_filename

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
            list_of_sasa = []

            M = bb.Molecule()
            M.import_pdb(path, include_hetatm=True)

            if self.include_modified: lys_coords, lys_idx = M.atomselect('*', (self.aa_properties['non_modified_codes'] + self.aa_properties['modified_codes']), self.aa_properties['atom_select_names_modified'], use_resname=True, get_index=True)
            else: lys_coords, lys_idx = M.atomselect('*', self.aa_properties['non_modified_codes'], self.aa_properties['atom_select_names_nomod'], use_resname=True, get_index=True)

            list_of_resid = list(M.data['resid'][lys_idx])
            list_of_chains = list(M.data['chain'][lys_idx])
            list_modified = list(a in self.aa_properties['modified_codes'] for a in list(M.data['resname'][lys_idx]))

            all_coords, all_idx = M.atomselect('*','*','*', get_index=True)

        except Exception as e:
            report_error_to_file('SASA 1', path, str(e), self.error_filename)
            print(f'SASA Calculation: 1 - Failed to extract lysine information from pdb file, error: {e}')

        for j, lys_coord in enumerate(lys_coords):

            coords_euc_dists = np.linalg.norm(all_coords - lys_coord, axis=1)
            list_close_points = [a for a, x in enumerate(coords_euc_dists < 15) if x]

            try:
                S = M.get_subset(idxs=list_close_points)
                chain = list_of_chains[j]
                resid = list_of_resid[j]

                #SASA is calculated for that lysine in the small molecule.
                # TODO GW 23.07.26 - update this to be general to other AA as element codes are LYS specific here
                pts_2, indx_2 = S.atomselect(chain, [resid], ["CB", "CG", "CD", "CE", "NZ", 'C03', 'C04', 'C05', 'C06', 'N07'],
                                            use_resname=False, get_index=True)

                x = bb.sasa(S, targets=indx_2, probe=1.4, n_sphere_point=960, threshold=0)
                list_of_sasa.append(x[0])

            except Exception as e:
                print(f'SASA Calculation: 2 - Error obtaining SASA at index value {str(j)} with error: {e}')
                list_of_sasa.append(None)
                report_error_to_file('SASA 2', path, f'Error obtaining SASA at index value {str(j)} with error: {e}', self.error_filename)
                continue

        try:
            df_sasa = pd.DataFrame({'Chain': list_of_chains,
                                'Resid': list_of_resid,
                                'sasa': list_of_sasa})
            if self.include_modified: df_sasa['Modified'] = list_modified

        except Exception as e:
            report_error_to_file('SASA 3', path, str(e), self.error_filename)
            print(f'SASA Calculation: 3 - Failed to write extracted sasa information to dataframe, error: {e}')

        return df_sasa


if __name__ == '__main__':
    sasa = SASA(include_modified=True)
    #print(sasa.calculate_sasa(path=f'1ubq.pdb'))
    #print(sasa.calculate_sasa(path=f'result{os.sep}curated{os.sep}1UBQ-alt-1.pdb'))
    print(sasa.calculate_sasa(path=f'result{os.sep}curated{os.sep}6XZ7-alt1A.pdb'))
    #print(sasa.calculate_sasa(path=f'1nsk_AmberMod0000.pdb'))
