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

# Depth packages
try:
    from Bio.PDB import PDBParser
    from Bio.PDB.ResidueDepth import min_dist, get_surface, residue_depth
except Exception as e:
    print(f"biopython and msms unavailable. Unable be able to calculate residue depth. Error: {e}")


class Depth():
    '''
    Class to house the different methods for calculating depth values for structures
    '''

    def __init__(self, include_modified=False, aa_properties = {'non_modified_codes': ['LYS', 'LYSN'],
                                                                'modified_codes': ['LYE', 'KCX'],
                                                                'atom_select_names_nonmod': ['NZ'],
                                                                'atom_select_names_modified': ['NZ', 'N07']},
                 error_filename = 'measure_errors.txt'):
        '''
        Initialise the Depth class, include any global variables that are required from
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

    def calculate_depth(self, path):
        '''
        Calculate the depth of the lysine from the surface of the protein within
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
            if self.include_modified: idx_nz = M.atomselect('*', (self.aa_properties['non_modified_codes'] + self.aa_properties['modified_codes']), self.aa_properties['atom_select_names_modified'], use_resname=True, get_index=True)[1]
            else: idx_nz = M.atomselect('*', self.aa_properties['non_modified_codes'], self.aa_properties['atom_select_names_nomod'], use_resname=True, get_index=True)[1]
            lys_res_nums = list(M.data['resid'][idx_nz])
            list_chains = list(M.data['chain'][idx_nz])
            list_modified = list(a in self.aa_properties['modified_codes'] for a in list(M.data['resname'][idx_nz]))
        except Exception as e:
            report_error_to_file('Depth 1', path, str(e), self.error_filename)
            raise Exception(f">> DEPTH error: could not find NZ atoms within atomic structure - {e}")

        try:
            parser = PDBParser()
            structure = parser.get_structure('structure', path)
            surface = get_surface(structure[0])
        except Exception as e:
            report_error_to_file('Depth 2', path, str(e), self.error_filename)
            raise Exception(f">> DEPTH error: could not get biopython structure - {e}")

        depth_results = []
        for i, idx in enumerate(idx_nz):
            mychain = structure[0][list_chains[i]]
            myres = mychain[int(lys_res_nums[i])]
            try:
                #rd = min_dist(pos[i], surface)  # NZ atom depth
                rd = residue_depth(myres, surface)  # average atom depth for all heavy atoms in residue of interest
            except Exception as e:
                report_error_to_file('Depth 3', path, str(e), self.error_filename)
                raise Exception(f">> DEPTH error: failed getting min_dist - {e}")
            depth_results.append(rd)

        df_depth = pd.DataFrame(columns=["Chain", "Resid", "depth"])
        try:
            df_depth['Chain'] = list_chains
            df_depth['Resid'] = lys_res_nums
            df_depth['depth'] = depth_results
            if self.include_modified: df_depth['Modified'] = list_modified
        except Exception as e:
            report_error_to_file('Depth 4', path, str(e), self.error_filename)
            print(f'Depth Calculation: 4 - Failed to create datafame to append to the overall dataframe: {e}')
        return df_depth

if __name__ == '__main__':
    depth = Depth(include_modified=True)
    print(depth.calculate_depth(path=f'1ubq.pdb'))
    #print(depth.calculate_depth(path=f'tmp_checking_pdb.pdb'))
    #print(depth.calculate_depth(path=f'1nsk_AmberMod0000.pdb'))