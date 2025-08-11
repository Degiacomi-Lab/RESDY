import os
import io
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

# Frustration packages
try:
    import frustratometer
except Exception as e:
    print(f"frustratometer unavailable. Unable to calculate frustration. Error: {e}")


class Frustration():
    '''
    Class to house the different methods for calculating frustration metric values
    for structures.
    '''

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
            _report_error_to_file('Frustratometer 1', path, str(e))
            return df_frustration

        # Frustratometer 2 - use model to calculate outputs
        try:
            single_residue_awsem_frustration = model_single_resids.frustration(kind='singleresidue')
            resid_densities = model_single_resids.rho_r
        except Exception as e:
            _report_error_to_file('Frustratometer 2', path, str(e))
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
            _report_error_to_file('Frustratometer 3', path, str(e))
            print(f'Frustratometer calculation 3 - failed to append data to return dataframe: {e}')

        return df_frustration
