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
from features.error_reporting import report_error_to_file

# Frustration packages
try:
    import frustratometer
except Exception as e:
    print(f"frustratometer unavailable. Unable to calculate frustration. Error: {e}")

pd.set_option('display.max_rows', 200)
class Frustration():
    '''
    Class to house the different methods for calculating frustration metric values
    for structures.
    '''

    def __init__(self, include_modified = False,
                 aa_properties = {'non_modified_codes': ['LYS', 'LYSN'],
                                'modified_codes': ['LYE', 'KCX'],
                                'atom_select_names_nonmod': ['NZ'],
                                'atom_select_names_modified': ['NZ', 'N07']},
                 error_filename = 'measure_errors.txt'):
        '''
        Initialise the Frustration class, include any global variables that are required from
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
            resid chain  Modified  frustration   density
        5       6     A     False    -1.180647  4.474738
        10     11     A     False    -1.277059  1.951710
        26     27     A     False    -0.584062  4.216964
        28     29     A     False    -0.721478  3.747545
        32     33     A     False    -0.706590  2.437673
        47     48     A     False    -0.974550  3.083861
        62     63     A     False    -0.532853  1.358850
        '''
        # temp - TODO modules required to add into readme - openmm, pdbfixer
        # Frustratometer 1 - create structure and AWSEM model
        df_frustration = pd.DataFrame()
        try:
            M = bb.Molecule()
            M.import_pdb(path, include_hetatm=True)
            M_ca = M.get_subset(M.atomselect('*', '*', 'CA', use_resname=True, get_index=True)[1])

            if self.include_modified: idx_res_interest = M_ca.atomselect('*', (self.aa_properties['non_modified_codes'] + self.aa_properties['modified_codes']), 'CA', use_resname=True, get_index=True)[1]
            else: idx_res_interest = M_ca.atomselect('*', self.aa_properties['non_modified_codes'], 'CA', use_resname=True, get_index=True)[1]

            list_modified = list(a in self.aa_properties['modified_codes'] for a in list(M_ca.data['resname']))
            df_frustration = M_ca.data[['resid', 'chain', 'resname']]
            df_frustration = df_frustration.assign(**{'Modified': list_modified})

            out_print_trap = io.StringIO()
            with redirect_stdout(out_print_trap):
                frust_struc = frustratometer.Structure(path)
                model_single_resids = frustratometer.AWSEM(frust_struc, min_sequence_separation_contact=2)
            del out_print_trap
        except Exception as e:
            print(f'Frustratometer calculation 1 - failed to create frustratometer structure or AWSEM model with error: {e}')
            df_frustration['frustration'] = None
            df_frustration['density'] = None
            report_error_to_file('Frustratometer 1', path, str(e), self.error_filename)
            return df_frustration

        # Frustratometer 2 - use model to calculate outputs and sort output dataframe
        try:
            single_residue_awsem_frustration = model_single_resids.frustration(kind='singleresidue')
            resid_densities = model_single_resids.rho_r
            df_frustration = df_frustration.assign(**{'frustration': single_residue_awsem_frustration,
                                                        'density': resid_densities})
            df_frustration_res_interest = df_frustration.iloc[idx_res_interest].drop(columns=['resname'])
            if not self.include_modified: df_frustration_res_interest = df_frustration_res_interest.drop(columns=['Modified'])

            try:
                cleaned_code_to_remove = os.path.splitext(os.path.basename(path))[0] + '_cleaned.pdb'
                os.remove(cleaned_code_to_remove)
            except Exception as ef:
                print(f'Failed to remove cleaned pdb for frustratometer calculation with error {ef}')
        except Exception as e:
            report_error_to_file('Frustratometer 2', path, str(e), self.error_filename)
            print(f'Frustratometer calculation 2 - failed to extract frustratometer outputs or to append data to return dataframe: {e}')

        return df_frustration_res_interest.reset_index(drop=True)


if __name__ == '__main__':
    frust = Frustration(include_modified=True)
    print(frust.calculate_frustration(path=f'data{os.sep}curated{os.sep}1UBQ-alt-1.pdb'))
    #print(frust.calculate_frustration(path=f'data{os.sep}curated{os.sep}1NSK-alt-1.pdb'))
    #print(frust.calculate_frustration(path=f'2I1V-alt-1.pdb'))
