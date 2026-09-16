import os
import io
from contextlib import redirect_stdout
import pandas as pd
import biobox as bb
from .error_reporting import report_error_to_file

try:
    import frustratometer
    frustration_packages_available = True
except Exception as e:
    frustration_packages_available = False
    print(f"frustratometer unavailable. Unable to calculate frustration. Error: {e}")

pd.set_option('display.max_rows', 200)


class Frustration():
    '''
    Frustration metric values for structures.
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

        :param include_modified: Toggle to include residues which have been modified within the
            featurisation
        :type include_modified: bool
        :param aa_properties: Properties of the amino acid of interest to investigate modification
            sites for. Defaults to lysine for carbamylation. Properties are the 3 letter codes for
            non modified ('non_modified_codes') and modified ('modified_codes') and the atom names
            for non modified ('atom_select_names_nonmod') and modified
            ('atom_select_names_modified')
        :type aa_properties: dict
        :param error_filename: Name of the text file passed through from overall measures to write
            any errors from calculating features out to.
        :type error_filename: str
        '''
        self.include_modified = include_modified
        self.aa_properties = aa_properties
        self.error_filename = error_filename

        if not frustration_packages_available:
            raise ImportError('>> Packages required for calculating frustation or density (Frustratometer) '
                              'are not available. Frustration or Density will be removed from feature list.')

    def calculate_frustration(self, path):
        '''
        Use the Frustratometer package to identify the frustration metric for the lysines of
        interest.

        .. rubric:: Method

        Load in protein structure into frustratometer before using AWSEM to create a model for this
        with desired parameters. Use this model to calculate the

        :param path: The path of the pdb file that the frustration metric is being calculated for.
        :type path: str
        :returns:
            Dataframe with information on chain, residue number and desired output from
            Frustratometer. Outline for all features::

                Chain   Resid   frustration   density
                x           x             x         x
        :rtype: pandas.DataFrame
        '''
        df_frustration = pd.DataFrame()
        try:
            M = bb.Molecule()
            M.import_pdb(path, include_hetatm=True)
            M_ca = M.get_subset(M.atomselect('*', '*', 'CA', use_resname=True, get_index=True)[1])

            if self.include_modified:
                idx_res_interest = M_ca.atomselect('*',
                                                   (self.aa_properties['non_modified_codes'] + self.aa_properties['modified_codes']),
                                                   'CA', use_resname=True, get_index=True)[1]
            else:
                idx_res_interest = M_ca.atomselect('*',
                                                   self.aa_properties['non_modified_codes'],
                                                   'CA', use_resname=True, get_index=True)[1]

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
            return df_frustration.rename(columns={'chain': 'Chain', 'resid': 'Resid'})

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

        return df_frustration_res_interest.rename(columns={'chain': 'Chain', 'resid': 'Resid'}).reset_index(drop=True)


if __name__ == '__main__':
    frust = Frustration(include_modified=True)
    print(frust.calculate_frustration(path=f'data{os.sep}curated{os.sep}1UBQ-alt-1.pdb'))
    #print(frust.calculate_frustration(path=f'data{os.sep}curated{os.sep}1NSK-alt-1.pdb'))
    #print(frust.calculate_frustration(path=f'2I1V-alt-1.pdb'))
