import os
import pandas as pd
from .error_reporting import report_error_to_file


class STRUCTURE():
    '''
    Structural measurement values for structures.
    '''

    def __init__(self, melodia_features=['all'],
                 include_modified=False,
                 aa_properties = {'non_modified_codes': ['LYS', 'LYSN'],
                                'modified_codes': ['LYE', 'KCX'],
                                'atom_select_names_nonmod': ['NZ'],
                                'atom_select_names_modified': ['NZ', 'N07']},
                 error_filename = 'measure_errors.txt'):
        '''
        Initialise the Structure class using melodia-py.

        :param melodia_features: List of features which are calculated through melodia which has
            been requested when the Measure class is initialised. Default is set to ['all'].
        :type melodia_features: list
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

        try:
            import melodia_py as mel
            self.mel = mel
        except ImportError as e:
            raise ImportError(f'>> Packages required for melodia calculations (melodia_py) are '
                            f'not available, melodia will be removed from features to calculate. '
                            f'Error: {str(e)}') from e

        self.melodia_features = melodia_features
        self.include_modified = include_modified
        self.aa_properties = aa_properties
        self.error_filename = error_filename
        if self.error_filename != 'no_record':
            self.record_errors = True
        else:
            self.record_errors = False
        if self.melodia_features == ['all']:
            self.melodia_features = ['curvature', 'writhing', 'torsion', 'arc_length', 'phi', 'psi']


    def calculate(self, path):
        '''
        Call the Melodia package to calculate data for the following structural features of the
        lysines of interest within the structure: curvature, arc-length, phi, psi

        .. rubric:: Method

        Call melodia on the path of the pdb file that has been passed to the function and create
        dataframe from the results. Process the dataframe to remove any of the calculated features
        that were not asked for.

        :param path: The path of the pdb file that the structural features are being calculated for.
        :type path: str
        :returns:
            Dataframe with information on chain, residue number and desired melodia output.
            Outline for all features::

                Chain   Resid   curvature   writhing    torsion   arc-length  phi psi
                x           x           x          x          x            x    x   x
        :rtype: pandas.DataFrame
        '''
        # Melodia 1 - Calculating geometry using melodia-py
        try:
            melodia_results = self.mel.geometry_from_structure_file(path)
            if isinstance(melodia_results, pd.Series):
                melodia_results = melodia_results.to_frame().T

        except Exception as e:
            if self.record_errors: report_error_to_file('Melodia 1', path, str(e), self.error_filename)
            print(f'Melodia 1: Error processing input file - {path} with error: {e}')
            return pd.DataFrame(columns=['Chain', 'Resid', 'curvature'])

        # Melodia 2 - Formatting and filtering
        try:
            melodia_results.rename({"chain": "Chain", "order": "Resid"}, axis="columns", inplace = True)
            lys_results = melodia_results['name'].isin((self.aa_properties['non_modified_codes'] +
                                                        self.aa_properties['modified_codes']))
            df_melodia = melodia_results[lys_results].copy()
            df_melodia.reset_index(inplace=True, drop=True)
            list_modified = list(a in self.aa_properties['modified_codes'] for a in list(df_melodia['name']))
            cols_to_drop = ['code', 'id', 'model', 'curvature', 'writhing',
                            'torsion', 'phi', 'psi', 'name', 'arc_length']
            cols_to_drop = [col for col in cols_to_drop if col not in self.melodia_features]
            df_melodia.drop(labels=cols_to_drop, axis = 'columns', inplace=True)
            if self.include_modified:
                df_melodia['Modified'] = list_modified

        except Exception as e:
            if self.record_errors:
                report_error_to_file('Melodia 2', path, str(e), self.error_filename)
            print(f'Melodia 2: Unable to reformat melodia output correctly for '
                  f'input {path} with error: {e}')
            return pd.DataFrame(columns=['Chain', 'Resid', 'curvature'])

        return df_melodia


if __name__ == '__main__':
    struc = STRUCTURE(include_modified=True)
    print(struc.calculate(path=f'result{os.sep}curated{os.sep}1UBQ-alt-1.pdb'))
