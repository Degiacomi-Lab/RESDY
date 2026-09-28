import os
import numpy as np
import pandas as pd
import biobox as bb
from .error_reporting import report_error_to_file


class FLEXIBILITY():
    '''
    Flexibility parameters for lysines within the protein structures.
    '''

    def __init__(self, include_modified=False,
                 aa_properties = {'non_modified_codes': ['LYS', 'LYSN'],
                                  'modified_codes': ['LYE', 'KCX'],
                                  'atom_select_names_nonmod': ['NZ'],
                                  'atom_select_names_modified': ['NZ', 'N07']},
                 error_filename = 'measure_errors.txt',
                 remove_outliers=False):
        '''
        Initialise the Flexibility class, include any global variables that are required from
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
        :param remove_outliers: Toggleable option to remove any beta values which are likely
            outliers from the calculations to avoid bias in normalisation. If set to true, values
            which are outliers will be set to NaN. Default is False and will treat all Beta values
            as valid.
        type remove_outliers: bool
        '''
        self.include_modified = include_modified
        self.aa_properties = aa_properties
        self.error_filename = error_filename
        self.remove_outliers = remove_outliers
        if self.error_filename != 'no_record':
            self.record_errors = True
        else:
            self.record_errors = False

    def calculate(self, path):
        '''
        Calculate the flexibility of the lysine of interest within the protein structure. This uses
        B-factor values for atoms within the lysine and reports the average of all atoms within the
        lysine to form the final scalar quantity.

        .. rubric:: Method

        Read in the column for the B-factor values, take values which correspond to the lysine of
        interest and calculate the average to report.

        :param path: The path of the pdb file that the flexibility is being calculated for.
        :type path: str
        :returns:
            Dataframe with information on chain, residue number and flexibility output.
            Outline::

                Chain   Resid   Flexibility
                x       x       x
        :rtype: pandas.DataFrame
        '''

        try:
            M = bb.Molecule()
            M.import_pdb(path, include_hetatm=True)
            M.atomselect('*', '*', ['C', 'CA', 'N', 'O'], get_index=True)

            if self.include_modified:
                idx_nz = M.atomselect('*',
                                      (self.aa_properties['non_modified_codes'] + self.aa_properties['modified_codes']),
                                      self.aa_properties['atom_select_names_modified'],
                                      use_resname=True, get_index=True)[1]
                # due to wider selection criteria, possible to get more than 1 hit per residue of interest, remove duplicates
                key_res_chain = zip(list(M.data['resid'].values[idx_nz]), list(M.data['chain'].values[idx_nz]))
                pairs_seen, keep_pos = set(), []
                for pair, pos in zip(key_res_chain, range(len(idx_nz))):
                    if pair not in pairs_seen:
                        pairs_seen.add(pair)
                        keep_pos.append(pos)
                idx_nz = idx_nz[keep_pos]

            else:
                idx_nz = M.atomselect('*',
                                      self.aa_properties['non_modified_codes'],
                                      self.aa_properties['atom_select_names_nonmod'],
                                      use_resname=True, get_index=True)[1]

            lys_res_nums = list(M.data['resid'][idx_nz])
            list_chains = list(M.data['chain'][idx_nz])
            list_modified = list(a in self.aa_properties['modified_codes'] for a in list(M.data['resname'][idx_nz]))
        except Exception as e:
            if self.record_errors: report_error_to_file('Flex 1', path, str(e), self.error_filename)
            print(f'Flex Calculation: 1 - Could not load and identify targets within the lysines for calculations: {e}')
            return pd.DataFrame(columns=['Chain', 'Resid', 'Flexibility'])

        af_struc = os.path.basename(path).startswith('AF-')
        if af_struc:
            print(f'>> {path} is an AF structure which means the Beta factor column is the PLDDT '
                f'value, therefore B-factors set to NaN for this file')
            avg_beta_output = [np.nan] * len(lys_res_nums)
        else:
            def _variety_in_resid(c):
                return M.data.groupby(['chain', 'resid'])[c].nunique().max() > 1

            beta_col  ='beta'
            if not _variety_in_resid('beta') and _variety_in_resid('occupancy'):
                print(f'>> Data in the beta column is constant within residues, while occupancy '
                      f'is not, reading beta values from occupancy instead. This is due to '
                      f'correcting for known bug in biobox. File: {path}')
                beta_col = 'occupancy'

            avg_beta_output = []
            mean_beta = M.data[beta_col].mean()
            med_beta = M.data[beta_col].median()
            std_beta = M.data[beta_col].std()

            if self.remove_outliers:
                mad = np.median([np.sqrt((a - med_beta)**2) for a in list(M.data[beta_col])])
                M.data['MAD'] = (0.6745*(M.data[beta_col] - med_beta)) / mad
                mean_beta = M.data.loc[M.data['MAD'] <= 3.5, beta_col].mean()
                std_beta = M.data.loc[M.data['MAD'] <= 3.5, beta_col].std()

            if std_beta == 0 or np.isnan(std_beta):
                print(f'>> B-factor column of file: {path} has no spread on flexibility, '
                      f'all normalised flexibility values set to NaN')
                M.data['normalised_beta'] = np.nan
            else:
                M.data['normalised_beta'] = (M.data[beta_col] - mean_beta) / std_beta

            if self.remove_outliers:
                M.data.loc[M.data['MAD'] > 3.5, 'normalised_beta'] = np.nan

            for lys_res, lys_chain in zip(lys_res_nums, list_chains):
                try:
                    tmp_lys_data = M.data[(M.data['resid'] == lys_res) & (M.data['chain'] == lys_chain)]
                    tmp_lys_beta_vals = list(tmp_lys_data['normalised_beta'])
                    avg_beta = sum(tmp_lys_beta_vals) / len(tmp_lys_beta_vals)
                    avg_beta_output.append(avg_beta)
                except Exception as e:
                    if self.record_errors: report_error_to_file('Flex 2', path, str(e), self.error_filename)
                    print(f'Flex Calculation: 2 - Failed to obtain the beta value for residue at position '
                          f'{lys_res}, set to NaN; error: {str(e)}')
                    avg_beta_output.append(np.nan)

        if all(b == 0 for b in avg_beta_output):
            print(f'>> All Beta factors are 0 in {path}, flexibility set to NaN for all residues of interest')
            avg_beta_output = [np.nan] * len(avg_beta_output)

        try:
            df_flex = pd.DataFrame(columns=['Chain', 'Resid', 'flexibility'])
            df_flex['Chain'] = list_chains
            df_flex['Resid'] = lys_res_nums
            df_flex['flexibility'] = avg_beta_output
            if self.include_modified: df_flex['Modified'] = list_modified
        except Exception as e:
            if self.record_errors: report_error_to_file('Flex 3', path, str(e), self.error_filename)
            print(f'Flex Calculation: 3 - Failed to create dataframe to append to '
                  f'overall measures dataframe: {e}')
            return pd.DataFrame(columns=['Chain', 'Resid', 'flexibility'])

        return df_flex


if __name__ == '__main__':
    flex = FLEXIBILITY(include_modified=True, remove_outliers=True)
    print(flex.calculate(path=f'demo{os.sep}curated{os.sep}1A6M-alt1A.pdb'))
