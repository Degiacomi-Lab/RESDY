import os
import pandas as pd
import biobox as bb
from src.features.error_reporting import report_error_to_file


class Flexibility():
    '''
    Flexibility parameters for lysines within the protein structures.
    '''

    def __init__(self, include_modified, aa_properties = {'non_modified_codes': ['LYS', 'LYSN'],
                                                        'modified_codes': ['LYE', 'KCX'],
                                                        'atom_select_names_nonmod': ['NZ'],
                                                        'atom_select_names_modified': ['NZ', 'N07']},
                 error_filename = 'measure_errors.txt'):
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
        '''
        self.include_modified = include_modified
        self.aa_properties = aa_properties
        self.error_filename = error_filename
        if self.error_filename != 'no_record': self.record_errors = True
        else: self.record_errors = False

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

        .. rubric:: Example

        ::

            >>> print(calculate_flexibility(1ubq.pdb))
              Chain  Resid  flexibility
            0     A      6    10.776667
            1     A     11    15.058889
            2     A     27     7.253333
            3     A     29    13.685556
            4     A     33    20.076667
            5     A     48    13.066667
            6     A     63    15.998889
        '''

        try:
            M = bb.Molecule()
            M.import_pdb(path, include_hetatm=True)

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

        try:
            avg_beta_output = []
            for lys_res, lys_chain in zip(lys_res_nums, list_chains):
                tmp_lys_data = M.data[(M.data['resid'] == lys_res) & (M.data['chain'] == lys_chain)]
                tmp_lys_beta_vals = list(tmp_lys_data['beta'])
                avg_beta = sum(tmp_lys_beta_vals) / len(tmp_lys_beta_vals)
                avg_beta_output.append(avg_beta)
        except Exception as e:
            if self.record_errors: report_error_to_file('Flex 2', path, str(e), self.error_filename)
            print(f'Flex Calculation: 2 - Failed to obtain the beta values and create an average for the lysine of interest at position {lys_res}: {e}')

        try:
            df_flex = pd.DataFrame(columns=['Chain', 'Resid', 'flexibility'])
            df_flex['Chain'] = list_chains
            df_flex['Resid'] = lys_res_nums
            df_flex['flexibility'] = avg_beta_output
            if self.include_modified: df_flex['Modified'] = list_modified
        except Exception as e:
            if self.record_errors: report_error_to_file('Flex 3', path, str(e), self.error_filename)
            print(f'Flex Calculation: 3 - Failed to create dataframe to append to overall measures dataframe: {e}')

        return df_flex


if __name__ == '__main__':
    flex = Flexibility(include_modified=True)
    print(flex.calculate(path=f'result{os.sep}curated{os.sep}1UBQ-alt-1.pdb'))
