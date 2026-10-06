import os
import pandas as pd
import numpy as np
import biobox as bb
from .error_reporting import report_error_to_file

class FEATURE():
    '''
    Example class for adding your own features into the codebase. To allow the measuring parent
    script to pick it up, please ensure that the class name is the same as the filename.
    '''

    def __init__(self, include_modified=False,
                 aa_properties = {'non_modified_codes': ['LYS', 'LYSN'],
                                  'modified_codes': ['LYE', 'KCX'],
                                  'atom_select_names_nonmod': ['NZ'],
                                  'atom_select_names_modified': ['NZ', 'N07']},
                 error_filename = 'measure_errors.txt'):
        '''
        Initialise the feature class here, note custom variables will not be taken through to the
        measurements unless defined as default

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
        Calculate the feature

        :param path: The path of the pdb file that the feature is being calculated for.
        :type path: str
        :returns:
            Dataframe with information on chain, residue number and feature output. Outline::

                Chain   Resid   feature
                x       x       x
        :rtype: pandas.DataFrame
        '''

        try:
            M = bb.Molecule()
            M.import_pdb(path, include_hetatm=True)
            A = M.get_subset(indices=M.atomselect('*', '*', 'CA', use_resname=True, get_index=True)[1])
            df_feature = A.data[['resname', 'chain', 'resid']]

            if self.include_modified:
                idx_atom_interest = M.atomselect('*',
                                      (self.aa_properties['non_modified_codes'] + self.aa_properties['modified_codes']),
                                      self.aa_properties['atom_select_names_modified'],
                                      use_resname=True, get_index=True)[1]
                # due to wider selection criteria, possible to get more than 1 hit per residue of interest, remove duplicates
                key_res_chain = zip(list(M.data['resid'].values[idx_atom_interest]), list(M.data['chain'].values[idx_atom_interest]))
                pairs_seen, keep_pos = set(), []
                for pair, pos in zip(key_res_chain, range(len(idx_atom_interest))):
                    if pair not in pairs_seen:
                        pairs_seen.add(pair)
                        keep_pos.append(pos)
                idx_atom_interest = idx_atom_interest[keep_pos]

            else:
                idx_atom_interest = M.atomselect('*',
                                      self.aa_properties['non_modified_codes'],
                                      self.aa_properties['atom_select_names_nonmod'],
                                      use_resname=True, get_index=True)[1]

            list_modified = list(a in self.aa_properties['modified_codes'] for a in list(M.data['resname']))
            if self.include_modified:
                df_feature = df_feature.assign(**{'Modified': list_modified})

            # sort calculating data and assign to dataframe before cutting down to size
            feature_list = []
            df_feature = df_feature.assign(**{'feature': feature_list})
            df_feature = df_feature.iloc[idx_atom_interest]


        except Exception as e:
            if self.record_errors: report_error_to_file('FEATURE 1', path, str(e), self.error_filename)
            print(f'FEATURE Calculation: 1 - could not load and identify the atoms of interest in the target residues of the structure: {e}')

        return df_feature.rename(columns={'chain': 'Chain', 'resid': 'Resid'}).reset_index(drop=True)


if __name__ == '__main__':
    feat = FEATURE(include_modified=True)
    print(feat.calculate(path=f'result{os.sep}curated{os.sep}1UBQ-alt-1.pdb'))
