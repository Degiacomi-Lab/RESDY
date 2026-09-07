import os
import pandas as pd
import numpy as np
import biobox as bb
from src.features.error_reporting import report_error_to_file

class Feature():
    '''
    Example class for adding your own features into the codebase. To allow the measuring parent script
    to pick it up, please enure that the class name if the same as the filename.
    '''

    def __init__(self, include_modified=False,
                 aa_properties = {'non_modified_codes': ['LYS', 'LYSN'],
                                  'modified_codes': ['LYE', 'KCX'],
                                  'atom_select_names_nonmod': ['NZ'],
                                  'atom_select_names_modified': ['NZ', 'N07']},
                 error_filename = 'measure_errors.txt'):
        '''
        Initialise the feauture class here, note custom variables will not be taken through
        to the measurements unless defined as default

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
        if self.error_filename != 'no_record': self.record_errors = True
        else: self.record_errors = False

    def calculate(self, path):
        '''
        Calculate the feature

        Parameters
        ----------
        path : string
            The path of the pdb file that DAS is being calculated for.

        Returns
        -------
        df_feature : dataframe
            Dataframe with information on chain, residue number and feature output. Outline:
            Chain   Resid   das
            x       x       x
        '''

        # 1: Load in the structure and locate all the NZ atoms within the lysines, calculate the list of chains and list of resids to go with this
        try:
            M = bb.Molecule()
            M.import_pdb(path, include_hetatm=True)

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

            res_interest_nums = list(M.data['resid'][idx_atom_interest])
            list_chains = list(M.data['chain'][idx_atom_interest])
            list_modified = list(a in  self.aa_properties['modified_codes'] for a in list(M.data['resname'][idx_atom_interest]))
        except Exception as e:
            if self.record_errors: report_error_to_file('DAS 1', path, str(e), self.error_filename)
            print(f'DAS Calculation: 1 - could not load and identify the atoms of interest in the target residues of the structure: {e}')

        # 2: Calculate the feature values for each of the residues of interest identified
        feature_output = []

        # 3: Create dataframe to return
        df_feature = pd.DataFrame(columns=["Chain", "Resid", "das"])
        try:
            df_feature['Chain'] = list_chains
            df_feature['Resid'] = res_interest_nums
            df_feature['das'] = feature_output
            if self.include_modified: df_feature['Modified'] = list_modified
        except Exception as e:
            if self.record_errors: report_error_to_file('DAS 3', path, str(e), self.error_filename)
            print(f'DAS Calculation: 3 - Failed to create datafame to append to the overall measurements dataframe: {e}')

        return df_feature


if __name__ == '__main__':
    feat = Feature(include_modified=True)
    print(feat.calculate(path=f'result{os.sep}curated{os.sep}1UBQ-alt-1.pdb'))
