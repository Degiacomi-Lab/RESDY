import os
from datetime import datetime
import pandas as pd
import biobox as bb
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
from .error_reporting import report_error_to_file

class SECONDARYSTRUCTURE():
    '''
    Example class for adding your own features into the codebase. To allow the measuring parent
    script to pick it up, please ensure that the class name is the same as the filename.
    '''

    def __init__(self, include_modified=False,
                 numbers_or_letters='letters',
                 aa_properties = {'non_modified_codes': ['LYS', 'LYSN'],
                                  'modified_codes': ['LYE', 'KCX'],
                                  'atom_select_names_nonmod': ['NZ'],
                                  'atom_select_names_modified': ['NZ', 'N07']},
                 error_filename = 'measure_errors.txt'):
        '''
        Initialise the Secondary Structure class.

        :param numbers_or_letters: Option to change the output format between either numbers (0-8
            numbers correspond the the 8 secondary structure classes) or direct letters for these
            8 secondary structure classes
        :type numbers_or_words: str
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
        if self.error_filename != 'no_record':
            self.record_errors = True
        else:
            self.record_errors = False
        if numbers_or_letters not in ['numbers', 'letters']:
            raise ValueError(f'>> Incorrect input given to the option numbers_or_letters: '
                             f'{numbers_or_letters}. Options are \'numbers\' or \'letters\'')
        self.numbers_or_letters = numbers_or_letters

    def calculate(self, path):
        '''
        Calculate the secondary structure for the protein given in path. Option to either report
        the output as numbers of as letters corresponding to the secondary structure. Options:
        - 0 or H
        - 1 or B
        - 2 or E
        - 3 or G
        - 4 or I
        - 5 or T
        - 6 or S
        - 7 or -

        :param path: The path of the pdb file that the feature is being calculated for.
        :type path: str
        :returns:
            Dataframe with information on chain, residue number and feature output. Outline::

                Chain   Resid   feature
                x       x       x
        :rtype: pandas.DataFrame
        '''

        ss_conv_dict = {'H': 0,
                        'G': 1,
                        'I': 2,
                        'P': 3,
                        'E': 4,
                        'B': 5,
                        'T': 6,
                        'S': 7,
                        '-': 8,
                        'C': 8}

        # 1: Load in the structure and locate all the NZ atoms within the lysines, calculate the list of chains and list of resids to go with this
        try:
            #create temporary pdb file with headers to satisfy dssp
            with open(file=path, mode='r') as orig_pdb:
                orig_lines = orig_pdb.readlines()

            if not orig_lines[0].startswith('HEADER'):
                orig_lines.insert(0, f"HEADER    TEMPORARY                               "
                                f"{datetime.today().strftime('%d-%m-%Y')}  {os.path.basename(path)[0]}")

                tmp_file_name = f"{path.split('.')[0]}_tmpdssp.pdb"
                with open(file=tmp_file_name, mode='w') as new_pdb:
                    new_pdb.writelines(orig_lines)
            else:
                tmp_file_name = path

            #run analysis
            M = bb.Molecule()
            M.import_pdb(tmp_file_name, include_hetatm=True)
            M = M.get_subset(idxs=M.atomselect('*', '*', 'CA', use_resname=True, get_index=True)[1])

            p = PDBParser()
            struc = p.get_structure(id=os.path.basename(tmp_file_name), file=tmp_file_name)
            model = struc[0]
            dssp = DSSP(model=model, in_file=tmp_file_name, file_type='PDB')
            sequence = ''
            sec_structure_eight_bit = ''
            for z in range(len(dssp)):
                a_key = list(dssp.keys())[z]
                sequence += dssp[a_key][1]
                sec_structure_eight_bit += dssp[a_key][2]

            df_ss = M.data[['resname', 'chain', 'resid']]

            if self.include_modified:
                idx_atom_interest = M.atomselect('*',
                                    (self.aa_properties['non_modified_codes'] + self.aa_properties['modified_codes']),
                                    'CA',
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
                                      'CA',
                                      use_resname=True, get_index=True)[1]

            if self.numbers_or_letters == 'numbers':
                sec_structure_list = [ss_conv_dict[a] for a in list(sec_structure_eight_bit)]
            else:
                sec_structure_list = list(sec_structure_eight_bit)

            list_modified = list(a in self.aa_properties['modified_codes'] for a in list(M.data['resname']))
            if self.include_modified:
                df_ss = df_ss.assign(**{'Modified': list_modified})

            df_ss = df_ss.assign(**{'secondarystructure': sec_structure_list})
            df_ss = df_ss.iloc[idx_atom_interest]

            os.remove(path=tmp_file_name)

        except Exception as e:
            if os.path.exists(tmp_file_name):
                os.remove(tmp_file_name)
            if self.record_errors: report_error_to_file('SECONDARYSTRUCTURE 1', path, str(e), self.error_filename)
            print(f'SECONDARYSTRUCTURE Calculation: 1 - could not calculate secondary structure for file {path}: {e}')
            return pd.DataFrame(columns=['Chain', 'Resid', 'secondarystructure'])

        return df_ss.rename(columns={'chain': 'Chain', 'resid': 'Resid'}).reset_index(drop=True)


if __name__ == '__main__':
    SS = SECONDARYSTRUCTURE(include_modified=True)
    print(SS.calculate(path=f'demo{os.sep}curated{os.sep}AF-P40616-F1-model_v6.pdb'))
    print(SS.calculate(path=f'demo{os.sep}curated{os.sep}1UPT-alt-1.pdb'))
