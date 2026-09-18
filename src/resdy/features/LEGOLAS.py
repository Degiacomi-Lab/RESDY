import sys
import os
import shutil
import subprocess
import pandas as pd
import biobox as bb
from dotenv import load_dotenv
from .error_reporting import report_error_to_file

try:
    from ase import Atoms
    import torch
    import torchani
except Exception as e:
    print(f'Packages required for legolas calculation are not available, '
          f'will not be able to calculate legolas data. Error: {e}')


class LEGOLAS():
    '''
    15N NMR values for structures.
    '''

    def __init__(self, outdir, include_modified=False,
                 aa_properties = {'non_modified_codes': ['LYS', 'LYSN'],
                                'modified_codes': ['LYE', 'KCX'],
                                'atom_select_names_nonmod': ['NZ'],
                                'atom_select_names_modified': ['NZ', 'N07']},
                 error_filename = 'measure_errors.txt'):
        '''
        Initialise the legolas nmr calculation class

        :param outdir: The output directory that measurements will be saved to.
        :type outdir: str
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

        self.outdir = outdir
        self.include_modified = include_modified
        self.aa_properties = aa_properties
        self.error_filename = error_filename
        if self.error_filename != 'no_record':
            self.record_errors = True
        else:
            self.record_errors = False

        load_dotenv()
        self.legolas_prog = os.getenv('LEGOLAS_PATH')
        if self.legolas_prog is None:
            raise ImportError('Searching for LEGOLAS path returned None, set LEGOLAS_PATH '
                            'or pass legolas_path=, legolas will be removed from features.')
        elif not os.path.isfile(self.legolas_prog):
            raise ImportError(f'legolas.py not found at {self.legolas_prog}, set LEGOLAS_PATH '
                              f'or pass legolas_path=, legolas will be removed from features.')
        else:
            print(f'>> LEGOLAS_PATH read from .env file as: {self.legolas_prog}')

        self.legolas_output_path = os.path.join(self.outdir, 'legolas')
        os.makedirs(self.legolas_output_path, exist_ok=True)


    def calculate_legolas(self, path):
        '''
        Calculate 15N nmr data using legolas

        .. rubric:: Method

        Use biobox to extract the positions of the lysines within the protein structure given in
        the path. Then change directory to the path of legolas and run legolas.py on the desired
        protein structure.

        :param path: The path of the pdb file that legolas is being calculated for.
        :type path: str
        '''
        # 1: Load in the structure and locate all the NZ atoms within the lysines, calculate the list of chains and list of resids to go with this
        try:
            M = bb.Molecule()
            M.import_pdb(path, include_hetatm=True)

            result_filename = os.path.join(self.legolas_output_path, path.split(f'{os.sep}')[-1].split('.')[0] + '_cs.csv')
            if os.path.exists(result_filename):
                already_exists = True
                print(f'Legolas: The file {result_filename} already exists, using this file')
            else:
                already_exists = False
            modified_struc = False
            if not already_exists:
                if 'HIE' in list(M.data['resname']):
                    print('modifying the structure to change HIE to HIS')
                    M.data.loc[M.data['resname'] == 'HIE', 'resname'] = 'HIS'
                    tmp__file_stem = f'tmp_legolas_{os.path.basename(path).split(".")[0]}'
                    M.write_pdb(f'{tmp__file_stem}.pdb')
                    modified_struc = True

            M_n = M.get_subset(M.atomselect('*', '*', 'N', use_resname=True, get_index=True)[1])
            if self.include_modified:
                idx_n_res_interest = M_n.atomselect('*',
                                                    (self.aa_properties['non_modified_codes'] + self.aa_properties['modified_codes']),
                                                    'N', use_resname=True, get_index=True)[1]
            else:
                idx_n_res_interest = M_n.atomselect('*',
                                                    self.aa_properties['non_modified_codes'],
                                                    'N', use_resname=True, get_index=True)[1]

            list_modified = list(a in self.aa_properties['modified_codes'] for a in list(M_n.data['resname']))
            df_legolas = M_n.data[['resid', 'chain', 'resname']]
            if self.include_modified: df_legolas = df_legolas.assign(**{'Modified': list_modified})

        except Exception as e:
            print(f'Legolas: 1 - could not load and identify the NZ atoms within the lysines of the structure: {e}')
            if self.record_errors: report_error_to_file('LEGOLAS 1', path, str(e), self.error_filename)
            return pd.DataFrame(columns=['Chain', 'Resid', 'legolas'])

        # 2: change location to legolas directory and run the legolas program on the specified pdb before changing back to working directory
        try:
            if modified_struc:
                pdb_absolute_path = f'{tmp__file_stem}.pdb'
                result_filename = f'{tmp__file_stem}_cs.csv'
            else:
                pdb_absolute_path = os.path.join(os.getcwd(), path)
                result_filename = path.split(f'{os.sep}')[-1].split('.')[0] + '_cs.csv'
            if not already_exists:
                subprocess.run([sys.executable, self.legolas_prog, pdb_absolute_path, '-atype', 'N'], check=True)
            else:
                result_filename = os.path.join(self.legolas_output_path, result_filename)
            df_nmr = pd.read_csv(result_filename)
            
            legolas_res_key = {'0': 'ALA', '1': 'ARG', '2': 'ASN', '3': 'ASP', '4': 'CYS',
                               '5': 'GLU', '6': 'GLN', '7': 'GLY', '8': 'HIS', '9': 'ILE',
                               '10': 'LEU', '11': 'LYS', '12': 'MET', '13': 'PHE', '14': 'PRO',
                               '15': 'SER', '16': 'THR', '17': 'TRP', '18': 'TYR', '19': 'VAL',
                               '20': 'HOH', '21': 'DOD'}

            if not already_exists:
                os.remove(result_filename.split('.')[0] + '.parquet')
                os.rename(result_filename, path.split(f'{os.sep}')[-1].split('.')[0] + '_cs.csv')
                result_filename = path.split(f'{os.sep}')[-1].split('.')[0] + '_cs.csv'
                shutil.move(result_filename, self.legolas_output_path)
                if modified_struc:
                    os.remove(f'{tmp__file_stem}.pdb')

            df_legolas = df_legolas.assign(**{'_of_interest':[i in set(idx_n_res_interest) for i in range(len(df_legolas))]})

            if len(df_legolas) != len(df_nmr):
                rows_to_drop = []
                for i, r in df_legolas.iterrows():
                    # LEGOLAS only calculates backbone 15N NMR => exactly 1 measurement per residue, any mismatches are gaps in sequence
                    if r['resname'] != legolas_res_key[int(df_nmr['RESIDUE_ID'].iloc[i])]:
                        rows_to_drop.append(i)
                df_legolas = df_legolas.drop(index=rows_to_drop).reset_index(drop=True)
            df_legolas = pd.concat([df_legolas, df_nmr.loc[:, ~df_nmr.columns.str.contains('^Unnamed')]], axis=1)

            df_legolas = df_legolas[df_legolas['_of_interest']].drop(columns=['_of_interest', 'resname', 'ATOM_TYPE', 'RESIDUE_ID', 'CHEMICAL_SHIFT_STD'])
            df_legolas = df_legolas.rename(columns={'CHEMICAL_SHIFT': 'legolas', 'chain': 'Chain', 'resid': 'Resid'})

        except Exception as e:
            print(f'Legolas 2: Failed to run the legolas program and extract the 15N nmr shifts for the protein: {e}')
            if self.record_errors: report_error_to_file('LEGOLAS 2', path, str(e), self.error_filename)
            return pd.DataFrame(columns=['Chain', 'Resid', 'legolas'])

        return df_legolas.reset_index(drop=True)


if __name__ == '__main__':
    nmr = LEGOLAS(outdir='result', include_modified=False)
    print(nmr.calculate_legolas(path=f'result{os.sep}curated{os.sep}1UBQ-alt-1.pdb'))
