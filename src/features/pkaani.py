import re
import os
import subprocess
import pandas as pd
from src.features.error_reporting import report_error_to_file


class PKAANI():
    '''
    pKa values for structures, calculated with pKaANI.
    '''

    def __init__(self, outdir, calc_method='propka', include_modified=False,
                 aa_properties = {'non_modified_codes': ['LYS', 'LYSN'],
                                'modified_codes': ['LYE', 'KCX'],
                                'atom_select_names_nonmod': ['NZ'],
                                'atom_select_names_modified': ['NZ', 'N07']},
                 error_filename = 'measure_errors.txt'):
        '''
        Setup the PKA class as required. Take input on which method to use and if to include
        modified structures.

        :param outdir: The output directory of the measures calculations. Needed to create the
            propkaoutput directory to store the output files from the PROPKA calculations.
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
        self.include_modified = include_modified
        self.aa_properties = aa_properties
        self.error_filename = error_filename
        if self.error_filename != 'no_record': self.record_errors = True
        else: self.record_errors = False


    def calculate(self, path):
        '''
        A second method for calculating the pKa of the NZ atom of the lysines within the protein
        structure, this time using the external program, pKaANI.

        .. rubric:: Method

        Call a subprocess to open the pkaani program with the desired pdb file given through path.
        Find the log file produced from running this and search this to find the values produced for
        LYS. Translate the data found within the log file to the dataframe to be returned.

        :param path: The path of the pdb file that pKa is being calculated for with pKaANI.
        :type path: str
        :returns:
            Dataframe with information on chain, residue number and pkaani output. Outline::

                Chain   Resid   pkaani
                x       x       x
        :rtype: pandas.DataFrame

        .. rubric:: Example

        ::

            >>> print(calculate_pkaani(1ubq.pdb))
              Chain  Resid pkaani
            0     A      6      x
            1     A     11      x
            2     A     27      x
            3     A     29      x
            4     A     33      x
            5     A     48      x
            6     A     63      x
        '''
        code_for_df = os.path.basename(path).split(".")[0]
        pdb_path = path.split(".")[0]
        test_path = pdb_path + '_pka.log'
        if not os.path.isfile(test_path):
            try:
                _ = subprocess.run(['pkaani', '-i', path])
            except Exception as e:
                print(e)
                if "[Errno 2] No such file or directory: 'pkaani'" == str(e):
                    if self.record_errors: report_error_to_file('pkaANI 1a', path, str(e), self.error_filename)
                    raise Exception(f'Error: Failed to obtain pkaANI data: {e}. Is pkaANI installed correctly? If so, reload environment and try again.')
                else:
                    if self.record_errors: report_error_to_file('pkaANI 1b', path, str(e), self.error_filename)
                    raise Exception(f"Failed to obtain pkaANI data through running pkaANI, error: {e}")

        try:
            log_file = pdb_path + '_pka.log'
            propres = open(log_file)
        except Exception as e:
            if self.record_errors: report_error_to_file('pkaANI 2', path, str(e), self.error_filename)
            raise Exception(f'Failed to find pkaANI log file: {e}') from e

        res_interest_resids = []
        pkas = []
        chains = []
        try:
            for line in propres:
                if re.search(f'^{self.aa_properties["non_modified_codes"][0]}', line):
                    line = line[4:]
                    info = line.split()

                    res_interest_resids.append(int(info[0]))
                    chains.append(info[1])
                    pkas.append(float(info[2]))

            propres.close()
        except Exception as e:
            propres.close()
            if self.record_errors: report_error_to_file('pkaANI 3', path, str(e), self.error_filename)
            raise Exception(f'Failure parsing pkaANI log file for: {code_for_df}_pka.log. {e}') from e

        try:
            df_pkaani = pd.DataFrame({'Resid':res_interest_resids,
                                      'Chain': chains,
                                      'pkaANI':pkas})

            df_pkaani.sort_values(by=['pkaANI'], inplace=True)
            df_pkaani = df_pkaani.dropna()
            df_pkaani = df_pkaani.drop_duplicates(subset=None, keep='first', inplace=False, ignore_index=True)
        except Exception as e:
            if self.record_errors: report_error_to_file('pkaANI 4', path, str(e), self.error_filename)
            raise Exception(f'Failed to construct pkaANI dataframe, error: {e}') from e

        return df_pkaani


if __name__ == '__main__':
    P = PKAANI(outdir='result',
              include_modified=False)
    #print(P._parse_propka_errors(path=f'data{os.sep}propkaoutput{os.sep}1A0F-alt-1_propka_errors.txt'))
    print(P.calculate(path=f'result{os.sep}curated{os.sep}1UBQ-alt-1.pdb'))
