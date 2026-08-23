import re
import os
import sys
import shutil
import subprocess
import pandas as pd
import numpy as np
import biobox as bb
from features.error_reporting import report_error_to_file


class PKA():
    '''
    Class to house the different methods for calculating pKa values for structures
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

        Parameters
        ----------
        outdir : string
            The output directory of the measures calculations. Needed to create the
            propkaoutput directory to store the output files from the PROPKA calculations.
        calc_method : string
            The method to use for caclulating pKa values. Default is set to 'propka'.
            Current options are: propka, pKaANI
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
        self.calc_method = calc_method
        self.include_modified = include_modified
        self.aa_properties = aa_properties
        self.error_filename = error_filename
        if self.error_filename != 'no_record': self.record_errors = True
        else: self.record_errors = False
        if self.calc_method == 'propka':
            self.pka_outdir = os.path.join(outdir, "propkaoutput")
            if not os.path.exists(self.pka_outdir):
                os.makedirs(self.pka_outdir)

    def calculate_pka(self, path):
        '''
        Take the preferred method of calculating the pKa values and call the appropriate
        function from the method options available. Current methods available are:
        - propka: Uses PROPKA3 program
        - pKaANI: Uses the pKaANI program
        
        Parameters
        ----------
        path : string
            The path of the pdb file that pKa is being calculated for
        '''
        match self.calc_method:
            case 'propka':
                return self.calculate_propka(path)
            case 'pkaANI':
                return self.calculate_pkaANI(path)
            case _:
                raise ValueError('Unknown pKa calculation method.' \
                                 'Current options are: propka, pKaANI.')

    def calculate_propka(self, path):
        '''
        Call PROPKA to calculate the pKa of a file, parse the .pka file to extract lysine data
        parse errors, and return a dataframe containing all measurements not yielding an error.
        
        Method
        ------
        Check if propka has been run before on this protein, otherwise run PROPKA3 on the given
        pdb file. Use the function _parse_propka_errors() to identify any lysines within the
        structure that did not caclulate correctly before searching the output file, extracting
        the pka values produced and writing them to df_propka to output.

        Parameters
        ----------
        path : string
            The path of the pdb file that pKa is being calculated for with PROPKA3.

        Returns
        -------
        df_propka : dataframe
            Dataframe with information on chain, residue number and propka output. Outline:
            Chain   Resid   propka
            x       x       x

        Example
        -------
        >> print(calculate_propka(1ubq.pdb))
          Chain  Resid  propka
        0     A      6       x
        1     A     11       x
        2     A     27       x
        3     A     29       x
        4     A     33       x
        5     A     48       x
        6     A     63       x
        '''
        # Get modified residues that aren't titratable by PROPKA
        if self.include_modified:
            df_mod = pd.DataFrame()
            try:
                M = bb.Molecule()
                M.import_pdb(path, include_hetatm=True)
    
                idx_mod_res = M.atomselect('*', self.aa_properties['modified_codes'], 'CA', use_resname=True, get_index=True)[1]
                lys_res_nums = list(M.data['resid'][idx_mod_res])
                list_chains = list(M.data['chain'][idx_mod_res])
                list_modified = list(a in self.aa_properties['modified_codes'] for a in list(M.data['resname'][idx_mod_res]))
                df_mod = pd.DataFrame({'Resid': lys_res_nums,
                                       'Chain': list_chains,
                                       'Modified': list_modified,
                                       'propka': ['NonTitratable'] * len(lys_res_nums)})
            except Exception as e:
                print(f'>> Failed to identify modified residues for path:{path}. '
                      f'If there are any modified residues present, they may be listed as NaN for this structure.')
        code_for_df = os.path.basename(path).split(".")[0]
        propka_error_file_name = os.path.join(self.pka_outdir, f"{code_for_df}_propka_errors.txt")

        propka_exist_output = os.path.join(self.pka_outdir, (code_for_df + '.pka'))
        if not os.path.isfile(propka_exist_output):
            try:
                f = open(propka_error_file_name, 'w')
                process = subprocess.Popen([sys.executable, '-m', 'propka', path],
                                    stdout=f, stderr=f)
                stdout, stderr = process.communicate()
                f.close()
                shutil.move(code_for_df, propka_exist_output)

            except Exception as e:
                f.close()

                try:
                    shutil.move(code_for_df, propka_exist_output)
                except:
                    pass

                if self.record_errors: report_error_to_file('PROPKA 1', path, str(e), self.error_filename)
                raise Exception(f'Failed to obtain pKa data (PROPKA): {e}') from e
        else:
            print(f'>> Previous PROPKA file found for pdb file {path}, using this instead of re-running.')

        try:
            propka_res_interest_fails = self._parse_propka_errors(propka_error_file_name)
        except Exception as e:
            if self.record_errors: report_error_to_file('PROPKA 2', path, str(e), self.error_filename)
            raise Exception(f"Failed extracting PROPKA errors from output file. {e}")

        try:
            propka_outfile = open(propka_exist_output)
        except Exception as e:
            if self.record_errors: report_error_to_file('PROPKA 3', path, str(e), self.error_filename)
            raise Exception(f'Failed to find {propka_file} output file to read') from e

        res_interest_numbers = []
        pka_vals = []
        chain = []
        try:
            for line in propka_outfile:
                # Only searching for non-modified canonical amino acid code as these codes only form part of the titratable subset which can be calculated by PROPKA
                if re.search(rf'^\s{{3}}(?:{self.aa_properties["non_modified_codes"][0]})\s' , line):
                    try:
                        line = line.split()
                        # reject adding entries associated with errors in structure (as per logfile)
                        if len(propka_res_interest_fails)>0:
                            idx = np.where((propka_res_interest_fails["Chain"] == line[2]) & (propka_res_interest_fails["Resid"].astype(int) == int(line[1])))
                            if len(idx[0])>0:
                                print(f'Error associated with lysine {int(line[1])} on chain {line[2]} found in the log file at index {idx[0]} for pKa calculation. Value not taken through to measurements.')
                                continue

                        res_interest_numbers.append(int(line[1]))
                        chain.append(line[2])
                        pka_vals.append(float(line[3]))

                    except Exception as e:
                        report_error_to_file('PROPKA 4', path, str(e), self.error_filename)
                        print(f"> Error {e}")
                        continue

            propka_outfile.close()

        except Exception as e:
            propka_outfile.close()
            if self.record_errors: report_error_to_file('PROPKA 5', path, str(e), self.error_filename)
            raise Exception(f'Failure parsing {code_for_df}.pka, error: {e}') from e

        try:
            df_propka = pd.DataFrame({'Resid': res_interest_numbers,
                                      'Chain': chain,
                                      'propka': pka_vals})

            df_propka = df_propka.dropna()
            df_propka = df_propka.drop_duplicates(subset=None, keep='first', inplace=False, ignore_index=True)
            
            if self.include_modified:
                df_propka = df_propka.assign(**{'Modified': [False] * len(df_propka)})
                if len(df_mod) > 0:
                    df_propka = df_propka.set_index(['Resid', 'Chain']).combine_first(df_mod.set_index(['Resid', 'Chain'])).reset_index()

        except Exception as e:
            if self.record_errors: report_error_to_file('PROPKA 6', path, str(e), self.error_filename)
            raise Exception(f'Failed to construct pKa (PROPKA) dataframe. {e}') from e

        return df_propka.sort_values(by='Resid').reset_index(drop=True)


    def _parse_propka_errors(self, path):
        '''
        parse the PROPKA output file and appends unique chain and resid of any lysines
        mentioned a DataFrame. This list is returned to main and later the residues in it
        are removed from the df.
        
        Method
        ------
        Go over the propka errors output file that is produced when running. Identify the lines
        which contain information about the lysines within the protein structure analysed that
        have errors associated with them. Extract the chain and resid number from this and append
        to a dataframe to return which contains a set of data on the lysines to remove from the
        read output pka values.

        Parameters
        ----------
        path : string
            The path of the errors output file from the PROPKA analysis of the pdb file of interest.

        Returns
        -------
        dataframe
            Dataframe with information on chain, residue number for lysines with calculation errors.
            Outline:
            Chain   Resid
            x       x
        '''
        if not os.path.exists(path):
            return []

        f = open(path, 'r')
        cnt = 0
        df_propka_errors = pd.DataFrame(columns=['Chain', 'Resid'])
        for line in f:
            cnt += 1
            res_interest_raw = re.findall(rf'{self.aa_properties["non_modified_codes"][0]}[\s]+[\d]+[\s]+[\w]+', line)
            if res_interest_raw:
                for line in res_interest_raw:
                    words = line.split()
                    resid = words[1]
                    chain = words[2]
                    if resid != '' and chain != '':
                        df_propka_errors = pd.concat([df_propka_errors, pd.DataFrame([{'Chain': chain, 'Resid': resid}])], axis=0, ignore_index=True)
            else:
                res_interest_raw_2 = re.findall(rf'[\d]*-{self.aa_properties["non_modified_codes"][0]} \(\w\)', line)
                for line in res_interest_raw_2:
                    words = line.split()
                    resid = line.split('-')[0]
                    chain = words[1].strip(' ()')
                    if resid != '' and chain != '':
                        df_propka_errors = pd.concat([df_propka_errors, pd.DataFrame([{'Chain': chain, 'Resid': resid}])], axis=0, ignore_index=True)

        f.close()

        # if the file is completely empty, let's just wipe it!
        if cnt == 0:
            os.remove(path)

        if len(df_propka_errors) == 0:
            return []
        else:
            return df_propka_errors.drop_duplicates()


    def calculate_pkaANI(self, path):
        '''
        A second method for calculating the pKa of the NZ atom of the lysines within the protein
        structure, this time using the external program, pKaANI.

        Method
        ------
        Call a subprocess to open the pkaani program with the desired pdb file given through path.
        Find the log file produced from running this and search this to find the values produced
        for LYS. Translate the data found within the log file to the dataframe to be returned.

        Parameters
        ----------
        path : string
            The path of the pdb file that pKa is being calculated for with pKaANI.

        Returns
        -------
        df_pkaani : dataframe
            Dataframe with information on chain, residue number and pkaani output. Outline:
            Chain   Resid   pkaani
            x       x       x

        Example
        -------
        >> print(calculate_pkaani(1ubq.pdb))
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
    pka = PKA(outdir='result',
              calc_method='propka',
              include_modified=False)
    #print(pka._parse_propka_errors(path=f'data{os.sep}propkaoutput{os.sep}1A0F-alt-1_propka_errors.txt'))
    print(pka.calculate_pka(path=f'result{os.sep}curated{os.sep}1UBQ-alt-1.pdb'))
