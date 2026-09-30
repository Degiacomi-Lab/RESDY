import re
import os
import sys
import shutil
import subprocess
import pandas as pd
import numpy as np
import biobox as bb
from .error_reporting import report_error_to_file


class PROPKA():
    '''
    pKa values for structures, calculated with PROPKA3.
    '''

    #: Residues PROPKA3 reports a pKa for. :class:`Measure <resdy.measure.Measure>` drops this
    #: feature when the residue of interest is not one of them, since it would otherwise
    #: contribute a column of NaN.
    SUPPORTED_RESIDUES = {'ASP', 'GLU', 'HIS', 'CYS', 'TYR', 'LYS', 'ARG'}

    #: Why the feature is unavailable, quoted back to the user when it is dropped.
    UNSUPPORTED_REASON = ('PROPKA3 reports a pKa only for ASP, GLU, HIS, CYS, TYR, LYS '
                          'and ARG')

    def __init__(self, outdir,
                 include_modified=False,
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
        if self.error_filename != 'no_record':
            self.record_errors = True
        else:
            self.record_errors = False
        self.pka_outdir = os.path.join(outdir, "propkaoutput")
        if not os.path.exists(self.pka_outdir):
            os.makedirs(self.pka_outdir)


    def calculate(self, path):
        '''
        Call PROPKA to calculate the pKa of a file, parse the .pka file to extract lysine data parse
        errors, and return a dataframe containing all measurements not yielding an error.

        .. rubric:: Method

        Check if propka has been run before on this protein, otherwise run PROPKA3 on the given pdb
        file. Use the function _parse_propka_errors() to identify any lysines within the structure
        that did not calculate correctly before searching the output file, extracting the pka values
        produced and writing them to df_propka to output.

        :param path: The path of the pdb file that pKa is being calculated for with PROPKA3.
        :type path: str
        :returns:
            Dataframe with information on chain, residue number and propka output. Outline::

                Chain   Resid   propka
                x       x       x
        :rtype: pandas.DataFrame
        '''
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

        propka_write_output = f'{code_for_df}.pka'
        propka_exist_output = os.path.join(self.pka_outdir, (code_for_df + '.pka'))

        reuse_propka = (os.path.isfile(propka_exist_output) and
                        os.path.getmtime(propka_exist_output) >= os.path.getmtime(path))

        if not reuse_propka:
            if os.path.isfile(propka_exist_output):
                print(f'>> PROPKA output for {code_for_df} is older than the curated structure '
                      f'for {path}, recalculating propka values')

            f = None
            try:
                f = open(propka_error_file_name, 'w')
                process = subprocess.Popen([sys.executable, '-m', 'propka', path],
                                    stdout=f, stderr=f)
                stdout, stderr = process.communicate()
                f.close()
                shutil.move(propka_write_output, propka_exist_output)

            except Exception as e:
                if f is not None and not f.closed: f.close()

                if os.path.isfile(propka_write_output):
                    try:
                        shutil.move(propka_write_output, propka_exist_output)
                    except Exception:
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
            raise Exception(f'Failed to find {propka_exist_output} output file to read') from e

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
        Parse the PROPKA output file and append the unique chain and resid of any lysines mentioned in
        it to a DataFrame. This list is returned to main and later the residues in it are removed
        from the df.

        .. rubric:: Method

        Go over the propka errors output file that is produced when running. Identify the lines
        which contain information about the lysines within the protein structure analysed that have
        errors associated with them. Extract the chain and resid number from this and append to a
        dataframe to return which contains a set of data on the lysines to remove from the read
        output pka values.

        :param path: The path of the errors output file from the PROPKA analysis of the pdb file of
            interest.
        :type path: str
        :returns:
            Dataframe with information on chain, residue number for lysines with calculation
            errors. Outline::

                Chain   Resid
                x       x
        :rtype: pandas.DataFrame
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


if __name__ == '__main__':
    P = PROPKA(outdir='result',
              include_modified=False)
    #print(pka._parse_propka_errors(path=f'data{os.sep}propkaoutput{os.sep}1A0F-alt-1_propka_errors.txt'))
    print(P.calculate(path=f'result{os.sep}curated{os.sep}1UBQ-alt-1.pdb'))
