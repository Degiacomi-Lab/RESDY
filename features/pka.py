import re
import os
import io
import logging
import datetime
import shutil
import subprocess
import glob
import time
from datetime import date
from multiprocessing import cpu_count
from multiprocessing import Manager
from multiprocessing.pool import Pool
from contextlib import redirect_stdout
from ast import literal_eval
import pandas as pd
import numpy as np
import biobox as bb
import matplotlib.pyplot as plt
#import dill


# AEV packages
try:
    from ase import Atoms
    import torch
    import torchani
except Exception as e:
    print(f'Packages required for AEV calculation are not available, will not be able to calculate AEVs. Error: {e}')

try:
    from Bio.PDB import PDBParser
    from Bio.PDB.ResidueDepth import min_dist, get_surface, residue_depth
except Exception as e:
    print(f"biopython and msms unavailable. Unable be able to calculate residue depth. Error: {e}")

# Frustration packages
try:
    import frustratometer
except Exception as e:
    print(f"frustratometer unavailable. Unable to calculate frustration. Error: {e}")

# Melodia packages
try:
    import melodia_py as mel
except Exception as e:
    print(f"melodia unavailable. Unable to calculate melodia. Error: {e}")

class PKA():
    '''
    Class to house the different methods for calculating pKa values for structures
    '''

    def calculate_pka_propka(self, path):
        '''
        Call PROPKA to calculate the pKa of a file, parse the .pka file to extract lysine data
        parse errors, and return a dataframe containing all measurements not yielding an error.
        
        Method
        ------
        Check if propka has been run before on this protein, otherwise run PROPKA3 on the given
        pdb file. Use the function _parse_propka_errors() to identify any lysines within the structure
        that did not caclulate correctly before searching the output file, extracting the pka
        values produced and writing them to df_propka to output.

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
        0     A      6   x
        1     A     11   x
        2     A     27   x
        3     A     29   x
        4     A     33   x
        5     A     48   x
        6     A     63   x
        '''

        code_for_df = os.path.basename(path).split(".")[0]
        error_file_name = os.path.join(self.pkaoutdir, f"{code_for_df}_propka_errors.txt")

        test_path = code_for_df + '.pka'
        if not os.path.isfile(test_path):
            try:   
                f = open(error_file_name, 'w')
                process = subprocess.Popen(['python', '-m', 'propka', path],
                                    stdout=f, stderr=f)
                stdout, stderr = process.communicate()
                f.close()

            except Exception as e:
                f.close()

                try:
                    shutil.move(code_for_df, os.path.join(self.pkaoutdir, code_for_df))
                except:
                    pass

                if self.report_errors: self._report_error_to_file('PROPKA 1', path, str(e))
                raise Exception(f'Failed to obtain pKa data (PROPKA): {e}') from e

        try:
            propka_lys_fails = self._parse_propka_errors(error_file_name)
        except Exception as e:
            if self.report_errors: self._report_error_to_file('PROPKA 2', path, str(e))
            raise Exception(f"Failed extracting PROPKA errors from output file. {e}")

        try:
            pkafile = code_for_df + '.pka'
            propres = open(pkafile)
        except Exception as e:
            if self.report_errors: self._report_error_to_file('PROPKA 3', path, str(e))
            raise Exception(f'Failed to find {pkafile} output file to read') from e

        lys_number = list()
        pkas = list()
        chain = list()
        try:
            for line in propres:
                if re.search('^   LYS' , line):
                    try:

                        line = line[6:]
                        line = line.split()

                        # reject adding entries associated with errors in structure (as per logfile)
                        if len(propka_lys_fails)>0:
                            idx = np.where((propka_lys_fails["Chain"] == line[1]) & (propka_lys_fails["Resid"].astype(int) == int(line[0])))
                            if len(idx[0])>0:
                                print(f'Error associated with lysine {int(line[0])} on chain {line[1]} found in the log file at index {idx[0]} for pKa calculation. Value not taken through to measurements.')
                                continue

                        lys_number.append(int(line[0]))
                        chain.append(line[1])
                        pkas.append(float(line[2]))

                    except Exception as e:
                        if self.report_errors: self._report_error_to_file('PROPKA 4', path, str(e))
                        print(f"> Error {e}")
                        continue

            propres.close()
            shutil.move(pkafile, os.path.join(self.pkaoutdir, pkafile))

        except Exception as e:
            propres.close()
            shutil.move(pkafile, os.path.join(self.pkaoutdir, pkafile))
            if self.report_errors: self._report_error_to_file('PROPKA 5', path, str(e))
            raise Exception(f'Failure parsing {code_for_df}.pka, error: {e}') from e

        try:
            df_propka = pd.DataFrame({'Resid':lys_number,
                                      'Chain': chain,
                                      'propka':pkas})

            df_propka.sort_values(by=['propka'], inplace=True)
            df_propka = df_propka.dropna()
            df_propka = df_propka.drop_duplicates(subset=None, keep='first', inplace=False, ignore_index=True)

        except Exception as e:
            if self.report_errors: self._report_error_to_file('PROPKA 6', path, str(e))
            raise Exception(f'Failed to construct pKa (PROPKA) dataframe. {e}') from e

        return df_propka


    def _parse_propka_errors(self, path):
        '''
        parse the PROPKA output file and appends unique chain and resid of any lysines mentioned a DataFrame.
        This list is returned to main and later the residues in it are removed from the df.
        
        Method
        ------
        Go over the propka errors output file that is produced when running. Identify the lines which
        contain information about the lysines within the protein structure analysed that have errors
        associated with them. Extract the chain and resid number from this and append to a dataframe
        to return which contains a set of data on the lysines to remove from the read output pka values.

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

        f = open(path, 'r')
        list_remove = list()
        cnt = 0
        for line in f:
            cnt += 1
            lys_raw = re.findall(r'LYS [\d]*[\s][\w]*', line)

            for line in lys_raw:
                words = line.split(' ')
                resid = words[1]
                chain = words[2]
                if resid != "" and chain != "":
                    list_remove.append([chain, resid])

            lys_raw_2 = re.findall(r'[\d]*-LYS \(\w\)', line)

            for line in lys_raw_2:
                words = line.split()
                chain = (words[1])[1:-1]
                words_2 = line.split('-')
                resid = words_2[0]
                if resid != "" and chain != "":
                    list_remove.append([chain, resid])

        f.close()

        # if the file is completely empty, let's just wipe it!
        if cnt == 0:
            os.remove(path)

        if len(list_remove) == 0:
            return []
        else:
            return pd.DataFrame(np.array(list_remove), columns=["Chain", "Resid"]).drop_duplicates()


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
        Chain  Resid  sasa
        0     A      6   x
        1     A     11   x
        2     A     27   x
        3     A     29   x
        4     A     33   x
        5     A     48   x
        6     A     63   x
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
                    if self.report_errors: self._report_error_to_file('pkaANI 1', path, str(e))
                    raise Exception(f'Error: Failed to obtain pkaANI data: {e}. Is pkaANI installed correctly? If so, reload environment and try again.')
                else:
                    if self.report_errors: self._report_error_to_file('pkaANI 1', path, str(e))
                    raise Exception(f"Failed to obtain pkaANI data through running pkaANI, error: {e}")

        try:
            log_file = pdb_path + '_pka.log'
            propres = open(log_file)
        except Exception as e:
            if self.report_errors: self._report_error_to_file('pkaANI 2', path, str(e))
            raise Exception(f'Failed to find pkaANI log file: {e}') from e

        lys_number = []
        pkas = []
        chains = []
        try:
            for line in propres:
                if re.search('^LYS', line):
                    line = line[4:]
                    info = line.split()

                    lys_number.append(int(info[0]))
                    chains.append(info[1])
                    pkas.append(float(info[2]))

            propres.close()
        except Exception as e:
            propres.close()
            if self.report_errors: self._report_error_to_file('pkaANI 3', path, str(e))
            raise Exception(f'Failure parsing pkaANI log file for: {code_for_df}_pka.log. {e}') from e

        try:
            df_pkaani = pd.DataFrame({'Resid':lys_number,
                                      'Chain': chains,
                                      'pkaANI':pkas})

            df_pkaani.sort_values(by=['pkaANI'], inplace=True)
            df_pkaani = df_pkaani.dropna()
            df_pkaani = df_pkaani.drop_duplicates(subset=None, keep='first', inplace=False, ignore_index=True)
        except Exception as e:
            if self.report_errors: self._report_error_to_file('pkaANI 4', path, str(e))
            raise Exception(f'Failed to construct pkaANI dataframe, error: {e}') from e

        return df_pkaani
