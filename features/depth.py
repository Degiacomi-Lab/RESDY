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
from error_reporting import _report_error_to_file

# Depth packages
try:
    from Bio.PDB import PDBParser
    from Bio.PDB.ResidueDepth import min_dist, get_surface, residue_depth
except Exception as e:
    print(f"biopython and msms unavailable. Unable be able to calculate residue depth. Error: {e}")


class Depth():
    '''
    Class to house the different methods for calculating depth values for structures
    '''

    def calculate_depth(self, path):
        '''
        Calculate the depth of the NZ atom from the surface of the protein within
        the overall protein structure.

        Method
        ------
        Identify all NZ atoms within the protein structure through biobox. Use the PDBparser
        from biopython to calculate the surface of the protein. Loop over the identified positions
        of all of the NZ atoms and use the residue_depth function contained in biopython
        to extract the minimum distances from the surface for each of the lysines.


        Parameters
        ----------
        path : string
            The path of the pdb file that the depth of the lysines are being calculated for.

        Returns
        -------
        df_depth : dataframe
            Dataframe with information on chain, residue number and depth output. Outline:
            Chain   Resid   depth
            x       x       x

        Example
        -------
        >> print(calculate_depth(1ubq.pdb))
        Chain  Resid  depth
        0     A      6   x
        1     A     11   x
        2     A     27   x
        3     A     29   x
        4     A     33   x
        5     A     48   x
        6     A     63   x
        '''
        try:
            M = bb.Molecule(path)
            pos, idx = M.atomselect("*", "*", "NZ", get_index=True)
        except Exception as e:
            _report_error_to_file('Depth 1', path, str(e))
            raise Exception(f">> DEPTH error: could not find NZ atoms within atomic structure - {e}")

        try:
            parser = PDBParser()
            structure = parser.get_structure('structure', path)
            surface = get_surface(structure[0])
        except Exception as e:
            _report_error_to_file('Depth 2', path, str(e))
            raise Exception(f">> DEPTH error: could not get biopython structure - {e}")


        results = []
        for i, coord in enumerate(pos):
            chain = M.data.loc[idx[i], ["chain"]].values[0]
            resid = M.data.loc[idx[i], ["resid"]].values[0]
            mychain = structure[0][chain]
            myres = mychain[int(resid)]
            try:
                #dist = min_dist(pos[i], surface)
                rd = residue_depth(myres, surface)
            except Exception as e:
                _report_error_to_file('Depth 3', path, str(e))
                raise Exception(f">> DEPTH error: failed getting min_dist - {e}")

            results.append([chain, resid, rd])

        df_depth = pd.DataFrame(results, columns=["Chain", "Resid", "depth"])
        return df_depth
