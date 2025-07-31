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

class Charge():
    '''
    Class to house the different methods for calculating charge values for structures
    '''
    
    def calculate_seqcharge(self, path, num_add_aa=10):
        '''
        Calculate the Sequence Charge of the local sequence around a LYS of interest.
        This is a single value number representing the summation of the charges of the amino acids
        over the specified number of amino acids either side of the lysine.

        Method
        ------
        Take the PDB file and extract the overall sequence using BioBox. Identify all lysines
        within the structure and any shift that has taken place in the PDB file compared to
        the Uniprot sequence. Extract sequences for the lysine of interest and calculate a
        value for the charge based on the summation of charged residues within the sequence.

        Parameters
        ----------
        path : string
            The path of the pdb file that DAS is being calculated for.

        num_add_aa : int
            The number of amino acids to include either side of the lysine of interest.
            An optional parameter which is set to 10 by default.

        Returns
        -------
        df_seqcharge : dataframe
            Dataframe with information on chain, residue number and seqcharge output. Outline:
            Chain   Resid   seqcharge
            x           x           x

        Example
        -------
        >> print(self.calculate_seqcharge(1M2F-alt-1.pdb))
        Chain   Resid  seqcharge
        0     A      95         -3
        '''

        # 1: Extract the overall sequence for the protein given
        try:
            M = bb.Molecule(path)
            idx_nz = M.atomselect('*', 'LYS', 'NZ', use_resname=True, get_index=True)[1]
            lys_res_nums = list(M.data['resid'][idx_nz])
            list_chains = list(M.data['chain'][idx_nz])

            c_alpha_idxs = M.atomselect('*', '*', 'CA', use_resname=True, get_index=True)[1]

            # TODO GW 16.04.25 - eventually will need to add in ability to use  letter codes and charges
            #                    for the 3 lettter cases rather than the 1 letter cases which when using
            #                    modified residues may run into problems

            protein_letters_dict = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
                                    'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G',
                                    'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
                                    'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
                                    'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
                                    'HIE': 'H', 'HID': 'H', 'HIP': 'H', 'LYN': 'K',
                                    'ASX': 'B', 'GLX': 'Z', 'SEC': 'U', 'PYL': 'O',
                                    'XAA': 'X', 'XLE': 'J', 'PSER': 'p', 'PTHR': 't',
                                    'PTYR': 'y', 'MELYS': 'k', 'MEARG': 'r', 'ACLYS': 'k',
                                    'KCX': 'X', 'LYE': 'X'}  
            # KCX and LYE down as X so that they are not treated as positive K

            def _catch(func, *args, handle=lambda e : e, **kwargs):
                try:
                    return func(*args, **kwargs)
                except Exception as e:
                    print(f"Could not convert {e} to a 1 letter code: using 'X' instead")
                    return 'X'

            sequence = ''.join([_catch(lambda : protein_letters_dict[a.upper()]) for a in list(M.data['resname'][c_alpha_idxs])])


        except Exception as e:
            if self.report_errors: self._report_error_to_file('Seqcharge 1', path, str(e))
            print(f'SeqCharge Calculation: 1 - could not extract the sequence from the protein file given: {e}')
            return pd.DataFrame(columns=["Chain", "Resid", "seqcharge"])

        # 2: Extract local sequences based on the overall chain, calculate charge score and add to output
        seqcharge_output = []
        for i, lys_res_idx in enumerate(lys_res_nums):
            try:
                # get the index of the start position of the residues to work out the shift
                shift_val = int(M.data['resid'].iloc[0]) - 1
                seq_lys_index = lys_res_idx - shift_val - 1

                start_idx = seq_lys_index - num_add_aa
                end_idx = seq_lys_index + num_add_aa + 1
                start_null = 0
                end_null = 0

                if start_idx < 0:
                    start_null = - start_idx
                    start_idx = 0

                if end_idx >= len(sequence):
                    end_null = - end_idx
                    end_idx = len(sequence)

                seq = ('-' * start_null) + sequence[start_idx:end_idx] + ('-' * end_null)
                seq_split = list(seq)
                if seq_split[10] != 'K':
                    print(f'A lysine was not found at the desired position {lys_res_idx+1} read in for PDB file {path}; sequence -> {seq}')
                    continue

                pos_aa = ['K', 'H', 'R']
                neg_aa = ['D', 'E']
                count = 0
                for aa in seq:
                    if aa in pos_aa:
                        count += 1
                    elif aa in neg_aa:
                        count -= 1
                seqcharge_output.append(count)

            except Exception as e:
                if self.report_errors: self._report_error_to_file('Seqcharge 2', path, str(e))
                print(f'SeqCharge Calculation 2: Could not calculate a charge for lysine at position {lys_res_idx}, error: {e}')
                seqcharge_output.append(None)

        # 3: Create dataframe to return
        df_seqcharge = pd.DataFrame(columns=["Chain", "Resid", "seqcharge"])
        try:
            df_seqcharge['Chain'] = list_chains
            df_seqcharge['Resid'] = lys_res_nums
            df_seqcharge['seqcharge'] = seqcharge_output
        except Exception as e:
            if self.report_errors: self._report_error_to_file('Seqcharge 3', path, str(e))
            print(f'SeqCharge Calculation: 3 - Failed to create datafame to append to the overall dataframe: {e}')

        #print(df_seqcharge)
        return df_seqcharge