import re
import os, sys
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
import torch
import matplotlib.pyplot as plt
from collections import OrderedDict
#from features.error_reporting import report_error_to_file

try: 
    import esm
except Exception as e:
    print(f'>> Failed to import esm package required for esm calculations, will not be able to calculate sequence features based on esm')


class Ensemble():
    '''
    Class to house the several metrics for calculating ensemble based features for a
    set of given structures
    '''

    def __init__(self, df_proteins, include_modified = False):
        '''
        Initialise the Charge class, include any global variables that are required from
        measures in here.

        Parameters
        ----------
        df_prot : dataframe
            Dataframe including the information passed into measures from protein about which
            structures correspond to the uniprot codes. 
        include_modified : bool
            Toggle to include residues which have been modified within the featurisation
        '''
        self.df_proteins = df_proteins
        self.include_modified = include_modified
    
    def calculate_rmsf(self, path, df_subset):
        '''
        Calculate the root mean square fluctuation of the lysines within the protein
        over all the structures which have been curated for the uniprot code

        Method
        ------
        Find all structures, loop over structures aligning and calculating deviations,
        calculate RMSF and return values as dataframe

        Parameters
        ----------
        path : string
            The path of the pdb file that DAS is being calculated for.
        df_subset : Dataframe
            Subset of the protein dataframe which is given as input to the measure class
            which shows which pdb files are associated with the uniprot code of interest

        Returns
        -------
        df_rmsf : dataframe
            Dataframe with information on chain, residue number and seqcharge output. Outline:
            Chain   Resid   rmsf
            x           x      x

        Example
        -------
        >> print(self.calculate_rmsf(1M2F-alt-1.pdb))
          Chain   Resid     rmsf
        0     A      95       -3
        '''
        # check over measures to see if this has already been calculated as can just copy values due to being the same calculation each time
        file_header = f'{os.sep}'.join(path.split('/')[:-1])
        df_prot_info = df_subset[df_subset['PDB_Code'] == file_header]
        # search over the protein dataframe to extract subset dataframe
        df_tmp_uniprot = df_subset[df_subset['Uniprot_Entry']]

        #https://userguide.mdanalysis.org/stable/examples/analysis/alignment_and_rms/rmsf.html


    def calculate_esm(self, path, num_add_aa=20):
        '''
        Calculate the ESM LLM output for subset protein sequences obtained through searching
        over the path of the protein structure given.

        Parameters
        ----------
        path -> string
            Path of the pdb file of interest
        num_add_aa -> int
            Number of amino acids to go either side of the lysine of interest
        
        Example
        -------
        >> print(E.calculate_esm('1UBQ-alt-1.pdb', num_add_aa = 20))
          Chain Resid                                                esm
        0     A     6  [-0.027814002707600594, -0.03618992120027542, ...
        1     A    11  [-0.0581485778093338, 0.01636454090476036, -0....
        2     A    27  [0.06258092075586319, -0.024953410029411316, -...
        3     A    29  [0.010486193001270294, 0.03880579397082329, -0...
        4     A    33  [-0.016549628227949142, 0.061681147664785385, ...
        5     A    48  [-0.040656737983226776, -0.01545296423137188, ...
        6     A    63  [-0.05812466889619827, 0.03620311990380287, 0....
        '''
        # 1: Extract the sequence sections from the given protein structure
        try:
            M = bb.Molecule(path)
            if self.include_modified: idx_nz = M.atomselect('*', ['LYS', 'LYE', 'KCX'], ['NZ', 'N07'], use_resname=True, get_index=True)[1]
            else: idx_nz = M.atomselect('*', 'LYS', 'NZ', use_resname=True, get_index=True)[1]
            lys_res_nums = list(M.data['resid'][idx_nz])
            list_chains = list(M.data['chain'][idx_nz])
            list_modified = list(a in ['KCX', 'LYE'] for a in list(M.data['resname'][idx_nz]))

            c_alpha_idxs = M.atomselect('*', '*', 'CA', use_resname=True, get_index=True)[1]
            subset_data = M.data.iloc[c_alpha_idxs]

            # for purposes of esm - needs to be asigned to original 20 AA codes
            # similar residues listed as normal equivalent, anything else labelled as X
            # (may need <unk> tag for esm instead)
            protein_letters_dict = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
                                    'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G',
                                    'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
                                    'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
                                    'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
                                    'HIE': 'H', 'HID': 'H', 'HIP': 'H', 'LYN': 'K',
                                    'ASX': 'B', 'GLX': 'Z', 'SEC': 'U', 'PYL': 'O',
                                    'XAA': 'X', 'XLE': 'J', 'PSER': 'P', 'PTHR': 'T',
                                    'PTYR': 'Y', 'MELYS': 'K', 'MEARG': 'R', 'ACLYS': 'K',
                                    'KCX': 'K', 'LYE': 'K'}

            def _catch(func, *args, handle=lambda e : e, **kwargs):
                try:
                    return func(*args, **kwargs)
                except Exception as e:
                    print(f"Could not convert {e} to a 1 letter code: using 'X' instead")
                    return 'X'

            pdb_seqs = {}
            for chain in list(OrderedDict.fromkeys(subset_data['chain'])):
                tmp_data = subset_data[subset_data['chain'] == chain]
                pdb_seqs[chain] = ''.join([_catch(lambda : protein_letters_dict[a.upper()]) for a in list(tmp_data['resname'])])

        except Exception as e:
            #report_error_to_file('Ensemble - ESM 1', path, str(e))
            print(f'Ensemble ESM Calculation: 1 - could not extract the sequence from the protein file given: {e}')
            return pd.DataFrame(columns=["Chain", "Resid", "esm"])


        # 2: Extract local sequences based on the overall chain, submit sections to esm, record to dataframe
        df_esm = pd.DataFrame(columns=["Chain", "Resid", "esm"])
        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        esm2_model_650M, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
        esm2_model_650M = esm2_model_650M.to(device)
        batch_converter = alphabet.get_batch_converter()
        esm2_model_650M.eval()

        for idx, (lys_chain, lys_num) in enumerate(zip(list_chains, lys_res_nums)):
            try:
                seq = pdb_seqs[lys_chain]
                chain_shift_val = int(M.data[M.data['chain'] == lys_chain]['resid'].iloc[0]) - 1
                seq_lys_index = lys_num - chain_shift_val - 1

                start_idx = seq_lys_index - num_add_aa
                end_idx = seq_lys_index + num_add_aa + 1
                start_null = 0
                end_null = 0
                
                if start_idx < 0:
                    start_null = - start_idx
                    start_idx = 0

                if end_idx >= len(seq):
                    end_null = end_idx - len(seq)
                    end_idx = len(seq)

                seq = ('X' * start_null) + seq[start_idx:end_idx] + ('X' * end_null)

                data = [("1", seq)]
                batch_labels, batch_strings, batch_tokens = batch_converter(data)
                batch_tokens = batch_tokens.to(device)
                with torch.no_grad():
                    results = esm2_model_650M(batch_tokens, repr_layers=[33], return_contacts=True)
                seq_encode_tokens = results["representations"][33]
                #seq_encode_tokens = seq_encode_tokens.to('cpu')
                #print(seq_encode_tokens.shape)

                
                esm_proj_linear = torch.nn.Linear(1280, 512).to(device)
                esm_linear = esm_proj_linear(seq_encode_tokens)# [1, 43, 512]
                lstm_module = torch.nn.LSTM(512, 256, 2, batch_first=True, bidirectional=True).to(device)
                hidden_n = torch.zeros(4, esm_linear.size(0), 256).to(device)
                cell_n = torch.zeros(4, esm_linear.size(0), 256).to(device)
                seq_lstm, _ = lstm_module(esm_linear, (hidden_n, cell_n))
                seq_out = seq_lstm.mean(dim=1, keepdim=True)
                seq_out = seq_out.to('cpu')

                if self.include_modified:
                    df_esm = pd.concat([df_esm, pd.DataFrame([{'Chain': lys_chain, 'Resid': lys_num, 'esm': str(seq_out.tolist()[0][0]), 'Modified': list_modified[idx]}])], ignore_index=True)
                else:
                    df_esm = pd.concat([df_esm, pd.DataFrame([{'Chain': lys_chain, 'Resid': lys_num, 'esm': str(seq_out.tolist()[0][0])}])], ignore_index=True)

            except Exception as e:
                #report_error_to_file('Ensemble ESM 2', path, str(e))
                print(f'Ensemble ESM Calculation 2: Failed to encode sequence for uniprot: {path}, chain: {lys_chain}, resid: {lys_num}, error: {e}')
                if self.include_modified:
                    df_esm = pd.concat([df_esm, pd.DataFrame([{'Chain': lys_chain, 'Resid': lys_num, 'esm': None, 'Modified': list_modified[idx]}])], ignore_index=True)
                else:
                    df_esm = pd.concat([df_esm, pd.DataFrame([{'Chain': lys_chain, 'Resid': lys_num, 'esm': None}])], ignore_index=True)

        return df_esm


if __name__ == '__main__':
    
    df_prot = pd.DataFrame()
    E = Ensemble(df_proteins=df_prot, include_modified=False)
    df_esm = E.calculate_esm(path=f'1UBQ-alt-1.pdb')

    df_esm.to_pickle('testing_esm_csv.pkl')  #  pickle used if keeping tensors in output dataframe, csv fine if taking list through as can literal_eval
    df_esm_csv = pd.read_pickle('testing_esm_csv.pkl')
    print(df_esm_csv)
