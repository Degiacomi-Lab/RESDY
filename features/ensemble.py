import os, sys
import glob
import pandas as pd
import numpy as np
import biobox as bb
from collections import OrderedDict
import MDAnalysis as mda
from MDAnalysis.analysis import rms, align
from features.error_reporting import report_error_to_file

try:
    import torch
    import torch.nn as nn
    import esm
except Exception as e:
    print(f'>> Failed to import packages required for esm calculations, will not be able to calculate sequence features based on esm. Error: {e}')


class Ensemble():
    '''
    Class to house the several metrics for calculating ensemble based features for a
    set of given structures
    '''

    def __init__(self, df_proteins, include_modified = False,
                 aa_properties = {'non_modified_codes': ['LYS', 'LYSN'],
                                'modified_codes': ['LYE', 'KCX'],
                                'atom_select_names_nonmod': ['NZ'],
                                'atom_select_names_modified': ['NZ', 'N07']},
                 error_filename = 'measure_errors.txt'):
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
        aa_properties -> dict
            Properties of the amino acid of interest to investigate modification sites for.
            Defaults to lysine for carbamylation. Properties are the 3 letter codes for
            non modified ('non_modified_codes') and modified ('modified_codes') and the atom
            names for non modified ('atom_select_names_nonmod') and modified ('atom_select_names_modified')
        error_filename : str
            Name of the text file passed through from overall measures to write any errors from
            calculating features out to.
        '''
        self.df_proteins = df_proteins
        self.include_modified = include_modified
        self.aa_properties = aa_properties
        self.error_filename = error_filename
        if self.error_filename != 'no_record': self.record_errors = True
        else: self.record_errors = False
        
        self.model_loaded = False


    def _initialise_esm_model(self):
        '''
        Initialise global parameters and models in calculating ESM values
        '''

        self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        self.esm2_model_650M, self.alphabet = esm.pretrained.esm2_t33_650M_UR50D()
        self.esm2_model_650M = self.esm2_model_650M.to(device=self.device)
        self.batch_converter = self.alphabet.get_batch_converter()
        self.esm2_model_650M.eval()

        torch.manual_seed(25)
        if self.device == 'cuda':
            torch.cuda.manual_seed(25)
            torch.cuda.manual_seed_all(25)
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = False

        self.model_loaded = True


    def calculate_rmsf(self, path):
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
        try:
            # check over measures to see if this has already been calculated as can just copy values due to being the same calculation each time
            file_loc = os.path.dirname(path)
            files = glob.glob(os.path.join(file_loc, "*pdb"))
            code = os.path.splitext(os.path.basename(path))[0]
            if 'AF-' not in code:
                code = code.split('-')[0]
            uniprot_interest = list(self.df_proteins[self.df_proteins['PDB_Code'] == code]['Uniprot_Entry'])[0]
            prot_info = list(self.df_proteins[self.df_proteins['Uniprot_Entry'] == uniprot_interest]['PDB_Code'])
            prot_match_exists = [a for a in files if os.path.splitext(os.path.basename(a))[0].split('-')[0] in prot_info]

            M = bb.Molecule()
            M.import_pdb(path, include_hetatm=True)
            
            M_ca = M.get_subset(M.atomselect('*', '*', 'CA', use_resname=True, get_index=True)[1])
            if self.include_modified: idx_n_res_interest = M_ca.atomselect('*', (self.aa_properties['non_modified_codes'] + self.aa_properties['modified_codes']), 'CA', use_resname=True, get_index=True)[1]
            else: idx_n_res_interest = M_ca.atomselect('*', self.aa_properties['non_modified_codes'], 'CA', use_resname=True, get_index=True)[1]
            list_modified = list(a in self.aa_properties['modified_codes'] for a in list(M_ca.data['resname']))
            df_rmsf = M_ca.data[['resid', 'chain']]
            if self.include_modified: df_rmsf = df_rmsf.assign(**{'Modified': list_modified})

            M_df_compare = M.data[M.data['name'] == 'CA'][['resname', 'chain', 'resid']].reset_index(drop=True)
        except Exception as e:
            if self.record_errors: report_error_to_file('RMSF 1', path, str(e), self.error_filename)
            print(f'RMSF Calculation 1: Failed to find other protein structures and get reference for uniprot: {path}, error: {e}')
            return pd.DataFrame(columns=['Chain', 'Resid', 'rmsf'])

        try:
            prot_matches = []
            for file in prot_match_exists:
                T = bb.Molecule()
                T.import_pdb(file, include_hetatm=True)
                T_df_compare = T.data[T.data['name'] == 'CA'][['resname', 'chain', 'resid']].reset_index(drop=True)
                if M_df_compare.equals(T_df_compare):
                    prot_matches.append(file)
            
            if len(prot_matches) > 1:
                prot_conf_unv = mda.Universe(prot_matches[0], prot_matches, format='PDB', dt=1.0)
                aligner = align.AlignTraj(prot_conf_unv, prot_conf_unv, select='protein and name CA', in_memory=True, ref_frame=0).run()
                prot_c_alphas = prot_conf_unv.select_atoms('protein and name CA')
                rmsf_calculator = rms.RMSF(prot_c_alphas).run()
                rmsf_vals = rmsf_calculator.results.rmsf
                df_rmsf = df_rmsf.assign(**{'rmsf': rmsf_vals})

            else:
                df_rmsf = df_rmsf.assign(**{'rmsf': np.NaN})

            df_rmsf = df_rmsf.iloc[idx_n_res_interest]
        except Exception as e:
            if self.record_errors: report_error_to_file('RMSF 2', path, str(e), self.error_filename)
            print(f'RMSF Calculation 2: Failed to calculate RMSF for uniprot: {path}, error: {e}')
            return pd.DataFrame(columns=['Chain', 'Resid', 'rmsf'])

        return df_rmsf.reset_index(drop=True).rename(columns={'resid': 'Resid', 'chain': 'Chain'})


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

        if not self.model_loaded: self._initialise_esm_model()

        # 1: Extract the sequence sections from the given protein structure
        try:
            M = bb.Molecule()
            M.import_pdb(path, include_hetatm=True)

            if self.include_modified:
                idx_nz = M.atomselect('*',
                                      (self.aa_properties['non_modified_codes'] + self.aa_properties['modified_codes']),
                                      self.aa_properties['atom_select_names_modified'],
                                      use_resname=True, get_index=True)[1]
                # due to wider selection criteria, possible to get more than 1 hit per residue of interest, remove duplicates
                key_res_chain = zip(list(M.data['resid'].values[idx_nz]), list(M.data['chain'].values[idx_nz]))
                pairs_seen, keep_pos = set(), []
                for pair, pos in zip(key_res_chain, range(len(idx_nz))):
                    if pair not in pairs_seen:
                        pairs_seen.add(pair)
                        keep_pos.append(pos)
                idx_nz = idx_nz[keep_pos]

            else:
                idx_nz = M.atomselect('*',
                                      self.aa_properties['non_modified_codes'],
                                      self.aa_properties['atom_select_names_nonmod'],
                                      use_resname=True, get_index=True)[1]

            lys_res_nums = list(M.data['resid'][idx_nz])
            list_chains = list(M.data['chain'][idx_nz])
            list_modified = list(a in self.aa_properties['modified_codes'] for a in list(M.data['resname'][idx_nz]))

            c_alpha_idxs = M.atomselect('*', '*', 'CA', use_resname=True, get_index=True)[1]
            subset_data = M.data.iloc[c_alpha_idxs]

            # for purposes of esm - needs to be assigned to original 20 AA codes
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
                                    'KCX': 'K', 'LYE': 'K', 'LSYN': 'K'}

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
            if self.record_errors: report_error_to_file('Ensemble - ESM 1', path, str(e), self.error_filename)
            print(f'Ensemble ESM Calculation: 1 - could not extract the sequence from the protein file given: {e}')
            return pd.DataFrame(columns=['Chain', 'Resid', 'esm'])


        # 2: Extract local sequences based on the overall chain, submit sections to esm, record to dataframe
        df_esm = pd.DataFrame(columns=['Chain', 'Resid', 'esm'])

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
                batch_labels, batch_strings, batch_tokens = self.batch_converter(data)
                batch_tokens = batch_tokens.to(device=self.device)
                with torch.no_grad():
                    results = self.esm2_model_650M(batch_tokens, repr_layers=[33], return_contacts=False)
                    seq_encode_tokens = results["representations"][33]
                    seq_out = seq_encode_tokens[0, num_add_aa + 1, :]  # this bit of code can be used to extract a direct lysine representation without dimension reduction
                    seq_out = seq_out.to('cpu')

                if self.include_modified:
                    df_esm = pd.concat([df_esm, pd.DataFrame([{'Chain': lys_chain, 'Resid': lys_num, 'esm': str(seq_out.tolist()), 'Modified': list_modified[idx]}])], ignore_index=True)
                else:
                    df_esm = pd.concat([df_esm, pd.DataFrame([{'Chain': lys_chain, 'Resid': lys_num, 'esm': str(seq_out.tolist())}])], ignore_index=True)

            except Exception as e:
                if self.record_errors: report_error_to_file('Ensemble ESM 2', path, str(e), self.error_filename)
                print(f'Ensemble ESM Calculation 2: Failed to encode sequence for uniprot: {path}, chain: {lys_chain}, resid: {lys_num}, error: {e}')
                if self.include_modified:
                    df_esm = pd.concat([df_esm, pd.DataFrame([{'Chain': lys_chain, 'Resid': lys_num, 'esm': None, 'Modified': list_modified[idx]}])], ignore_index=True)
                else:
                    df_esm = pd.concat([df_esm, pd.DataFrame([{'Chain': lys_chain, 'Resid': lys_num, 'esm': None}])], ignore_index=True)

        return df_esm


if __name__ == '__main__':
    outdir = 'Demo'
    df_prot = pd.read_csv(f'{outdir}{os.sep}proteins_demo.csv')
    E = Ensemble(df_proteins=df_prot, include_modified=False)
    print(E.calculate_esm(path=f'result{os.sep}curated{os.sep}1UBQ-alt-1.pdb'))

    #print(E.calculate_rmsf(path=f'{outdir}{os.sep}curated{os.sep}1M2E-alt-1.pdb'))
    #df_esm.to_pickle('testing_esm_csv.pkl')  #  pickle used if keeping tensors in output dataframe, csv fine if taking list through as can literal_eval
    #df_esm_csv = pd.read_pickle('testing_esm_csv.pkl')
    #print(df_esm_csv)
