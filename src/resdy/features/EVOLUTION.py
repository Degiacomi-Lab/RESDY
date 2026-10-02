import os
from collections import OrderedDict
import pandas as pd
import biobox as bb
from .error_reporting import report_error_to_file


class EVOLUTION():
    '''
    ESM vectors for proteins.
    '''

    def __init__(self, include_modified = False,
                 aa_properties = {'non_modified_codes': ['LYS', 'LYSN'],
                                'modified_codes': ['LYE', 'KCX'],
                                'atom_select_names_nonmod': ['NZ'],
                                'atom_select_names_modified': ['NZ', 'N07']},
                 error_filename = 'measure_errors.txt'):
        '''
        Initialise the Evolution class, include any global variables that are required from measures
        in here.

        :param df_prot: Dataframe including the information passed into measures from protein about
            which structures correspond to the uniprot codes.
        :type df_prot: pandas.DataFrame
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

        try:
            import torch
            import torch.nn as nn
            import esm

            self.torch = torch
            self.nn = nn
            self.esm = esm

            self.esm_packages_available = True
        except ImportError as e:
            self.esm_packages_available = False
            raise ImportError(f'>> Failed to import packages required for esm calculations, '
                f'will not be able to calculate sequence features based on esm. Error: {e}') from e

        self.include_modified = include_modified
        self.aa_properties = aa_properties
        self.error_filename = error_filename
        if self.error_filename != 'no_record':
            self.record_errors = True
        else:
            self.record_errors = False

        self.model_loaded = False


    def check_esm_model_available(self):
        '''
        Check that the esm package is available, set model loaded to be False.
        '''
        if not self.esm_packages_available:
            raise ImportError(f'>> Failed to import the packages required (esm/torch) '
                              f'for esm calculations, esm will be removed from features.')
        self.model_loaded = False


    def _ensure_esm_model(self):
        '''
        Load in the model when running the process in the pool rather than before forking
        '''
        if self.model_loaded:
            return

        self.device = self.torch.device('cuda' if self.torch.cuda.is_available() else 'cpu')
        self.esm2_model_650M, self.alphabet = self.esm.pretrained.esm2_t33_650M_UR50D()
        self.esm2_model_650M = self.esm2_model_650M.to(device=self.device)
        self.batch_converter = self.alphabet.get_batch_converter()
        self.esm2_model_650M.eval()

        self.torch.manual_seed(25)
        if self.device.type == 'cuda':
            self.torch.cuda.manual_seed(25)
            self.torch.cuda.manual_seed_all(25)
        self.torch.backends.cudnn.deterministic = True
        self.torch.backends.cudnn.benchmark = False

        self.model_loaded = True


    def calculate(self, path, num_add_aa=20):
        '''
        Calculate the ESM LLM output for subset protein sequences obtained through searching over
        the path of the protein structure given.

        :param path: Path of the pdb file of interest
        :type path: str
        :param num_add_aa: Number of amino acids to go either side of the lysine of interest
        :type num_add_aa: int
        '''

        if not self.model_loaded:
            self.check_esm_model_available()

        self._ensure_esm_model()

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
                                    'KCX': 'K', 'LYE': 'K', 'LYSN': 'K',
                                    'ASH': 'D', 'GLH': 'E', 'CYX': 'C'}

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

            resid_to_index = {}
            for chain in pdb_seqs:
                tmp_data = subset_data[subset_data['chain'] == chain]
                for pos, res in enumerate(tmp_data['resid']):
                    resid_to_index[(chain, int(res))] = pos

        except Exception as e:
            if self.record_errors: report_error_to_file('Evolution - ESM 1', path, str(e), self.error_filename)
            print(f'Evolution ESM Calculation: 1 - could not extract the sequence from the protein file given: {e}')
            return pd.DataFrame(columns=['Chain', 'Resid', 'esm'])


        # 2: Extract local sequences based on the overall chain, submit sections to esm, record to dataframe
        df_esm = pd.DataFrame(columns=['Chain', 'Resid', 'evolution'])

        for idx, (lys_chain, lys_num) in enumerate(zip(list_chains, lys_res_nums)):
            try:
                seq = pdb_seqs[lys_chain]
                seq_lys_index = resid_to_index.get((lys_chain, int(lys_num)))
                if seq_lys_index is None:
                    print(f'> Res {lys_num} on chain {lys_chain} has no CA in {path}, skipped')
                    continue

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
                with self.torch.no_grad():
                    results = self.esm2_model_650M(batch_tokens, repr_layers=[33], return_contacts=False)
                    seq_encode_tokens = results["representations"][33]
                    seq_out = seq_encode_tokens[0, num_add_aa + 1, :]  # this bit of code can be used to extract a direct lysine representation without dimension reduction
                    seq_out = seq_out.to('cpu')

                if self.include_modified:
                    df_esm = pd.concat([df_esm, pd.DataFrame([{'Chain': lys_chain, 'Resid': lys_num, 'evolution': str(seq_out.tolist()), 'Modified': list_modified[idx]}])], ignore_index=True)
                else:
                    df_esm = pd.concat([df_esm, pd.DataFrame([{'Chain': lys_chain, 'Resid': lys_num, 'evolution': str(seq_out.tolist())}])], ignore_index=True)

            except Exception as e:
                if self.record_errors: report_error_to_file('Evolution ESM 2', path, str(e), self.error_filename)
                print(f'Evolution ESM Calculation 2: Failed to encode sequence for uniprot: {path}, chain: {lys_chain}, resid: {lys_num}, error: {e}')
                if self.include_modified:
                    df_esm = pd.concat([df_esm, pd.DataFrame([{'Chain': lys_chain, 'Resid': lys_num, 'evolution': None, 'Modified': list_modified[idx]}])], ignore_index=True)
                else:
                    df_esm = pd.concat([df_esm, pd.DataFrame([{'Chain': lys_chain, 'Resid': lys_num, 'evolution': None}])], ignore_index=True)

        return df_esm


if __name__ == '__main__':
    outdir = 'Demo'
    E = EVOLUTION(include_modified=False)
    print(E.calculate(path=f'result{os.sep}curated{os.sep}1UBQ-alt-1.pdb'))
