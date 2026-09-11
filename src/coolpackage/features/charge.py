import os
import pandas as pd
import biobox as bb
from .error_reporting import report_error_to_file
from collections import OrderedDict


class Charge():
    '''
    Charge values for structures.
    '''

    def __init__(self, include_modified = False,
                 aa_properties = {'non_modified_codes': ['LYS', 'LYSN'],
                                'modified_codes': ['LYE', 'KCX'],
                                'atom_select_names_nonmod': ['NZ'],
                                'atom_select_names_modified': ['NZ', 'N07']},
                 error_filename = 'measure_errors.txt'):
        '''
        Initialise the Charge class, include any global variables that are required from measures in
        here.

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
    
    def calculate(self, path, num_add_aa=10):
        '''
        Calculate the Sequence Charge of the local sequence around a LYS of interest. This is a
        single value number representing the summation of the charges of the amino acids over the
        specified number of amino acids either side of the lysine.

        .. rubric:: Method

        Take the PDB file and extract the overall sequence using BioBox. Identify all lysines within
        the structure and any shift that has taken place in the PDB file compared to the Uniprot
        sequence. Extract sequences for the lysine of interest and calculate a value for the charge
        based on the summation of charged residues within the sequence.

        :param path: The path of the pdb file that the sequence charge is being calculated for.
        :type path: str
        :param num_add_aa: The number of amino acids to include either side of the lysine of
            interest. An optional parameter which is set to 10 by default.
        :type num_add_aa: int
        :returns:
            Dataframe with information on chain, residue number and seqcharge output.
            Outline::

                Chain   Resid   seqcharge
                x           x           x
        :rtype: pandas.DataFrame

        .. todo::

           Add the ability to use 3 letter codes and their charges rather than the 1 letter codes, which
           may run into problems when modified residues are used. The change would be to stop converting
           to the classic 1 letter code sequence, and instead take the list of residues and map the
           charges onto it to sum (GW, 16.04.25).
        '''
        # 1: Extract the overall sequence for the protein given
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


            protein_letters_dict = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
                                    'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G',
                                    'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
                                    'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
                                    'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
                                    'HIE': 'H', 'HID': 'H', 'HIP': 'H', 'LYN': 'K',
                                    'ASX': 'B', 'GLX': 'Z', 'SEC': 'U', 'PYL': 'O',
                                    'XAA': 'X', 'XLE': 'J', 'PSER': 'p', 'PTHR': 't',
                                    'PTYR': 'y', 'MELYS': 'k', 'MEARG': 'r', 'ACLYS': 'k',
                                    'LYSN': 'K'}
            # KCX, LYE and LYSN (neutral lysine from gromacs) down as X so that they are not treated as positive K when calculations aren't including modified lysines
            # however when including modified, need to consider these K for checking purposes, gets tricky when considering near lysines that are all modified...
            if self.include_modified: protein_letters_dict['KCX'] = 'K'; protein_letters_dict['LYE'] = 'K'
            else: protein_letters_dict['KCX'] = 'X'; protein_letters_dict['LYE'] = 'X' 

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
            if self.record_errors: report_error_to_file('Seqcharge 1', path, str(e), self.error_filename)
            print(f'SeqCharge Calculation: 1 - could not extract the sequence from the protein file given: {e}')
            return pd.DataFrame(columns=["Chain", "Resid", "seqcharge"])

        # 2: Extract local sequences based on the overall chain, calculate charge score and add to output
        df_seqcharge = pd.DataFrame(columns=["Chain", "Resid", "seqcharge"])

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

                seq = ('-' * start_null) + seq[start_idx:end_idx] + ('-' * end_null)
                seq_split = list(seq)
                if seq_split[num_add_aa] != 'K':
                    print(f'A lysine was not found at the desired position {lys_num} on chain {lys_chain} '
                          f'read in for PDB file {path}; sequence -> {seq}')
                    continue

                pos_aa = ['K', 'H', 'R']
                neg_aa = ['D', 'E']
                count = 0
                for aa in seq:
                    if aa in pos_aa:
                        count += 1
                    elif aa in neg_aa:
                        count -= 1

                if self.include_modified:
                    df_seqcharge = pd.concat([df_seqcharge, pd.DataFrame([{'Chain': lys_chain, 'Resid': lys_num, 'seqcharge': count, 'Modified': list_modified[idx]}])], ignore_index=True)
                else:
                    df_seqcharge = pd.concat([df_seqcharge, pd.DataFrame([{'Chain': lys_chain, 'Resid': lys_num, 'seqcharge': count}])], ignore_index=True)

            except Exception as e:
                if self.record_errors: report_error_to_file('Seqcharge 2', path, str(e), self.error_filename)
                print(f'SeqCharge Calculation 2: Could not calculate a charge for lysine at position {lys_num}, error: {e}')
                if self.include_modified:
                    df_seqcharge = pd.concat([df_seqcharge, pd.DataFrame([{'Chain': lys_chain, 'Resid': lys_num, 'seqcharge': None, 'Modified': list_modified[idx]}])], ignore_index=True)
                else:
                    df_seqcharge = pd.concat([df_seqcharge, pd.DataFrame([{'Chain': lys_chain, 'Resid': lys_num, 'seqcharge': None}])], ignore_index=True)

        return df_seqcharge



if __name__ == '__main__':
    C = Charge(include_modified=True)
    print(C.calculate(path=f'result{os.sep}curated{os.sep}1UBQ-alt-1.pdb'))
