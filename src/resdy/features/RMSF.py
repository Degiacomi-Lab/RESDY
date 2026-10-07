import os
import glob
import pandas as pd
import numpy as np
import biobox as bb
from .error_reporting import report_error_to_file


class RMSF():
    '''
    RMSF values for curated proteins.
    '''

    def __init__(self, df_proteins,
                 align_type='backbone',
                 include_modified = False,
                 aa_properties = {'non_modified_codes': ['LYS', 'LYSN'],
                                'modified_codes': ['LYE', 'KCX'],
                                'atom_select_names_nonmod': ['NZ'],
                                'atom_select_names_modified': ['NZ', 'N07']},
                 error_filename = 'measure_errors.txt'):
        '''
        Initialise the RMSF class, include any global variables that are required from measures in
        here.

        :param df_proteins: Dataframe including the information passed into measures from protein
            about which structures correspond to the uniprot codes.
        :type df_proteins: pandas.DataFrame
        :param align_type: The type of alignment to perform when calculating the RMSF values.
            Options:
            - 'local' - Aligns to backbone of lysine to calculate RMSF value for
            - 'backbone' (DEFAULT) - Aligns to backbone of full structure
        :type align_type: str
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
        self.df_proteins = df_proteins
        self.align_type = align_type
        self.include_modified = include_modified
        self.aa_properties = aa_properties
        self.error_filename = error_filename
        if self.error_filename != 'no_record':
            self.record_errors = True
        else:
            self.record_errors = False

        self.model_loaded = False


    def calculate(self, path):
        '''
        Calculate the root mean square fluctuation of the lysines within the protein over all the
        structures which have been curated for the Uniprot code. This method uses biobox but doesn't
        adjust the size of the window used to gain more matches.

        .. rubric:: Method

        Find all structures, loop over structures aligning and calculating deviations, calculate
        RMSF and return values as dataframe

        :param path: The path of the pdb file that the RMSF is being calculated for.
        :type path: str
        :returns:
            Dataframe with information on chain, residue number and rmsf output.
            Outline::

                Chain   Resid   rmsf
                x           x      x
        :rtype: pandas.DataFrame
        '''
        try:
            # check over measures to see if this has already been calculated as can just copy values due to being the same calculation each time
            file_loc = os.path.dirname(path)
            files = glob.glob(os.path.join(file_loc, "*pdb"))
            code = os.path.splitext(os.path.basename(path))[0]
            # The minimised copy of a structure is "<stem>_relaxed.pdb", and only the stem is
            # listed in the proteins dataframe. Splitting on '-' happens to drop the suffix of a
            # PDB code ('1UPT-alt-1_relaxed' -> '1UPT'), but an AlphaFold code keeps its hyphens
            # and so has to have the suffix removed explicitly, otherwise every relaxed AlphaFold
            # structure fails the lookup below and the feature returns nothing for it.
            if code.endswith('_relaxed'):
                code = code[:-len('_relaxed')]
            if 'AF-' not in code:
                code = code.split('-')[0]

            if 'Uniprot_Entry' not in self.df_proteins.columns:
                prot_info = [code]
            else:
                matches = self.df_proteins[self.df_proteins['PDB_Code'] == code]['Uniprot_Entry']
                if matches.empty:
                    raise KeyError(f'{code} is not in the dataframe given to RMSF calculations')
                uniprot_interest = matches.iloc[0]
                prot_info = list(self.df_proteins[self.df_proteins['Uniprot_Entry'] == uniprot_interest]['PDB_Code'])

            prot_match_exists = [a for a in files if (os.path.splitext(os.path.basename(a))[0].split('-')[0] in prot_info)
                                    or (os.path.basename(a).split('.')[0] in prot_info)]

            M = bb.Molecule()
            M.import_pdb(path, include_hetatm=True)

            M_ca = M.get_subset(M.atomselect('*', '*', 'CA', use_resname=True, get_index=True)[1])
            if self.include_modified: _, idx_n_res_interest = M.atomselect('*', (self.aa_properties['non_modified_codes'] + self.aa_properties['modified_codes']), 'CA', use_resname=True, get_index=True)
            else: _, idx_n_res_interest = M.atomselect('*', self.aa_properties['non_modified_codes'], 'CA', use_resname=True, get_index=True)
            list_modified = list(a in self.aa_properties['modified_codes'] for a in list(M.get_subset(idx_n_res_interest).data['resname']))
            df_rmsf = M.data.loc[idx_n_res_interest, ['chain', 'resid']].rename(columns={'chain': 'Chain', 'resid': 'Resid'}).reset_index(drop=True)
            if self.include_modified: df_rmsf = df_rmsf.assign(**{'Modified': list_modified})

            M_df_compare = M.data[M.data['name'] == 'CA'][['resname', 'chain', 'resid']].reset_index(drop=True)
        except Exception as e:
            if self.record_errors: report_error_to_file('RMSF 1', path, str(e), self.error_filename)
            print(f'RMSF Calculation 1: Failed to find other protein structures and get reference for uniprot: {path}, error: {e}')
            return pd.DataFrame(columns=['Chain', 'Resid', 'rmsf'])

        try:
            P = bb.Molecule()
            prot_matches = []
            for file in prot_match_exists:
                T = bb.Molecule()
                T.import_pdb(file, include_hetatm=True)
                T_df_compare = T.data[T.data['name'] == 'CA'][['resname', 'chain', 'resid']].reset_index(drop=True)
                if M_df_compare.equals(T_df_compare):
                    prot_matches.append(file)

            if len(prot_matches) > 1:
                P = bb.Molecule()
                P.import_pdb(prot_matches[0])
                for f in prot_matches[1:]:
                    P2 = bb.Molecule()
                    P2.import_pdb(f)
                    P2_xyz = P2.get_xyz()
                    P.add_xyz(P2_xyz)

                df_out = pd.DataFrame()
                res_interest = self.aa_properties['non_modified_codes'] + self.aa_properties['modified_codes']

                # A 'backbone' alignment is on the whole structure, so neither the reference
                # atoms nor the superposition they drive depend on the residue being measured.
                # Both are therefore done once here rather than once per residue. Note that
                # rmsd_one_vs_all(align=True) rewrites the coordinates in place, so repeating it
                # was superposing an already superposed ensemble. A 'local' alignment is on the
                # backbone of the residue itself and does stay inside the loop.
                align_once = self.align_type != 'local'
                alignment_failed = False
                if align_once:
                    try:
                        _, idx_ref = P.atomselect('*', '*', ["C", "CA", "N", "O"], get_index=True)
                        P.rmsd_one_vs_all(0, indices=idx_ref, align=True)
                    except Exception:
                        # The superposition used to be attempted inside the residue loop, where a
                        # failure left every residue with a NaN rather than discarding the whole
                        # structure. That is preserved by flagging it here and letting each
                        # residue take the same path it did before.
                        alignment_failed = True

                for i, r in P.data.drop_duplicates(subset=['chain', 'resid']).loc[P.data['resname'].isin(res_interest), ['chain', 'resid']].iterrows():
                    chain = r['chain']
                    resid = r['resid']
                    try:
                        if alignment_failed:
                            raise RuntimeError('the whole-structure superposition failed')
                        if not align_once:
                            _, idx_ref = P.atomselect('*', resid, ["C", "CA", "N", "O"], get_index=True)
                            P.rmsd_one_vs_all(0, indices=idx_ref, align=True)
                        _, idx_target = P.atomselect(chain, resid, self.aa_properties['atom_select_names_nonmod'], get_index=True)

                        rmsf = P.rmsf(indices=idx_target)[0]

                        df_out = pd.concat([df_out, pd.DataFrame([{'Chain': chain, 'Resid': resid, 'rmsf': rmsf}])])
                    except Exception:
                        df_out = pd.concat([df_out, pd.DataFrame([{'Chain': chain, 'Resid': resid, 'rmsf': np.nan}])])
                        pass

                df_rmsf = df_rmsf.merge(df_out, how='left', on=['Chain', 'Resid'])

            else:
                print('>> The only matching structure is itself, therefore rmsf cannot be calculated')
                df_rmsf = df_rmsf.assign(**{'rmsf': np.nan})

        except Exception as e:
            if self.record_errors: report_error_to_file('RMSF 2', path, str(e), self.error_filename)
            print(f'RMSF Calculation 2: Failed to calculate RMSF for uniprot: {path}, error: {e}')
            return pd.DataFrame(columns=['Chain', 'Resid', 'rmsf'])

        return df_rmsf.reset_index(drop=True).rename(columns={'resid': 'Resid', 'chain': 'Chain'})


if __name__ == '__main__':
    outdir = 'demo'
    df_prot = pd.read_csv(f'{outdir}{os.sep}proteins_demo.csv')
    rmsf = RMSF(df_proteins=df_prot, include_modified=False)
    print(rmsf.calculate(path=f'{outdir}{os.sep}curated{os.sep}4DCN-alt-1.pdb'))
