import os
import glob
import hashlib
import pandas as pd
import numpy as np
import biobox as bb
from .error_reporting import report_error_to_file

# The RMSF of a residue is a property of an ensemble of structures, so every file of an
# ensemble has the same result. The three caches below make the work proportional to the
# number of files instead of its square. They are held by the module, so that a worker
# process keeps them from one structure to the next.

#: {folder: (modification time of the folder, the pdb files it holds)}
_FOLDER_FILES = {}
#: {file: (modification time, size, digest of its CA table)}
_CA_SIGNATURES = {}
#: {(files considered, digest of the CA table, settings): (result, True when the ensemble
#: holds a single structure)}
_ENSEMBLE_RESULTS = {}


def clear_cache():
    '''
    Forget every folder listing, CA signature and ensemble result held by this module.
    '''
    _FOLDER_FILES.clear()
    _CA_SIGNATURES.clear()
    _ENSEMBLE_RESULTS.clear()


def _folder_files(folder):
    '''
    The pdb files of a folder, listed again only when the folder has changed.

    :param folder: folder holding the curated structures.
    :type folder: str
    :returns: paths of the pdb files, as ``glob`` lists them.
    :rtype: list
    '''
    stamp = os.stat(folder or '.').st_mtime_ns
    cached = _FOLDER_FILES.get(folder)
    if cached is None or cached[0] != stamp:
        cached = (stamp, glob.glob(os.path.join(folder, "*pdb")))
        _FOLDER_FILES[folder] = cached
    return cached[1]


def _ca_table(M):
    '''
    Residue name, chain and residue number of every CA atom of a structure. Two files belong to
    the same ensemble when their tables are equal.

    :param M: the structure.
    :type M: biobox.Molecule
    :rtype: pandas.DataFrame
    '''
    return M.data[M.data['name'] == 'CA'][['resname', 'chain', 'resid']].reset_index(drop=True)


def _digest(ca_table):
    '''
    A digest that is equal for two CA tables exactly when the tables are equal.

    :param ca_table: output of :func:`_ca_table`.
    :type ca_table: pandas.DataFrame
    :rtype: str
    '''
    h = hashlib.sha1()
    h.update(repr((ca_table.shape, [str(t) for t in ca_table.dtypes])).encode())
    h.update(pd.util.hash_pandas_object(ca_table, index=False).to_numpy().tobytes())
    return h.hexdigest()


def _ca_signature(path, M=None):
    '''
    Digest of the CA table of a file. The file is read only if it has not been read before, or
    has changed since.

    :param path: pdb file.
    :type path: str
    :param M: the file already loaded with its heteroatoms, to save reading it.
    :type M: biobox.Molecule, None
    :rtype: str
    '''
    info = os.stat(path)
    stamp = (info.st_mtime_ns, info.st_size)
    cached = _CA_SIGNATURES.get(path)
    if cached is not None and cached[:2] == stamp:
        return cached[2]
    if M is None:
        M = bb.Molecule()
        M.import_pdb(path, include_hetatm=True)
    signature = _digest(_ca_table(M))
    _CA_SIGNATURES[path] = stamp + (signature,)
    return signature


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


    def _settings_key(self):
        '''
        The settings an RMSF result depends on, in a form that can key a dictionary.
        '''
        return (self.align_type, bool(self.include_modified),
                tuple(self.aa_properties['non_modified_codes']),
                tuple(self.aa_properties['modified_codes']),
                tuple(self.aa_properties['atom_select_names_nonmod']))


    def calculate(self, path, M=None):
        '''
        Calculate the root mean square fluctuation of the lysines within the protein over all the
        structures which have been curated for the Uniprot code. This method uses biobox but doesn't
        adjust the size of the window used to gain more matches.

        .. rubric:: Method

        Find all structures of the Uniprot code, keep those whose CA atoms carry the same residue
        names, chains and residue numbers as the file given (the ensemble), superpose them,
        calculate the RMSF and return the values as a dataframe. The result is the same for every
        file of an ensemble, so it is calculated for the first of them and returned from memory
        for the others, and each file is read once to establish which ensemble it belongs to.

        :param path: The path of the pdb file that the RMSF is being calculated for.
        :type path: str
        :param M: The structure of ``path`` already loaded with its heteroatoms, which saves
            reading the file again. It is read from ``path`` when None.
        :type M: biobox.Molecule, None
        :returns:
            Dataframe with information on chain, residue number and rmsf output.
            Outline::

                Chain   Resid   rmsf
                x           x      x
        :rtype: pandas.DataFrame
        '''
        try:
            file_loc = os.path.dirname(path)
            files = _folder_files(file_loc)
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

            # the file is read here only when it was not handed over and its ensemble is not
            # known yet
            cached = _CA_SIGNATURES.get(path)
            info = os.stat(path)
            if M is None and (cached is None or cached[:2] != (info.st_mtime_ns, info.st_size)):
                M = bb.Molecule()
                M.import_pdb(path, include_hetatm=True)
            signature = _ca_signature(path, M)

            ensemble_key = (tuple(prot_match_exists), signature, self._settings_key())
            known = _ENSEMBLE_RESULTS.get(ensemble_key)
            if known is not None:
                df_known, single_structure = known
                if single_structure:
                    print('>> The only matching structure is itself, therefore rmsf cannot be calculated')
                return df_known.copy()

            if M is None:
                M = bb.Molecule()
                M.import_pdb(path, include_hetatm=True)

            if self.include_modified: _, idx_n_res_interest = M.atomselect('*', (self.aa_properties['non_modified_codes'] + self.aa_properties['modified_codes']), 'CA', use_resname=True, get_index=True)
            else: _, idx_n_res_interest = M.atomselect('*', self.aa_properties['non_modified_codes'], 'CA', use_resname=True, get_index=True)
            list_modified = list(a in self.aa_properties['modified_codes'] for a in list(M.get_subset(idx_n_res_interest).data['resname']))
            df_rmsf = M.data.loc[idx_n_res_interest, ['chain', 'resid']].rename(columns={'chain': 'Chain', 'resid': 'Resid'}).reset_index(drop=True)
            if self.include_modified: df_rmsf = df_rmsf.assign(**{'Modified': list_modified})
        except Exception as e:
            if self.record_errors: report_error_to_file('RMSF 1', path, str(e), self.error_filename)
            print(f'RMSF Calculation 1: Failed to find other protein structures and get reference for uniprot: {path}, error: {e}')
            return pd.DataFrame(columns=['Chain', 'Resid', 'rmsf'])

        try:
            # the ensemble: the files considered whose CA table is that of the file given
            prot_matches = [file for file in prot_match_exists if _ca_signature(file) == signature]
            single_structure = len(prot_matches) <= 1

            if not single_structure:
                P = bb.Molecule()
                P.import_pdb(prot_matches[0])
                for f in prot_matches[1:]:
                    P2 = bb.Molecule()
                    P2.import_pdb(f)
                    P2_xyz = P2.get_xyz()
                    P.add_xyz(P2_xyz)

                rows = []
                res_interest = self.aa_properties['non_modified_codes'] + self.aa_properties['modified_codes']

                # A 'backbone' alignment is on the whole structure, so neither the reference
                # atoms nor the superposition they drive depend on the residue being measured.
                # Both are therefore done once here rather than once per residue. Note that
                # rmsd_one_vs_all(align=True) rewrites the coordinates in place, so repeating it
                # would superpose an already superposed ensemble. A 'local' alignment is on the
                # backbone of the residue itself and does stay inside the loop.
                align_once = self.align_type != 'local'
                alignment_failed = False
                if align_once:
                    try:
                        _, idx_ref = P.atomselect('*', '*', ["C", "CA", "N", "O"], get_index=True)
                        P.rmsd_one_vs_all(0, indices=idx_ref, align=True)
                    except Exception:
                        # A failed superposition leaves every residue with a NaN rather than
                        # discarding the whole structure: it is flagged here and each residue
                        # takes the path below.
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

                        rows.append({'Chain': chain, 'Resid': resid, 'rmsf': rmsf})
                    except Exception:
                        rows.append({'Chain': chain, 'Resid': resid, 'rmsf': np.nan})

                df_out = pd.concat([pd.DataFrame([row]) for row in rows]) if rows else pd.DataFrame()
                df_rmsf = df_rmsf.merge(df_out, how='left', on=['Chain', 'Resid'])

            else:
                print('>> The only matching structure is itself, therefore rmsf cannot be calculated')
                df_rmsf = df_rmsf.assign(**{'rmsf': np.nan})

        except Exception as e:
            if self.record_errors: report_error_to_file('RMSF 2', path, str(e), self.error_filename)
            print(f'RMSF Calculation 2: Failed to calculate RMSF for uniprot: {path}, error: {e}')
            return pd.DataFrame(columns=['Chain', 'Resid', 'rmsf'])

        df_rmsf = df_rmsf.reset_index(drop=True).rename(columns={'resid': 'Resid', 'chain': 'Chain'})
        _ENSEMBLE_RESULTS[ensemble_key] = (df_rmsf, single_structure)
        return df_rmsf.copy()


if __name__ == '__main__':
    outdir = 'demo'
    df_prot = pd.read_csv(f'{outdir}{os.sep}proteins_demo.csv')
    rmsf = RMSF(df_proteins=df_prot, include_modified=False)
    print(rmsf.calculate(path=f'{outdir}{os.sep}curated{os.sep}4DCN-alt-1.pdb'))
