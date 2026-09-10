import os
import re
import io
import glob
import shutil
import time
from contextlib import redirect_stdout
from multiprocessing import cpu_count
from multiprocessing import Manager
from multiprocessing.pool import Pool
from collections import OrderedDict
import requests
import pandas as pd
import numpy as np
from Bio.Align import PairwiseAligner, substitution_matrices
import biobox as bb
from . import alphafold as af
from . import patcher
from .helper import get_download_tool, ShutUp


class PDB(object):
    '''
    Class for taking an input dataframe of the proteins required to investigate, downloading the
    proteins required from both the PDB and alphafold and cleaning the structures up to ensure that
    these are all ready to be used within measurements for the prediction model input. Cleaning is
    currently done using Modeller.
    '''

    def __init__(self, outdir="result", gap=10, parallel = False,
                 PDB_only=False, include_hetatm=False,
                 resnames_of_interest = ['LYS'], minimise_af=True):
        '''
        Initialise the PDB class.

        :param outdir: The directory in which files should be downloaded and curated within
        :type outdir: str
        :param gap: The maximum gap that is allowed in the sequence for a structure that has been
            downloaded that patching will be done on. For structures with a gap in the sequence
            greater than this, the structure will be removed.
        :type gap: int
        :param parallel: Toggle to run the curation of pdb files in parallel rather than series
            (default False)
        :type parallel: bool
        :param PDB_only: Toggle for if you want to download a system from a list of PDB files (True)
            or from a Uniprot dataframe (False) created from the Uniprot class.
        :type PDB_only: bool
        :param include_hetatm: Toggleable option to allow hetatms to pass through biobox
        :type include_hetatm: bool
        :param resnames_of_interest: List of residues to investigate, only used for curating list of
            PLDDT values for the residues of interest here.
        :type resnames_of_interest: list
        :param minimise_af: Toggleable option to re-minimise AF structures in implicit solvent rather
            than vacuum as standard AF structures are. Default is set to True. 
        :type minimise_af: bool
        '''

        self.outdir = outdir
        self.PDB_only = PDB_only
        self.parallel = parallel
        self.include_hetatm = include_hetatm
        self.resnames_of_interest = resnames_of_interest
        self.minimise_af = minimise_af

        # create folder of curated protein structures
        self.curated_dir = os.path.join(outdir, "curated")
        if not os.path.exists(self.curated_dir):
            os.makedirs(self.curated_dir)

        # create working folder
        self.raw_dir = os.path.join(outdir, "conformations")
        if not os.path.exists(self.raw_dir):
            os.makedirs(self.raw_dir)

        # dataframe storing data
        if self.PDB_only:
            columns = ['PDB_Code']
            self.df = pd.DataFrame(columns=columns)
        else:
            columns = ['Uniprot_Entry', 'PDB_Code', 'Method', 'Resolution', 'Chains', 'Largest_Gap']
            self.df = pd.DataFrame(columns=columns)

        self.gap = gap


    def save_state(self, outname='proteins.csv'):
        '''
        Take the current dataframe of proteins that have been gathered (pdb.df) and write this to a
        csv file in the output directory for later use.

        :param outname: the desired name for the pdb.df csv file to be written as; default:
            'proteins.csv'
        :type outname: str

        .. rubric:: Example

        ::

            pdb.save_state()
        '''
        self.df.to_csv(os.path.join(self.outdir, outname), index_label=False, index=False)


    def load_state(self, fname, outdir="", gap=10, PDB_only=False):
        '''
        Allows part resuming of a partially run curation. Loads in the proteins.csv file curated in
        a previous run into the self.df dataframe. While it doesn't stop the programme trying to
        curate everything again, this isn't a problem as it won't re-curate data anyway if
        skip_if_true is set to True. Sets other necessary parameters for matching the previous run.
        If no output directory is given, the folder containing the csv file is used.

        :param fname: filename for the previous run protein.csv file
        :type fname: str
        :param outdir: path to the output directory to work in, if not output directory set, it
            tries to take the one where the proteins.csv file was saved
        :type outdir: str
        :param gap: max gap allowed in sequences to be able to be patched. Default 10
        :type gap: int
        :param PDB_only: run the curation based on pdb files alone without curating for full uniprot
            codes
        :type PDB_only: str

        .. todo::

           Consider inferring the PDB_only state from the format of the measures dataframe
           (GW, 30.07.26).
        '''
        if outdir == '':
            outdir = os.path.dirname(fname) if self.outdir == '' else self.outdir

        self.outdir = outdir
        self.curated_dir = os.path.join(outdir, 'curated')
        self.raw_dir = os.path.join(outdir, 'conformations')
        os.makedirs(self.curated_dir, exist_ok=True)
        os.makedirs(self.raw_dir, exist_ok=True)

        self.PDB_only = PDB_only
        self.gap = gap

        try:
            self.df = pd.read_csv(fname)
            # remove duplicates from the dataframe to avoid extra unnecessary calculations
            self.df = self.df.drop_duplicates()

        except Exception as e:
            print(f"Could not load csv file, error: {str(e)}")


    def gather_proteins(self, uniprot_df, skip_if_found=True):
        '''
        Parallel method for gathering protein structures from the given dataframe containing list of
        uniprot codes. This turns the input data into a list for input into the parallel pool,
        format of the list depends on if running as PDB_only or not. Output from the parallel run is
        concat into the overall output dataframe self.df.

        :param uniprot_df: Dataframe containing all the uniprot codes of interest, can also include
            the residues of interest along with this. Columns: 'Uniprot_Entry', 'PDB_Code' etc This
            comes from the output of uniprot.py
        :type uniprot_df: pandas.DataFrame
        :param skip_if_found: Set to true to not curate another structure for the given pdb if a
            curated file is found, if set to false it will ignore previously curated files and start
            again
        :type skip_if_found: bool
        '''
        if self.PDB_only:
            try:
                # when the inputs are PDB codes only, convert them to a dataframe
                dic = {'PDB_Code':uniprot_df}
                uniprot_df = pd.DataFrame(dic)
            except Exception as e:
                return e

        if self.parallel:
            n_cores_to_use = max(1, int(round(cpu_count() * 0.75)))
        else:
            n_cores_to_use = 1

        with Manager() as manager:
            lock = manager.Lock()
            items = []
            for i, r in uniprot_df.iterrows():
                try:
                    if not self.PDB_only:
                        method_obtained = r.get('Method', '')
                        resolution = r.get('Resolution', '')
                        chains = r.get('Chains', '')
                        items.append([[r['Uniprot_Entry'], r['PDB_Code'], method_obtained, resolution, chains], skip_if_found, lock])
                    else:
                        items.append([[r['PDB_Code']], skip_if_found, lock])
                except Exception as e:
                    print(f'>> In adding rows to the parallel pool, failed to add row {i} '
                          f'and therefore will not be created; error: {e}')
                    broken_prot = pd.DataFrame([{'Uniprot': r.get('Uniprot_Entry', ''),
                                                 'PDB': r.get('PDB_Code', ''),
                                                 'Error': f'unreadable file row setting parallel list; error: {str(e)}'}])
                    broken_path = os.path.join(self.outdir, 'uncuratable_pdb_files.csv')
                    broken_prot.to_csv(broken_path, mode='a', index=False, header=not os.path.exists(broken_path))


            with Pool(n_cores_to_use, maxtasksperchild=20) as pool:
                results = [t for t in pool.starmap(self._curate_row, items) if t is not None]

            if not results:
                print(f'>> No structures of the {len(uniprot_df)} uniprot codes could be curated, '
                      f'for more details, see uncuratable_pdb_files.csv')
                return
            df_tmp = pd.concat(results, ignore_index=True)

            self.df = pd.concat([self.df, df_tmp], ignore_index=True).reset_index(drop=True)


    def _curate_row(self, row_details, skip_if_found, lock):
        '''
        Function for use within the parallel gather proteins function for gathering pdb information
        for the given uniprot code within the row details. Checks if the pdb code contains 'AF-', if
        it does, follow the alphafold curation path, download the file and record the plddt scores
        for it. If 'AF-' is not in it, follow route for pdb file code. If skip_if_true set to True,
        it will check if any curated structures exist for that pdb and if there are, check alignment
        of sequence and return. If not it will call the functions to clean and split the pdb into
        curated structures.

        :param row_details:
        :type row_details: list
        :param skip_if_found: Set to true to not curate another structure for the given pdb if a
            curated file is found, if set to false it will ignore previously curated files and start
            again
        :type skip_if_found: bool
        :param lock: Lock provided from the multiprocessing manager to use to stop child processes
            writing to files at the same time.
        :type lock: lock
        '''
        print_statements = []
        start_time = time.time()
        if not self.PDB_only:
            uniprot_code, pdb_code, method_obtained, resolution, chains = row_details
        else:
            pdb_code, = row_details
            uniprot_code = 'Running PDB only'
            method_obtained, resolution, chains = '', '', ''

        print_statements.append(f'>> Start file curation for uniprot: {uniprot_code}, pdb code: {pdb_code}')

        if not self.PDB_only:
            print_statements.append(f"UNIPROT: {uniprot_code}, PDB: {pdb_code}")
        else:
            print_statements.append(f"PDB: {pdb_code}")

        def finish_curate_jobs(failed=False):
            '''
            Generalised function for finishing off the _curate_row() function output to terminal,
            placed into function as this will be used at several points in the _curate_row()
            function.

            :param failed: Determines the final statement to print on if the process failed (True)
                or not (False); default is False
            :type failed: bool
            '''

            time_taken = round((time.time() - start_time), 2)
            if not failed:
                print_statements.append(f'>> Finished; Uniprot: {uniprot_code}, PDB: {pdb_code}; file curated in {time_taken}s\n\n')
            else:
                print_statements.append(f'>> Finished; Uniprot: {uniprot_code}, PDB: {pdb_code}; file attempted curation in {time_taken}s\n\n')
            with lock:
                for statement in print_statements:
                    print(statement)

        if pdb_code[:2] == "AF":

            if skip_if_found:
                files=[os.path.basename(c).split(".")[0] for c in glob.glob(os.path.join(self.curated_dir, "*pdb"))]
                if pdb_code in files:
                    print_statements.append(f">> curated {pdb_code} PDB found, continuing...")
                    if not self.PDB_only:
                        data = {'Uniprot_Entry': uniprot_code, 'PDB_Code': pdb_code, 'Method': 'Predicted', 'Resolution': np.nan, 'Chains': "A"}
                    else:
                        data = {'PDB_Code': pdb_code}

                    finish_curate_jobs(failed=False)
                    return pd.DataFrame.from_records(data, index=[0])

            try:
                out_print_trap_af_download = io.StringIO()
                with redirect_stdout(out_print_trap_af_download):
                    try:
                        af.download_AF_struc(pdb_code,outfolder=self.outdir)
                    except Exception as e:
                        print_statements.append(f'>> Failed to download AF structure calling af.download_AF_struc() with error: {str(e)}')
                print_statements.append(out_print_trap_af_download.getvalue())
            except Exception as e:
                print_statements.append(f">> FAILED on calling download_AF_struc: {str(e)}")
                finish_curate_jobs(failed=True)
                return None

            # check if the AlphaFold file contains ATOM statements
            af_filename = os.path.join(self.curated_dir, f"{pdb_code}.pdb")
            if not os.path.exists(af_filename):
                print_statements.append(f'Failed: no AF structure was written for {af_filename}')
                finish_curate_jobs(failed=True)
                return None

            fin = open(af_filename, "r")
            pdb_data_present = False
            for line in fin:
                if line.startswith("ATOM"):
                    pdb_data_present = True
                    break
            fin.close()

            if pdb_data_present:
                if not self.PDB_only:
                    data = {'Uniprot_Entry': uniprot_code, 'PDB_Code': pdb_code, 'Method': 'Predicted', 'Resolution': np.nan, 'Chains': "A"}
                else:
                    data = {'PDB_Code': pdb_code}

                try:
                    # find the PLDDT codes for AF structures
                    out_print_trap_af_plddt = io.StringIO()
                    with redirect_stdout(out_print_trap_af_plddt):
                        with lock:
                            af.find_af_plddt(pdb_code,outfolder=self.outdir, resnames=self.resnames_of_interest)
                    print_statements.append(out_print_trap_af_plddt.getvalue())
                except Exception as e:
                    print(f'> Failed to calculate PLDDT values for Af structure: {pdb_code}; Error: {str(e)}')

                try:
                    if self.minimise_af:
                        # apply minimisation to AF structure
                        out_print_trap_af_minimisation = io.StringIO()
                        with redirect_stdout(out_print_trap_af_minimisation):
                            af.apply_minimisation(pdb_code, outfolder=self.outdir)
                        print_statements.append(out_print_trap_af_minimisation.getvalue())
                except Exception as e:
                    print(f'>> Failed on af.apply_minimisation() for pdb {pdb_code} with error: {str(e)}')
                    finish_curate_jobs(failed=True)
                    return None

                finish_curate_jobs(failed=False)
                return pd.DataFrame.from_records(data, index=[0])

            else:
                print_statements.append(">> FAILED: structure not found in AlphaFold database")
                try:
                    os.remove(af_filename)
                except Exception as e:
                    pass
                finish_curate_jobs(failed=True)

        else:
            if isinstance(chains, str):
                chains = [c for c in chains.split('/') if c]
            if skip_if_found:
                files=[os.path.basename(c).split("-")[0] for c in glob.glob(os.path.join(self.curated_dir, "*pdb"))]
                if pdb_code in files:
                    print_statements.append(f">> curated {pdb_code} PDB found, continuing...")
                    if not self.PDB_only:
                        matched_curated_files = [os.path.basename(a) for a in glob.glob(os.path.join(self.curated_dir, "*pdb")) if pdb_code.upper() == os.path.splitext(os.path.basename(a))[0].split('-')[0]]
                        for file in matched_curated_files:
                            out_print_trap_check = io.StringIO()
                            with redirect_stdout(out_print_trap_check):
                                try:
                                    self._check_curated_structure(os.path.join(self.curated_dir, file), uniprot_code, chains)
                                except Exception as e:
                                    print_statements.append(f'>> Failed checking structure on previously curated structure with error {str(e)}')
                            print_statements.append(out_print_trap_check.getvalue())
                        data = {'Uniprot_Entry': uniprot_code, 'PDB_Code': pdb_code, 'Method': method_obtained,
                                'Resolution': resolution, 'Chains': '/'.join(chains), 'Largest_Gap': np.NaN}
                    else:
                        data = {'PDB_Code': pdb_code}

                    finish_curate_jobs(failed=False)
                    return pd.DataFrame.from_records(data, index=[0])

            # load, clean, and split it in alternate conformations
            out_print_trap = io.StringIO()
            data = None
            try:
                with redirect_stdout(out_print_trap):
                    if not self.PDB_only:
                        largest_gap = self.clean_and_split_pdb(pdb_code, uniprot_code, chains)
                        data = {'Uniprot_Entry': uniprot_code, 'PDB_Code': pdb_code, 'Method': method_obtained,
                                'Resolution': resolution, 'Chains': '/'.join(chains), 'Largest_Gap': largest_gap}
                    else:
                        largest_gap = self.clean_and_split_pdb(pdb_code)
                        data = {'PDB_Code': pdb_code, 'Largest_Gap': largest_gap}
            except Exception as e:
                print_statements.append(out_print_trap.getvalue())
                broken_prot = pd.DataFrame([{'Uniprot': uniprot_code, 'PDB': pdb_code, 'Error': f'{type(e).__name__}: {str(e)}'}])
                broken_path = os.path.join(self.outdir, 'uncuratable_pdb_files.csv')
                with lock:
                    broken_prot.to_csv(broken_path, mode='a', index=False, header=not os.path.exists(broken_path))
                print_statements.append(f'>> Failed clean_and_split_pdb() with error: {str(e)}')
                finish_curate_jobs(failed=True)
                return None

            print_statements.append(out_print_trap.getvalue())
            finish_curate_jobs(failed=False)
            return pd.DataFrame.from_records(data, index=[0])


    def clean_and_split_pdb(self, pdb, uniprot_code = '', chains=[]):
        '''
        Download the pdb file from rcsb, clean the structure, split it based on alternative
        conformations, find if there are gaps in the sequences, if more than the max gap set in
        class definition, remove files, if smaller, patch the curated files. The results are saved
        into files: ``[outfolder]/conformations/[PDB code]-clean.pdb``.

        :param pdb: code of the pdb file to curate
        :type pdb: str
        :param uniprot_code: corresponding uniprot code for the pdb file, not required but used for
            realignment. Default ''
        :type uniprot_code: str
        :param chains: list of chains for the pdb interested in, not required but used if pdb file
            goes over more than 1 uniprot code. Default: []
        :type chains: list

        .. todo::

           Renumbering is only possible for structures that have a Uniprot code associated with them, so
           that a sequence to align to is available. Find a way of making it work with PDB_only
           (24.03.26).
        '''
        try:
            #download and clean the structure
            self.download_pdb(pdb)
            replacement_dict = self.clean(pdb, uniprot_code=uniprot_code)

            #splits into all alternative conformations into independent structures
            self.split_struc_nmr(pdb)
            self.split_struc_alt_aa(pdb)

        except Exception as e:
            raise Exception(f'Error cleaning {pdb}: {str(e)}') from e


        files = glob.glob(os.path.join(self.raw_dir, f"*{pdb}*pdb"))
        test = False
        largest_gap = np.NaN
        for cnt, f in enumerate(files):
            mypath = os.path.split(f)[0]
            fasta_loc = os.path.join(mypath, f"{pdb}.fasta")

            try:

                fname, largest_gap = patcher.curate(f, fasta_loc, outdir=self.curated_dir,
                                       gap=self.gap, include_hetatm=self.include_hetatm)
                if len(replacement_dict) > 0:
                    reverse_replacement_dict = dict((v,k) for k,v in replacement_dict.items())
                    self.replace_chains(fname, reverse_replacement_dict)

                if not self.PDB_only:
                    try:
                        self._align_resnum_uniprot(uniprot_code, fname, chains)
                    except Exception as e:
                        print(f'>> Renumbering skipped for conformer {cnt} with error: {str(e)}')

                test = True

            except Exception as e:
                tmp_name = f'tmp_{os.path.splitext(os.path.basename(f))[0]}'
                shutil.rmtree(os.path.join(self.curated_dir, tmp_name), ignore_errors=True)
                print(f">> Patching failed for conformer {cnt}. Error: {str(e)}")
                continue

        if not test:
            raise Exception("Patching failed for all conformers")

        return largest_gap


    def download_pdb(self, pdb):
        '''
        Download the required PDB file for the PDB code specified. This is used for the coordinates
        of the protein to extract the featurised data for each protein from.

        .. rubric:: Method

        Identify the download tool that is available for use. Check if PDB file has already been
        downloaded, if so skip. Otherwise download the PDB file from the RCSB website.

        :param pdb: the PDB code for the PDB file to be downloaded
        :type pdb: str

        .. rubric:: Example

        ::

            >>> self.download_pdb('1ubq')
        '''
        files=glob.glob(os.path.join(self.raw_dir, "*pdb"))
        test_file = os.path.join(self.raw_dir, f'{pdb}.pdb')
        if test_file not in files:
            print(f">> Downloading PDB {pdb}")
            try:
                response = requests.get(url=f'https://files.rcsb.org/download/{pdb}.pdb', timeout=20)
                response.raise_for_status()
                open(os.path.join(self.raw_dir, f'{pdb}.pdb'), 'wb').write(response.content)

            except Exception as e:
                print(f'Error downloading file. {str(e)}')
        else:
            print(f'PDB file ({pdb}) has previously been downloaded, using previous copy.')


    def download_fasta(self, pdb):
        '''
        Download the required fasta file for the PDB code specified. This is used to align the
        structures and check for missing residues and extract sequences for comparison.

        .. rubric:: Method

        Identify the download tool that is available for use. Check if fasta file has already been
        downloaded, if so skip. Otherwise download the fasta file from the RCSB website.

        :param pdb: the PDB code for the fasta file to be downloaded
        :type pdb: str

        .. rubric:: Example

        ::

            >>> self.download_fasta('1ubq')
        '''
        def get_data():
            try:
                print(f">> downloading FASTA for {pdb}")

                response = requests.get(url=f'https://www.rcsb.org/fasta/entry/{pdb}/download', timeout=20)
                response.raise_for_status()
                open(os.path.join(self.raw_dir, f'{pdb}.fasta'), 'wb').write(response.content)

            except Exception as e:
                print(f'>> Failed downloading FASTA sequence for chain name comparison for {pdb}, error: {str(e)}')
                raise Exception(f'Failed downloading FASTA sequence for chain name comparison.{str(e)}') from e
        files=[c.split(os.sep)[-1][:4] for c in glob.glob(os.path.join(self.raw_dir, "*.fasta"))]
        if pdb not in files:
            get_data()
        elif os.path.getsize(os.path.join(self.raw_dir, f'{pdb}.fasta')) == 0:
            get_data()
        else:
            print(f'Fasta file for {pdb} previously downloaded and has data, using previous copy.')


    def clean(self, pdb, uniprot_code = ''):
        '''
        Rename the protein's chains during the cleaning process so that they match the chain names
        given in the FASTA file. This is required as pdb files name their chains using the 'auth'
        name and fasta with the normal chain name. Therefore to avoid confusion we rename them all
        to what is used in the fasta file.

        Iterates over the lines of the file performing checks:

        - insertion codes -> make note and not patch structure if present
        - modified lysines (KCX) -> remove CO2 and rename to LYS
        - modified cysteines (SEC) -> replace SE with S and rename to CYS
        - modified methionine (MSE) -> replace SE with S and rename to MET
        - element codes -> if not present in pdb file, add these in based on guess from atomtype col
        - neglect HETATMs unless metal ions and hydrogens
        - check for modified residues in structure, convert back to

        :param pdb: pdb code for the structure of interest
        :type pdb: str
        :param uniprot_code: Uniprot code corresponding to structure of interest - default is '' to
            allow for pdb only cases
        :type uniprot_code: str
        :returns: Dictionary containing the mapping between original chains and renamed chains
            produced by get_chain_replacement(pdb)
        :rtype: dict
        '''
        try:
            replacement_dict = self.get_chain_replacement(pdb)
            path = os.path.join(self.raw_dir, f"{pdb}.pdb")
            if len(replacement_dict) > 0:
                self.replace_chains(path, replacement_dict)

        except Exception as e:
            raise Exception(f'Error renaming chains. {str(e)}')

        #Next it opens and starts reading the .pdb file and starts writing a new file with the ending '-clean.pdb'.
        list_of_metals = ['ZN', 'NI', 'CU', 'FE', 'MG', 'MN', 'NA', 'K', 'CA', 'CO', 'CL', 'MO']
        standard_resids = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C', 'GLN': 'Q',
                           'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
                           'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W',
                           'TYR': 'Y', 'VAL': 'V'}
        standard_resids_inv = {v: k for k, v in standard_resids.items()}
        fasta_chains = {}
        uniprot_fasta = ''
        list_prev_mod_resids = {}

        read_file_path = os.path.join(self.raw_dir, f"{pdb}.pdb")
        read_file = open(read_file_path)

        write_file_path = os.path.join(self.raw_dir, f"{pdb}-clean.pdb")
        write_file = open(write_file_path, 'w')

        # Write the clean file, including HETATMs (if they are metal ions),
        # all atoms and lines starting with TER and END.
        test_MSE = False
        test_SEC = False
        test_KCX = False

        res_insertion_codes = []
        modres_sites = {}

        for line in read_file:

            # check for modified residues
            if line[:6] == 'MODRES':
                # only add sites which have changed residue code (keeping mod residues with same code the same)
                if line[12:15] != line[24:27]:
                    chain = line[16]
                    res_num = line[18:22].strip()
                    modres_sites[chain + res_num] = line[24:27]

            # check for resid insertion code in position 26
            if not self.include_hetatm:
                if line[:4] == 'ATOM' and len(line) > 26 and line[26].isalpha():
                    res_insertion_codes.append(line[26])
            else:
                if (line[:4] == 'ATOM' or line[:6] == 'HETATM') and len(line) > 26 and line[26].isalpha():
                    res_insertion_codes.append(line[26])

            # replace selenomethionine with methionine
            if "MSE" in line and ('ATOM' in line or 'HETATM' in line):
                line = line.replace('HETATM', 'ATOM  ')
                line = line.replace('MSE', 'MET')
                line = line.replace('SE', ' S')
                test_MSE = True

            # replace selenocysteine with cysteine
            if "SEC" in line and ('ATOM' in line or 'HETATM' in line):
                #line = line.replace("HETATM", "ATOM  ")
                line = line.replace('SEC', 'CYS')
                line = line.replace('SE', ' S')
                test_SEC = True

            #transform carboxylated lysine into a normal lysine
            if line[17:20] == 'KCX':

                if line[12:16].strip() in ['CX', 'OQ1', 'OQ2']:
                    continue

                else:
                    line = line.replace('HETATM', 'ATOM  ')
                    line = line.replace('KCX', 'LYS')

                test_KCX = True

            element_col = line[76:78].strip().upper()
            if not element_col:
                # guess element based on atom type name
                element_col = line[12:16].strip().lstrip('0123456789')[:2].upper()

            #check that amino acid is one of the 20 standard amino acids
            if (line[:4] == 'ATOM' or line[:6] == 'HETATM') and (line[21] + line[22:26].strip()) in list(list_prev_mod_resids):
                if line[12:16].strip() not in ['C', 'N', 'O', 'CA']:
                    continue
                line = line.replace('HETATM', 'ATOM  ')
                line = line.replace(line[17:20], list_prev_mod_resids[line[21] + line[22:26].strip()])
            elif (line[:4] == 'ATOM' or line[:6] == 'HETATM') and line[17:20].upper() not in standard_resids:

                if (line[21] + line[22:26].strip()) not in modres_sites:
                    continue
                if line[12:16].strip() not in ['C', 'N', 'O', 'CA']:
                    continue

                aa_num = int(line[22:26].strip())
                old_mod_aa_code = line[17:20]

                if uniprot_code != '':
                    if uniprot_fasta == '':
                        tmp_url = f'https://rest.uniprot.org/uniprotkb/{uniprot_code}.fasta'
                        fasta_text = requests.get(tmp_url, timeout=20).text
                        uniprot_fasta = ''.join(fasta_text.split('\n')[1:])
                    
                    M = bb.Molecule()
                    M.import_pdb(pdb=read_file_path, include_hetatm=True)
                    M = M.get_subset(M.atomselect('*', '*', 'CA', use_resname=True, get_index=True)[1])
                    subset_data = M.data.iloc[M.atomselect('*', '*', 'CA', use_resname=True, get_index=True)[1]]
                    pdb_seqs = {}
                    for chain in list(OrderedDict.fromkeys(subset_data['chain'])):
                        tmp_data = subset_data[subset_data['chain'] == chain]
                        pdb_seqs[chain] = ''.join([standard_resids[a] if a in list(standard_resids.keys()) else 'X' for a in list(tmp_data['resname'])])
                    chain_res_list = sorted(list(subset_data.loc[subset_data['chain'] == line[21], 'resid']))
                    aligner = PairwiseAligner()
                    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
                    aligner.open_gap_score = -11
                    aligner.extend_gap_score = -11
                    aligner.target_end_gap_score = 0.0
                    alignment = aligner.align(uniprot_fasta, pdb_seqs[line[21]])[0]

                    res_mapper = {}

                    for (uni_start, uni_end), (pdb_start, pdb_end) in zip(alignment.aligned[0], alignment.aligned[1]):
                        for uni_pos, pdb_pos in zip(range(uni_start+1, uni_end+1), range(pdb_start, pdb_end)):
                            res_mapper[chain_res_list[pdb_pos]] = uni_pos

                    aa_interest = uniprot_fasta[res_mapper[aa_num] - 1]
                    if aa_interest in list(standard_resids.values()):
                        line = line.replace('HETATM', 'ATOM  ')
                        line = line.replace(old_mod_aa_code, standard_resids_inv[aa_interest])
                        list_prev_mod_resids[line[21] + str(aa_num)] = standard_resids_inv[aa_interest]

                        # edit the fasta file to change the position
                        fasta_name = os.path.join(self.raw_dir, f'{pdb.split("-")[0]}.fasta')
                        if os.path.getsize(fasta_name) == 0:
                            self.download_fasta(pdb=pdb)
                        try:
                            f = open(fasta_name, 'r')
                            fasta_headers = []
                            fasta_chain_letters = []
                            fasta_seqs = []
                            for fasta_line in f:
                                if ">" in fasta_line:
                                    fasta_headers.append(fasta_line)
                                    chain_raw_info = fasta_line.split("|")[1][6:].split(",")
                                    if len(chain_raw_info[0]) == 1:
                                        chain_info = chain_raw_info
                                    else:
                                        chain_info = [a.strip()[0] for a in chain_raw_info]
                                    fasta_chain_letters.append(chain_info)
                                    if "fasta_seq" in locals():
                                        fasta_seqs.append(fasta_seq)

                                    fasta_seq = []

                                else:
                                    fasta_seq.append(fasta_line)
    
                            fasta_seqs.append(fasta_seq)
                            f.close()

                            new_fasta = open(fasta_name, 'w')
                            for header, chain, seq in zip(fasta_headers, fasta_chain_letters, fasta_seqs):
                                new_fasta.write(header)
                                seq = seq[0].strip()
                                if chain[0] == line[21]:
                                    
                                    F = bb.Molecule()
                                    F.import_pdb(pdb=read_file_path, include_hetatm=True)
                                    F = F.get_subset(F.atomselect('*', '*', 'CA', use_resname=True, get_index=True)[1])
                                    subset_data = F.data.iloc[F.atomselect('*', '*', 'CA', use_resname=True, get_index=True)[1]]
                                    pdb_seqs = {}
                                    for chain in list(OrderedDict.fromkeys(subset_data['chain'])):
                                        tmp_data = subset_data[subset_data['chain'] == chain]
                                        pdb_seqs[chain] = ''.join([standard_resids[a] if a in list(standard_resids.keys()) else 'X' for a in list(tmp_data['resname'])])
                                    chain_res_list = sorted(list(subset_data.loc[subset_data['chain'] == line[21], 'resid']))
                                    aligner = PairwiseAligner()
                                    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
                                    aligner.open_gap_score = -11
                                    aligner.extend_gap_score = -11
                                    aligner.target_end_gap_score = 0.0
                                    alignment = aligner.align(seq, pdb_seqs[line[21]])[0]

                                    for (seq_start, seq_end), (pdb_start, pdb_end) in zip(alignment.aligned[0], alignment.aligned[1]):
                                        fasta_mapper = dict(zip(range(pdb_start, pdb_end), range(seq_start, seq_end)))
                                    
                                    seq = seq[:(fasta_mapper[aa_num] - 1)] + aa_interest + seq[fasta_mapper[aa_num]:]
                                seq = seq + '\n'
                                new_fasta.write(seq)
                            new_fasta.close()

                        except Exception as e:
                            print(f'FASTA file parsing failed. Error: {str(e)}')
                            if 'f' in locals():
                                f.close()
                            fasta_chains = {}

                else:
                    raise Exception(f'Error in patching file: could not get residue type for modified '
                                    f'residue to convert to standard; file: {pdb}')

            #neglect HETATM atoms, unless they are metal ions
            if line[:6] == 'HETATM':
                if element_col in list_of_metals:
                    write_file.write(line)
                    continue

            # write lines and ignore hydrogen atoms
            if line[:4] == 'ATOM':
                if element_col[:1] not in ['H', 'D']:
                    write_file.write(line)
                    continue

            #same terminal statements
            if line[:3] in ['END', 'TER'] or line[0:6] == 'ENDMDL':
                write_file.write(line)
                continue

        if test_MSE: print(">> mutated MSE to MET")
        if test_SEC:
            print('>> mutated SEC to CYS')
            self._adapt_fasta(pdb=pdb)
        if test_KCX: print(">> mutation added to remove a lysine carboxylation")

        write_file.close()
        read_file.close()

        if res_insertion_codes:
            os.remove(write_file_path)
            os.remove(read_file_path)
            raise Exception(f'>> Residue insertion codes found for pdb {pdb}, '
                            f'insertion codes found: {res_insertion_codes}, residue numbering therefore not unique')

        os.remove(read_file_path)

        return replacement_dict


    def split_struc_nmr(self, pdb):
        '''
        Write a new file for each alternate NMR structure. Takes the clean pdb file produced
        previously and identifies statements in the pdb file indicating different conformations, for
        each conformation writes a new file of the format f'{pdb}-alt-{i}.pdb'.

        :param pdb: pdb name for the protein to split
        :type pdb: str
        '''

        #define the name format which will be followed for each file.
        number = 1
        name = f'{pdb}-alt-{str(number)}.pdb'

        path = os.path.join(self.raw_dir, f"{pdb}-clean.pdb")

        # Open the cleaned file, search over the data to find 'ENDMDL' statements to work out the number of structures available
        f = open(path)
        endmdl_instances = []

        for line in f:
            if re.search('ENDMDL', line):
                endmdl_instances.append(line)
        f.close()

        #If there are no 'ENDMDL' statements the clean file is renamed to suit the new format.
        if not endmdl_instances:
            path_rename = os.path.join(self.raw_dir, name)
            if os.path.exists(path_rename):
                os.remove(path_rename)
            os.rename(path, path_rename)

        #If there are 'ENDMDL' statements a new file is written for each model.
        elif endmdl_instances:

            print(">> Alternate model(s) found. Splitting...")

            #First it opens the clean file in the conformations folder and opens a new folder to write in
            f = open(path)

            path_rename = os.path.join(self.raw_dir, name)
            f_write = open(path_rename, 'w')

            #Next it writes the new file.
            #It includes every line until it gets to 'ENDMDL', where it opens a new file to write in.
            #The process stops when it gets to 'MASTER'
            for line in f:
                try:
                    line = str(line)

                    if re.search('^END\s*$', line):
                        f.close()
                        f_write.write(line)
                        f_write.close()
                        os.remove(path)
                        os.remove(path_rename)
                        break

                    elif re.search('ENDMDL', line):
                        f_write.write(line)
                        f_write.close()
                        number = number + 1
                        name = pdb + '-alt-' + str(number) + '.pdb'

                        path_rename = os.path.join(self.raw_dir, name)
                        f_write = open(path_rename, 'w')

                    else:
                        f_write.write(line)
                except Exception as e:
                    continue

            f_write.close()
            f.close()

        return


    def split_struc_alt_aa(self, pdb):
        '''
        Writes a new file for each alternative amino acid conformation present. Finds all pdb files
        created previously for the given pdb file, for each of the files, open and search column 16
        for any alternate structure codes present in the file. If there are codes present, for each
        code identified goes over the file and splits it up to produce files of the format:
        f'{pdb}-alt{i}{y}.pdb' where i is the nmr structure and y is the alternate code

        :param pdb: pdb code of the protein to curate
        :type pdb: str
        '''
        #Get a list of all the .pdb files present in [outdir]/conformations
        #and select those that belong to the pdb we are interested in.
        list_of_files = glob.glob(os.path.join(self.raw_dir, f"*{pdb}*.pdb"))

        for f in list_of_files:

            #Next it checks if there are any alternate amino acid conformations present (i.e. if line[16 == A, B or C]).
            read_file = open(f, 'r')
            alt_loc_vals = sorted({line[16] for line in read_file if line[:4] == 'ATOM' and line[16] != ' '})
            read_file.close()

            # If there are none for this file it moves on to the next.
            if not alt_loc_vals:
                continue
            else:
                print(">> Alternate amino acid conformations found. Splitting...")

            # If there are it rewrites a new file for all alternate conformations
            for i, alt_loc_letter in enumerate(alt_loc_vals):
                try:
                    prev_name, ext = os.path.splitext(f)
                    pdb_name, _, ver = prev_name.rpartition('-')
                    write_path = f'{pdb_name}{ver}{alt_loc_letter}{ext}'
                    f_write = open(write_path, 'w')

                    target_letter = alt_loc_letter
                    non_target_letters = []
                    for letter in alt_loc_vals:
                        if letter != target_letter:
                            non_target_letters.append(letter)

                    read = open(f)
                    for line in read:

                        if (line[:4] == 'ATOM' or line[:6] == 'HETATM') and (line[16] == target_letter):
                            newline = line[:16] + ' ' + line[17:]
                            f_write.write(newline)
                            continue

                        if (line[:4] == 'ATOM' or line[:6] == 'HETATM') and (line[16] in non_target_letters):
                            continue

                        if (line[:4] == 'ATOM') or (line[:6] == 'HETATM') or (line[:3] == 'TER'):
                            f_write.write(line)

                    f_write.close()
                    read.close()

                except Exception as e:

                    if "read_file" in locals():
                        read_file.close()

                    if "f_write" in locals():
                        f_write.close()

                    raise Exception(f'{str(e)}')

            #The original is then removed if it has been replaced.
            os.remove(f)

        return


    def _align_resnum_uniprot(self, uniprot_code, pdb_code, chains):
        '''
        Function to align the residue numbers within the pdb file with the canonical sequence
        available from the Uniprot website. This will perform an alignment using the biopython
        pairwise aligner with BLOSUM62 matrix, biased towards matching at the start of the uniprot
        sequence. Sequence is derived for the pdb using biobox. After alignment performs a mapping
        on the biobox instance and rewrites the pdb with updated chain numbering.

        :param uniprot_code: Uniprot code for pdb file of interest
        :type uniprot_code: str
        :param pdb_code: Code for the pdb structure of interest
        :type pdb_code: str
        :param chains: List of the chains within the pdb structure that match the uniprot code
            within the pdb file
        :type chains: list
        '''
        tmp_url = f'https://rest.uniprot.org/uniprotkb/{uniprot_code}.fasta'
        fasta_text = requests.get(tmp_url, timeout=20).text
        uniprot_fasta = ''.join(fasta_text.split('\n')[1:])
        M = bb.Molecule()
        M.import_pdb(pdb_code, include_hetatm=self.include_hetatm)
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
                                'KCX': 'K', 'LYE': 'K', 'LYSN': 'K'}

        pdb_seqs = {}
        for chain in list(OrderedDict.fromkeys(subset_data['chain'])):
            tmp_data = subset_data[subset_data['chain'] == chain]
            pdb_seqs[chain] = ''.join([protein_letters_dict[a] if a in list(protein_letters_dict.keys()) else 'X' for a in list(tmp_data['resname'])])

        if isinstance(chains, str):
            chains = [c for c in chains.split('/') if c]

        failed = []
        for chain in chains:
            try:
                aligner = PairwiseAligner()
                aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
                aligner.open_gap_score = -11
                aligner.extend_gap_score = -11
                aligner.target_end_gap_score = 0.0
                alignment = aligner.align(uniprot_fasta, pdb_seqs[chain])[0]

                res_mapper = {}
                chain_res_list = sorted(list(subset_data.loc[subset_data['chain'] == chain, 'resid']))

                for (uni_start, uni_end), (pdb_start, pdb_end) in zip(alignment.aligned[0], alignment.aligned[1]):
                    for uni_pos, pdb_pos in zip(range(uni_start+1, uni_end+1), range(pdb_start, pdb_end)):
                        res_mapper[chain_res_list[pdb_pos]] = uni_pos

                mapped_res = M.data['resid'].map(res_mapper)
                sel_chain = M.data['chain'] == chain
                if mapped_res[sel_chain].isna().any():
                    failed.append(chain)
                    continue
                M.data.loc[sel_chain, 'resid'] = mapped_res[sel_chain].astype(M.data['resid'].dtype)

            except Exception as e:
                print(f'Failed alignment of pdb {pdb_code}, chain {chain}, with error: {str(e)}')
                failed.append(chain)

        if not failed:
            M.write_pdb(pdb_code)
            print(f'>> Chains aligned to canonical uniprot sequence for pdb code: {pdb_code}')

        else:
            print(f'>> Alignment failed for pdb: {pdb_code}, not writing new file as not all chains matched properly')


    def get_chain_replacement(self, pdb_code):
        '''
        Identify the names of the chains given from the FASTA file that has been downloaded for the
        specific PDB code.

        .. rubric:: Method

        Download and open the corresponding FASTA file for the PDB code and parse the data to
        identify the chain names for the protein. If a chain has two different names (i.e. an auth
        name and the name given by the RCSB) the chain name will have the following format in the
        fasta file: Chains K[auth M], L[auth N]. Therefore the code here puts the auth name and name
        given by RCSB into a dictionary which is later used to replace the auth names. e.g. the two
        examples here would be appended into the dictionary as {M:K, N:L}. If the auth name and RCSB
        name are the same nothing is appended to the dictionary.

        :param pdb_code: PDB code of the file to get the chain information for
        :type pdb_code: str
        :returns: Dictionary of the mapping from auth name to rcsb chain name for all chains in
            structure
        :rtype: dict
        '''
        try:
            self.download_fasta(pdb_code)
        except Exception as e:
            raise Exception(f"Could not download FASTA file. Error: {str(e)}")

        name = os.path.join(self.raw_dir, pdb_code + '.fasta')

        try:
            #append chain information into the chains_raw list
            f = open(name, 'r')
            list_of_chains = []
            chains_raw =[]
            replacement = []
            need_replacing = []
            replacement_dict = dict()
            for line in f:
                m = re.findall(r'Chain[ a-z , 0-9, A-Z \[\]]*|', line)
                for entry in m:
                    if len(entry) > 0:
                        chains_raw.append(entry)

            f.close()

        except Exception as e:
            print(f'FASTA file parsing failed. Error: {str(e)}')
            f.close()
            return {}

        for entry in chains_raw:
            try:
                chains_raw_2 = entry.split(',')
                for entry in chains_raw_2:
                    if re.search('auth', entry):
                        chains_wrong = entry.split('auth')
                        for entry in chains_wrong:
                            select = entry.split(' ')
                            for entry in select:
                                if entry[-1:] == '[':
                                    replacement = entry[:-1]

                                if entry[-1:] == ']':
                                    need_replacing = entry[:-1]
                        replacement_dict[need_replacing] = replacement

                    else:
                        chain = entry.split(' ')
                        for entry in chain:
                            if 'Chain' not in entry and (len(entry) > 0):
                                list_of_chains.append(entry)

            except Exception as e:
                print(f'Failed renaming chains. Error: {str(e)}')
                continue

        #dictionary used to replace the auth chain names with the RCSB chain names
        return replacement_dict


    def _check_curated_structure(self, pdb, uniprot='', chains=[]):
        '''
        Call the required functions that may be needed if a curated pdb file is already found.
        Potential for this to be needed if a pdb is curated for a uniprot code but pdb file also
        contains chains from another uniprot code which may not have been corrected fully. This is
        non-exhaustive and as other checks are added, may need to be included here.

        :param pdb: file name for the pdb file to check
        :type pdb: str
        :param uniprot: uniprot code for the file to check. Used to perform the alignment of the
            structure to the uniprot sequence
        :type uniprot: str
        :param chains: list of chains within the structure to check over
        :type chains: list
        '''
        self._align_resnum_uniprot(uniprot, pdb, chains)


    def replace_chains(self, path, replacement_dict):
        '''
        Replace chain names in a structure. Using the replacement dict given as a parameter, it
        loads in a biobox instance of the pdb file and replaces the chain names in the 'chain'
        column. It then writes a temporary pdb file before inserting the new data into the old file
        to keep the additional lines present.

        :param path: The path of the pdb file that the work is being done on
        :type path: str
        :param replacement_dict: The dictionary which contains the information about which chains
            need replacing
        :type replacement_dict: dict
        '''

        #The relevant file in conformations is then opened and rewritten.
        try:
            M = bb.Molecule()
            M.import_pdb(path, include_hetatm=self.include_hetatm)
            # find the indices of the atoms in the pdb file
            indices = M.atomselect('*', '*', '*', True, False)[1]

            M.data['chain'] = M.data['chain'].replace(replacement_dict)

            pdb = path.split(os.sep)[-1].split('.')[0]
            path_temp = os.path.join(self.raw_dir, f"{pdb}_temp.pdb")
            M.write_pdb(path_temp, index=indices, split_struc=False)

            # take the new written file and insert in place where it would sit in the overall pdb file
            lines = open(path, 'r').readlines()
            first_lines = []
            for line in lines:
                if line[:4]  == "ATOM" or line[:6] == 'HETATM' or line[:3] == 'TER' or line[:5] == 'MODEL':
                    break
                first_lines.append(line)

            new_atom_lines = open(path_temp, 'r').readlines()
            new_file_output = first_lines + new_atom_lines + ['END\n']

            cleaned_file = open(path, 'w')
            cleaned_file.writelines(new_file_output)
            cleaned_file.close()

            os.remove(path_temp)


        #If the protein fails the replacement process it is removed from the conformations folder.
        except Exception as e:
            os.remove(path)
            raise Exception(f'Failed replacing chains. Error: {str(e)}')

        return

    def _adapt_fasta(self, pdb):
        '''
        Short function called if selenocysteine found in the structure to adapt the fasta file in
        order for modeller to be able to be called to fix the structure. This writes over the fasta
        file with the correct sequence.

        :param pdb: pdb code of the fasta file to change the data
        :type pdb: str
        '''
        try:
            fasta = os.path.join(self.raw_dir, f"{pdb}.fasta")
            fin = open(fasta, "r")
            fasta_headers = []
            fasta_seqs = []
            for line in fin:
                if ">" in line:
                    fasta_headers.append(line)
                else:
                    temp_sequence = line
                    for i, res in enumerate(temp_sequence):
                        if res == 'U':
                            temp_sequence = temp_sequence[:i] + 'C' + temp_sequence[(i+1):]
                    fasta_seqs.append(temp_sequence)

            new_fasta = open(fasta, 'w')
            for header, sequence in zip(fasta_headers, fasta_seqs):
                new_fasta.write(header)
                new_fasta.write(sequence)
            new_fasta.close()

        except Exception as e:
            print(f'Failed to replace selenocysteine in fasta file, error: {str(e)}')
            raise Exception(f'Failed to replace selenocysteine in fasta file, error: {str(e)}')


    def rewrite_pdb(self, path):
        '''
        Short function to take a pdb file, load it into biobox as a molecule and rewrite a pdb file
        from the molecule class. If it fails it remove the path of the file to rewrite.

        :param path: the path to the pdb file to rewrite
        :type path: str

        .. rubric:: Example

        ::

            >>> self.rewrite_pdb(path)
        '''
        try:
            M = bb.Molecule()
            M.import_pdb(path, include_hetatm=self.include_hetatm)
            indices = M.atomselect('*', '*', '*', True, False)[1]
            M.write_pdb(path, index=indices, split_struc=False)

        except Exception as e:
            os.remove(path)
            print(f'Failed rewriting pdb file with error: {str(e)}')



if __name__ == "__main__":

    PDB = PDB(outdir='result')
    #PDB.clean_and_split_pdb('13LD', 'P10724') # test KCX to LYS mutation
    #PDB.clean_and_split_pdb('1PAE', 'P22887') # test SEC to CYS mutation
    #PDB.clean_and_split_pdb('6XZ7', 'P60422') # test MSE to MET mutation
    #PDB.clean_and_split_pdb('2MBH', 'Q13351') # test splitting of models
    #PDB.clean_and_split_pdb('1A6M', 'P02185') # test splitting rotamers
    #PDB.clean_and_split_pdb('4WNC', 'P04406') # test splitting rotamers
    #PDB.clean_and_split_pdb('3DBJ', 'P50030', chains=['A', 'C', 'E', 'G']) # test renumbering residues with canonical uniprot sequence
    PDB.clean_and_split_pdb('2MWS', 'P0CG48') #  test removal of modified residue
