#OVERALL STRUCTURE:
#- download file into [outdir]/conformations folder
#- clean files (i.e. removes heteroatoms that aren't metal ions)
#- write files for alternate conformations (usually from NMR ensembles)
#- write files for alternate amino acid conformations

import os
import subprocess
import re
import glob
import fileinput
import requests
import Bio
import pandas as pd
import numpy as np
import biobox as bb
import alphafold as af # to load alphafold data
import patcher # to patch PDB structures with missing regions
from helper import get_download_tool, ShutUp
from collections import OrderedDict
from Bio.Align import PairwiseAligner, substitution_matrices


class PDB(object):
    '''
    Class for taking an input dataframe of the proteins required to investigate, downloading
    the proteins required from both the PDB and alphafold and cleaning the structures up
    to ensure that these are all ready to be used within measurements for the prediction
    model input. Cleaning is currently done using Modeller.
    '''

    def __init__(self, outdir="result", gap=10, PDB_only=False):
        '''
        Initialise the PDB class.
        
        Parameters
        ----------
        outdir : string
            The directory in which files should be downloaded and curated within
        
        gap : int
            The maximum gap that is allowed in the sequence for a structure that has been
            downloaded that patching will be done on. For structures with a gap in the
            sequence greater than this, the structure will be removed.
        
        PDB_only : bool
            Toggle for if you want to download a system from a list of PDB files (True)
            or from a Uniprot dataframe (False) created from the Uniprot class. 
        '''
        self.outdir = outdir
        self.PDB_only = PDB_only

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
            columns = ['Uniprot_Entry', 'PDB_Code', 'Method', 'Resolution', 'Chains']
            self.df = pd.DataFrame(columns=columns)

        #gap to consider as small enough to justify patching
        self.gap = gap


    def save_state(self, outname="proteins.csv"):
        '''
        Take the current dataframe of proteins that have been gathered (pdb.df) and write this
        to a csv file in the output directory for later use.

        Parameters
        ----------
        outname -> string
            the desired name for the pdb.df csv file to be written as

        Example
        -------
        pdb.save_state()
        '''
        self.df.to_csv(os.path.join(self.outdir, outname), index_label=False, index=False)


    def load_state(self, fname, outdir="", gap=10, PDB_only=False):
        '''
        initialize DataFrame from csv file.
        If no output directory is given, the folder containing the csv file is used.
        '''
        # TODO GW 30.07.26 - could we infer the PDB_only state from the measures df format
        if outdir == "":
            outdir = os.path.dirname(fname)

        self.PDB_only = PDB_only
        self.outdir = outdir
        self.gap = gap

        try:
            self.df = pd.read_csv(fname)
            # remove duplicates from the dataframe to avoid extra unneccessary calculations
            self.df = self.df.drop_duplicates()
            
        except Exception as e:
            print(f"Could not load csv file, error: {str(e)}")


    def gather_proteins(self, uniprot_df, skip_if_found=True, gap=10):
        '''
        Iterate over lines of a DataFrame containing PDB structure information.
        Download every stucture, and curate it if necessary
        '''
        if self.PDB_only:
            try:
                # when the inputs are PDB codes only, convert them to a dataframe
                dic = {'PDB_Code':uniprot_df}
                uniprot_df = pd.DataFrame(dic)
            except Exception as e:
                return e

        for idx, row in uniprot_df.iterrows():

            pdb_code = row["PDB_Code"]

            if not self.PDB_only:
                uniprot_code = row["Uniprot_Entry"]
                print(f"\nUNIPROT: {uniprot_code}, PDB: {pdb_code}")
            else:
                print(f"\nPDB: {pdb_code}")

            if pdb_code[:2] == "AF":

                if skip_if_found:
                    files=[os.path.basename(c).split(".")[0] for c in glob.glob(os.path.join(self.curated_dir, "*pdb"))]
                    if pdb_code in files:
                        print(f">> curated {pdb_code} PDB found, continuing...")
                        if not self.PDB_only:
                            data = {'Uniprot_Entry': uniprot_code, 'PDB_Code': pdb_code, 'Method': 'Predicted', 'Resolution': np.nan, 'Chains': "A"}
                        else:
                            data = {'PDB_Code': pdb_code}

                        self.df = pd.concat([self.df, pd.DataFrame.from_records(data, index=[0])], ignore_index=True)
                        continue

                try:
                    af.download_AF_struc(pdb_code,outfolder=self.outdir)
                except Exception as e:
                    print(f">> FAILED on calling download_AF_struc: {str(e)}")
                    continue

                # check if the AlphaFold file contains ATOM statements
                af_filename = os.path.join(self.curated_dir, f"{pdb_code}.pdb")
                fin = open(af_filename, "r")
                test = False
                for line in fin:
                    if line.split()[0] == "ATOM":
                        test = True
                        break
                fin.close()

                if test:
                    if not self.PDB_only:
                        data = {'Uniprot_Entry': uniprot_code, 'PDB_Code': pdb_code, 'Method': 'Predicted', 'Resolution': np.nan, 'Chains': "A"}
                    else:
                        data = {'PDB_Code': pdb_code}
                    self.df = pd.concat([self.df, pd.DataFrame.from_records(data, index=[0])], ignore_index=True)

                    # find the PLDDT codes for AF structures - not fully sure where to put this
                    af.find_af_plddt(pdb_code,outfolder=self.outdir)

                else:
                    print(">> FAILED: structure not found in AlphaFold database")
                    try:
                        os.remove(af_filename)
                    except Exception as e:
                        pass

            else:
                try:
                    if not self.PDB_only:
                        method_obtained = row["Method"]
                        resolution = row["Resolution"]
                        chains = row["Chains"]
                        if isinstance(chains, str):
                            chains = [c for c in chains.split('/') if c]
                except Exception as e:
                    method_obtained, resolution, chains = '', '', ''
                    pass
                if skip_if_found:
                    files=[os.path.basename(c).split("-")[0] for c in glob.glob(os.path.join(self.curated_dir, "*pdb"))]
                    if pdb_code in files:
                        print(f">> curated {pdb_code} PDB found, continuing...")
                        if not self.PDB_only:
                            matched_curated_files = [os.path.basename(a) for a in glob.glob(os.path.join(self.curated_dir, "*pdb")) if pdb_code.upper() == os.path.splitext(os.path.basename(a))[0].split('-')[0]]
                            for file in matched_curated_files:
                                self._check_curated_structure(os.path.join(self.curated_dir, file), uniprot_code, chains)
                            data = {'Uniprot_Entry': uniprot_code, 'PDB_Code': pdb_code, 'Method': method_obtained, 'Resolution': resolution, 'Chains': '/'.join(chains)}
                        else:
                            data = {'PDB_Code': pdb_code}
                        self.df = pd.concat([self.df, pd.DataFrame.from_records(data, index=[0])], ignore_index=True)
                        continue

                # load, clean, and split it in alternate conformations
                try:
                    if not self.PDB_only:
                        self.clean_and_split_pdb(pdb_code, uniprot_code, chains)
                        data = {'Uniprot_Entry': uniprot_code, 'PDB_Code': pdb_code, 'Method': method_obtained, 'Resolution': resolution, 'Chains': '/'.join(chains)}
                    else:
                        self.clean_and_split_pdb(pdb_code)
                        data = {'PDB_Code': pdb_code}
                    self.df = pd.concat([self.df, pd.DataFrame.from_records(data, index=[0])], ignore_index=True)

                except Exception as e:
                    print(f">> FAILED clean and splitting of PDB entry file: {str(e)}")
                    continue


    def clean_and_split_pdb(self, pdb, uniprot_code = '', chains=[]):
        '''
        Download a pdb, and return a collection of cleaned and splitted alternative conformations
        The results are saved into files: [outfolder]/conformations/*PDB code*-clean.pdb.
        '''
        try:
            #download and clean the structure
            self.download_pdb(pdb)
            replacement_dict = self.clean(pdb)

            #splits into all alternative conformations into independent structures
            self.split_struc_nmr(pdb)
            self.split_struc_alt_aa(pdb)

        except Exception as e:
            raise Exception(f'Error cleaning {pdb}: {str(e)}') from e


        files = glob.glob(os.path.join(self.raw_dir, f"*{pdb}*pdb"))
        test = False
        for cnt, f in enumerate(files):
            mypath = os.path.split(f)[0]
            fasta_loc = os.path.join(mypath, f"{pdb}.fasta")
            print(fasta_loc)

            try:

                fname = patcher.curate(f, fasta_loc, outdir=self.curated_dir, gap=self.gap)
                if len(replacement_dict) > 0:
                    reverse_replacement_dict = dict((v,k) for k,v in replacement_dict.items())
                    #self.replace_chains(fname, reverse_replacement_dict)
                    self.new_replace_chains(fname, reverse_replacement_dict)

                # TODO 24.03.26 - this is only possible with structures that have a uniprot code associated in order to get the sequence to align to - find way to make work with pdb_only
                if not self.PDB_only:
                    self._align_resnum_uniprot(uniprot_code, fname, chains)

                test = True

            except Exception as e:
                print(f">> Patching failed for conformer {cnt}. {str(e)}")
                continue

        if not test:
            raise Exception("Patching failed for all conformers")

        return


    def download_pdb(self, pdb):
        '''
        Download the required PDB file for the PDB code specified. This is used for the coordinates
        of the protein to extract the featurised data for each protein from.
        
        Method
        ------
        Identify the download tool that is available for use. Check if PDB file has already been
        downloaded, if so skip. Otherwise download the PDB file from the RCSB website. 
        
        Parameters
        ----------
        pdb -> string
            the PDB code for the PDB file to be downloaded
            
        Example
        -------
        >> self.download_pdb('1ubq')
        '''
        files=glob.glob(os.path.join(self.raw_dir, "*pdb"))
        test_file = os.path.join(self.raw_dir, f'{pdb}.pdb')
        if test_file not in files:
            print(f">> Downloading PDB {pdb}")
            tool = get_download_tool()
            try:
                if tool == "curl":
                    line = f"curl -s -f -o {pdb}.pdb https://files.rcsb.org/download/{pdb}.pdb"
                elif tool == "wget":
                    line = f"wget https://files.rcsb.org/download/{pdb}.pdb"

                else:
                    raise RuntimeError("You don't have a commandline tool for downloading files")

                response = requests.get(url=f'https://files.rcsb.org/download/{pdb}.pdb', timeout=20)
                response.raise_for_status()
                subprocess.check_call(line, shell=True)
                os.rename(os.path.join(os.getcwd(), f'{pdb}.pdb'), os.path.join(self.raw_dir, f'{pdb}.pdb'))

            except Exception as e:
                print(f'Error downloading file. {str(e)}')
        else:
            print(f'PDB file ({pdb}) has previously been downloaded, using previous copy.')


    def download_fasta(self, pdb):
        '''
        Download the required fasta file for the PDB code specified. This is used to
        align the structures and check for missing residues and extract sequences for comparison.
        
        Method
        ------
        Identify the download tool that is available for use. Check if fasta file has already been
        downloaded, if so skip. Otherwise download the fasta file from the RCSB website. 
        
        Parameters
        ----------
        pdb -> string
            the PDB code for the fasta file to be downloaded
            
        Example
        -------
        >> self.download_fasta('1ubq')
        '''
        tool = get_download_tool()
        files=[c.split(os.sep)[-1][:4] for c in glob.glob(os.path.join(self.raw_dir, "*.fasta"))]
        if pdb not in files:
            try:
                print(f">> downloading FASTA for {pdb}")

                if tool == "curl":
                    line = f"curl -s -f -o {pdb}.fasta https://www.rcsb.org/fasta/entry/{pdb}/download"
                elif tool == "wget":
                    line = f"wget -O {pdb}.fasta https://www.rcsb.org/fasta/entry/{pdb}/download"
                else:
                    raise RuntimeError("You don't have a commandline tool for downloading files")

                response = requests.get(url=f'https://www.rcsb.org/fasta/entry/{pdb}/download', timeout=20)
                response.raise_for_status()
                subprocess.check_call(line, shell=True)
                os.rename(os.path.join(os.getcwd(), f'{pdb}.fasta'), os.path.join(self.raw_dir, f'{pdb}.fasta'))

            except Exception as e:
                print(f'>> Failed downloading FASTA sequence for chain name comparison for {pdb}, error: {str(e)}')
                raise Exception(f'Failed downloading FASTA sequence for chain name comparison.{str(e)}') from e
        else:
            print(f'Fasta file for {pdb} previously downloaded, using previous copy.')


    def clean(self, pdb):
        '''
        Rename the protein's chains during the cleaning process so that they match
        the chain names given in the FASTA file.
        This is required as pdb files name their chains using the 'auth' name and
        fasta with the normal chain name.
        Therefore to avoid confusion we rename them all to what is used in the fasta file.
        '''
        try:
            replacement_dict = self.get_chain_replacement(pdb)
            path = os.path.join(self.raw_dir, f"{pdb}.pdb")
            if len(replacement_dict) > 0:
                self.new_replace_chains(path, replacement_dict)
                #self.replace_chains(path, replacement_dict)

        except Exception as e:
            raise Exception(f'Error renaming chains. {str(e)}')

        #Next it opens and starts reading the .pdb file and starts writing a new file with the ending '-clean.pdb'.
        list_of_metals = ['ZN', 'NI', 'CU', 'FE', 'MG', 'MN', 'NA', 'K', 'CA', 'CO', 'CL', 'MO']

        read_file_path = os.path.join(self.raw_dir, f"{pdb}.pdb")
        read_file = open(read_file_path)

        write_file_path = os.path.join(self.raw_dir, f"{pdb}-clean.pdb")
        write_file = open(write_file_path, 'w')

        # Write the clean file, including HETATMs (if they are metal ions),
        # all atoms and lines starting with TER and END.
        test_MSE = False
        test_SEC = False
        test_KCX = False
        for line in read_file:

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

            #neglect HETATM atoms, unless they are metal ions
            if line[:6] == 'HETATM':
                if line[76:78].strip().upper() in list_of_metals:
                    write_file.write(line)
                    continue

            # write lines and ignore hydrogen atoms
            if line[:4] == 'ATOM':
                if line[76:78].strip() != 'H':
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

        #self.rewrite_pdb(write_file_path)
        #remove the original pdb file as it is not needed anymore.
        os.remove(read_file_path)

        return replacement_dict


    def split_struc_nmr(self, pdb):
        '''
        Write a new file for each alternate NMR structure.
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
        Writes a new file for each alternative amino acid conformation present.
        '''
        #Get a list of all the .pdb files present in [outdir]/conformations
        #and select those that belong to the pdb we are interested in.
        list_of_files = glob.glob(os.path.join(self.raw_dir, f"*{pdb}*.pdb"))

        for f in list_of_files:

            #Next it checks if there are any alternate amino acid conformations present (i.e. if line[16 == A, B or C]).
            read_file = open(f, 'r')
            alt_loc_vals = sorted({line[16] for line in read_file if line[:4] == 'ATOM' and line[16] != ' '})

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

                        if (line[:4] == 'ATOM') and (line[16] == target_letter):
                            newline = line[:16] + ' ' + line[17:]
                            f_write.write(newline)
                            continue

                        if (line[:4] == 'ATOM') and (line[16] in non_target_letters):
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
            read_file.close()
            os.remove(f)

        return


    def _align_resnum_uniprot(self, uniprot_code, pdb_code, chains):
        '''
        Function to align the residue numbers within the pdb file with the canonical
        sequence available from the Uniprot website

        Parameters
        ----------
        uniprot_code -> string
            Uniprot code for pdb file of interest
        pdb_code -> string
            Code for the pdb structure of interest
        chains -> list
            List of the chains within the pdb structure that match the uniprot code
            within the pdb file
        '''
        tmp_url = f'https://rest.uniprot.org/uniprotkb/{uniprot_code}.fasta'
        fasta_text = requests.get(tmp_url).text
        uniprot_fasta = ''.join(fasta_text.split('\n')[1:])
        M = bb.Molecule(pdb_code)
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
            align_shift_dict = {}
            try:
                aligner = PairwiseAligner()
                aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
                aligner.open_gap_score = -11
                aligner.extend_gap_score = -1
                aligner.target_end_gap_score = 0.0
                alignment = aligner.align(uniprot_fasta, pdb_seqs[chain])[0]
                dash_locs = [i for i, aa in enumerate(alignment[1]) if aa =='-']
                for idx, loc in enumerate(dash_locs):
                    align_shift_dict[loc - (idx)] = idx + 1
                old_res_count = 1
                old_res_curr_num = -1
                curr_chain = 'XXXX'
                for i, r in M.data.iterrows():
                    if r['chain'] == chain:
                        if old_res_curr_num == -1: old_res_curr_num = r['resid']
                        if curr_chain == 'XXXX': curr_chain = r['chain']
                        shift = 0
                        res_num = r['resid']
                        if len(list(align_shift_dict.keys())) != 0:
                            for bound in list(align_shift_dict.keys()):
                                if res_num > bound:
                                    shift = align_shift_dict[bound]
                        if old_res_curr_num != r['resid']:
                            old_res_count += 1
                            old_res_curr_num = r['resid']
                        if curr_chain != r['chain']:
                            curr_chain = r['chain']
                            old_res_count = 1
                        M.data.at[i, 'resid'] = old_res_count + shift
            except Exception as e:
                print(f'Failed alignment of pdb {pdb_code}, chain {chain}, with error: {str(e)}')
                failed.append(chain)

        if failed:
            #raise RuntimeError(f'>> Alignment failed for pdb: {pdb_code}, not writing new file as not all chains matched properly')
            print(f'>> Alignment failed for pdb: {pdb_code}, not writing new file as not all chains matched properly')
        
        else:
            M.write_pdb(pdb_code)
            print(f'Chains aligned to canonical uniprot sequence for pdb code: {pdb_code}')


    def get_chain_replacement(self, pdb_code):
        '''
        Identify the names of the chains given from the FASTA file that has been
        downloaded for the specific PDB code. 

        Method
        ------
        Open the corresponding FASTA file for the PDB code and parse the data to
        identify the chain names for the protein. 
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

            #The fasta file is then removed
            f.close()
            
        except Exception as e:
            print(f'FASTA file parsing failed. Error: {str(e)}')
            f.close()
            #os.remove(name)
            return {}

        #If a chain has two different names (i.e. an auth name and the name given by the RCSB)
        # the chain name will have the following format in the fasta file:
        #Chains K[auth M], L[auth N]
        #Therefore the code here puts the auth name and name given by RCSB into a dictionary which is later used to replace the auth names.
        #e.g. the two examples here would be appended into the dictionary as {M:K, N:L}
        #If the auth name and RCSB name are the same nothing is appended to the dictionary.
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


    def replace_chains(self, path, replacement_dict):
        '''
        Old replace chains, left here incase useful when implementing the new_replace chains method
        This now isnt called, left incase problems arise with new method 19.05.26
        '''
        #First, the auth chain names are put in a list.
        auth_list = list(replacement_dict)

        #The relevant file in conformations is then opened and rewritten.
        try:
            with fileinput.FileInput(path, inplace = True) as f:
                for line in f:
                    try:
                        #On lines with 'ATOM', 'TER' or 'HETATM' if the chain is in the auth_list
                        #the auth chain name is replaced with the RCSB chain name.
                        if (line[:4] == 'ATOM') or (line[:3] == 'TER') or (line[:6] == 'HETATM'):
                            chain_name = line[21]
                            if chain_name in auth_list:
                                replacement_chain_name = replacement_dict.get(chain_name)
                                if line[25] != ' ':
                                    line = line[:21] + replacement_chain_name + line[22:]
                                else:
                                    line = line[:21] + replacement_chain_name + line[22:]

                            while line[31] == ' ':
                                line = line[:31] + line[32:]

                            print(line, end ='')

                        else:
                            print(line, end='')
                    except:
                        print(line, end='')

        #If the protein fails the replacement process it is removed from the conformations folder.
        except Exception as e:
            os.remove(path)
            raise Exception(f'Failed replacing chains. Error: {str(e)}')

        return


    def _check_curated_structure(self, pdb, uniprot='', chains=[]):
        '''
        Call the required functions that may be needed if a curated pdb file is already found.
        Potential for this to be needed if a pdb is curated for a uniprot code but pdb file also
        contains chains from another uniprot code which may not have been corrected fully. This
        is non-exhaustive and as other checks are added, may need to be included here.
        '''
        self._align_resnum_uniprot(uniprot, pdb, chains)


    def new_replace_chains(self, path, replacement_dict):
        '''
        New version of the replace chains which correctly produces pdb files afterwards
        for the patching as the original left some proteins in space. 

        Parameters
        ----------
        path -> string
            The path of the pdb file that the work is being done on

        replacement_dict -> dict
            The dictionary which contains the information about which chains
            need replacing 

        '''

        # First, the auth chain names are put in a list.
        auth_list = list(replacement_dict)

        #The relevant file in conformations is then opened and rewritten.
        try:
            M = bb.Molecule(path)
            # find the indices of the atoms in the pdb file
            indices = M.atomselect('*', '*', '*', True, False)[1]

            # transform to data lists and change the chain names
            for index, row in M.data.iterrows():
                if row['chain'] in auth_list:
                    replacement_chain_name = replacement_dict.get(row['chain'])
                    M.data.at[index, 'chain'] = replacement_chain_name

            pdb = path.split(os.sep)[-1].split('.')[0]
            path_temp = os.path.join(self.raw_dir, f"{pdb}_temp.pdb")
            M.write_pdb(path_temp, index=indices, split_struc=True)

            # take the new written file and insert in place where it would sit in the overall pdb file
            lines = open(path, 'r').readlines()
            start_atoms = False
            first_lines = []
            second_lines = []
            for line in lines:
                if line[:4]  == "ATOM" or line[:6] == 'HETATM' or line[:3] == 'TER':
                    start_atoms = True
                else:
                    if start_atoms:
                        second_lines.append(line)
                    else:
                        first_lines.append(line)

            new_atom_lines = open(path_temp, 'r').readlines()
            new_file_output = first_lines + new_atom_lines[1:-1] + second_lines

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
        Short function called if selenocysteine found in the structure to adapt the
        fasta file in order for modeller to be able to be called to fix the structure.

        Parameters
        ----------
        pdb -> string
            pdb code of the fasta file to change the data
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
        Short function to take a pdb file, load it into biobox as a molecule and
        rewrite a pdb file from the molecule class.

        Parameters
        ----------
        path -> string
            the path to the pdb file to rewrite

        Example
        -------
        >> self.rewrite_pdb(path)
        '''
        try:
            M = bb.Molecule(path)
            indices = M.atomselect('*', '*', '*', True, False)[1]
            M.write_pdb(path, index=indices, split_struc=False)

        except Exception as e:
            os.remove(path)
            print(f'Failed rewriting pdb file with error: {str(e)}')



if __name__ == "__main__":

    PDB = PDB(outdir='result')
    #PDB.clean_and_split_pdb('13LD') # test KCX to LYS mutation
    #PDB.clean_and_split_pdb('1PAE') # test SEC to CYS mutation
    #PDB.clean_and_split_pdb('6XZ7') # test MSE to MET mutation
    #PDB.clean_and_split_pdb('2MBH') # test splitting of models
    PDB.clean_and_split_pdb('1A6M') # test splitting rotamers
    #PDB.clean_and_split_pdb('4WNC', 'P04406') # test splitting rotamers
    #PDB.clean_and_split_pdb('3DBJ', 'P50030', chains=['A', 'C', 'E', 'G']) # test renumbering residues with canonical uniprot sequence


    #PDB._align_resnum_uniprot('P50030', f'result{os.sep}curated{os.sep}3DBJ-alt-1.pdb', chains=['A', 'C', 'E', 'G'])
