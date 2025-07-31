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
import shutil
import pandas as pd
import numpy as np
import biobox as bb
import alphafold as af # to load alphafold data
import patcher # to patch PDB structures with missing regions
from helper import get_download_tool, ShutUp


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

        if outdir == "":
            outdir = os.path.dirname(fname)

        self._setup(outdir, gap, PDB_only)

        try:
            self.df = pd.read_csv(fname)
            # remove duplicates from the dataframe to avoid extra uneccessary calculations
            self.df = self.df.drop_duplicates()
        except Exception as e:
            print(f"Could not load csv file, error: {e}")


    def gather_proteins(self, uniprot_df, skip_if_found=True, gap=10):
        '''
        iterate over lines of a DataFrame containing PDB structure information.
        download every stucture, and curate it if necessary
        '''
        if self.PDB_only:
            try:
                # when the inputs are PDB codes only, convert them to a dataframe
                dic = {'PDB_Code':uniprot_df}
                uniprot_df = pd.DataFrame(dic)
            except Exception as e:
                return e

        for index, row in uniprot_df.iterrows():

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
                    # possible bug fixed
                    af.download_AF_struc(pdb_code,outfolder=self.outdir)
                except Exception as e:
                    print(f">> FAILED: {e}")
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
                    except:
                        pass

            else:
                try:
                    method_obtained = row["Method"]
                    resolution = row["Resolution"]
                    chains = row["Chains"]
                except:
                    pass
                if skip_if_found:
                    files=[os.path.basename(c).split("-")[0] for c in glob.glob(os.path.join(self.curated_dir, "*pdb"))]
                    if pdb_code in files:
                        print(f">> curated {pdb_code} PDB found, continuing...")
                        if not self.PDB_only:
                            data = {'Uniprot_Entry': uniprot_code, 'PDB_Code': pdb_code, 'Method': method_obtained, 'Resolution': resolution, 'Chains': chains}
                        else:
                            data = {'PDB_Code': pdb_code}
                        self.df = pd.concat([self.df, pd.DataFrame.from_records(data, index=[0])], ignore_index=True)
                        continue

                # load, clean, and split it in alternate conformations
                try:
                    self.clean_and_split_pdb(pdb_code)
                    if not self.PDB_only:
                        data = {'Uniprot_Entry': uniprot_code, 'PDB_Code': pdb_code, 'Method': method_obtained, 'Resolution': resolution, 'Chains': chains}
                    else:
                        data = {'PDB_Code': pdb_code}
                    self.df = pd.concat([self.df, pd.DataFrame.from_records(data, index=[0])], ignore_index=True)

                except Exception as e:
                    print(f">> FAILED: {e}")
                    continue


    def clean_and_split_pdb(self, pdb):
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
            raise Exception(f'Error cleaning {pdb}: {e}') from e



        files = glob.glob(os.path.join(self.raw_dir, f"*{pdb}*pdb"))
        test = False
        for cnt, f in enumerate(files):
            mypath = os.path.split(f)[0]
            fasta = os.path.join(mypath, f"{pdb}.fasta")

            try:

                fname = patcher.curate(f, fasta, outdir=self.curated_dir, gap=self.gap)
                if len(replacement_dict) > 0:
                    reverse_replacement_dict = dict((v,k) for k,v in replacement_dict.items())
                    self.replace_chains(fname, reverse_replacement_dict)
                    # TODO 30.07.25 - test this with replacing to new_replace_chains

                test = True

            except Exception as e:
                print(f">> Patching failed for conformer {cnt}. {e}")
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
        cwd = os.getcwd()

        # check if the file has already been downloaded
        files=glob.glob(os.path.join(self.raw_dir, "*pdb"))
        test_file = os.path.join(self.raw_dir, f'{pdb}.pdb')
        if test_file not in files:
            #go into [[outfolder]/conformations and downloads the .pdb file.
            print(f"> downloading PDB {pdb}")
            os.chdir(self.raw_dir)
            tool = get_download_tool()
            try:
                if tool == "curl":
                    line = f"curl -s -o {pdb}.pdb https://files.rcsb.org/download/{pdb}.pdb"
                elif tool == "wget":
                    line = f"wget https://files.rcsb.org/download/{pdb}.pdb"

                else:
                    raise RuntimeError("You don't have a commandline tool for downloading files")

                subprocess.check_call(line, shell=True)
                os.chdir(cwd)

            except Exception as e:
                print(f'Error downloading file. {e}')
                os.chdir(cwd)
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
        cwd = os.getcwd()

        tool = get_download_tool()
        # check if the file has already been downloaded
        files=[c.split(os.sep)[-1][:4] for c in glob.glob(os.path.join(self.raw_dir, "*.fasta"))]
        if pdb not in files:
            try:

                print(f"> downloading FASTA for {pdb}")
                os.chdir(self.raw_dir)

                web_url = f'https://www.rcsb.org/fasta/entry/{pdb}/download'
                name = pdb + '.fasta'

                if tool == "curl":
                    line = f"curl -s -o {name} {web_url}"
                elif tool == "wget":
                    line = f"wget -O {name} {web_url}"
                else:
                    raise RuntimeError("You don't have a commandline tool for downloading files")

                subprocess.check_call(line, shell=True)
                os.chdir(cwd)

            except Exception as e:
                os.chdir(cwd)
                raise Exception(f'Failed downloading FASTA sequence for chain name comparison.{e}') from e
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
            self.replace_selenocysteine(path, pdb)
            if len(replacement_dict) > 0:
                self.new_replace_chains(pdb, path, replacement_dict)
                #self.replace_chains(path, replacement_dict)

        except Exception as e:
            raise Exception(f'Error renaming chains. {e}')

        #Next it opens and starts reading the .pdb file and starts writing a new file with the ending '-clean.pdb'.
        list_of_metals = ['ZN', 'NI', 'CU', 'FE', 'MG', 'MN', 'NA', 'K', 'CA', 'CO', 'CL', 'MO']

        read_file_path = os.path.join(self.raw_dir, f"{pdb}.pdb")
        read_file = open(read_file_path)

        write_file_path = os.path.join(self.raw_dir, f"{pdb}-clean.pdb")
        write_file = open(write_file_path, 'w')

        # Write the clean file, including HETATMs (if they are metal ions),
        # all atoms and lines starting with TER and END.
        test_MSE = False
        test_KCX = False
        for line in read_file:

            # replace selenomethionine with methionine
            if "MSE" in line:
                line = line.replace("HETATM", "ATOM  ")
                line = line.replace('MSE', 'MET')
                line = line.replace('SE', ' S')
                test_MSE = True

            #transform carboxylated lysine into a normal lysine
            if "KCX" in line:

                if ("CX" in line) or ("OQ1" in line) or ("OQ2" in line):
                    continue

                else:
                    line = line.replace("HETATM", "ATOM  ")
                    line = line.replace('KCX', 'LYS')  # TODO: change to investigate as changed 'MES' to 'LYS' - this should be 'LYS' but wasnt before for some reason

                test_KCX = True

            words = line.split()

            #neglect HETATM atoms, unless they are metal ions
            if line[:6] == 'HETATM':    
                if words[3] in list_of_metals:
                    write_file.write(line)
                    continue

            #ignore hydrogen atoms
            if line[:4] == 'ATOM':
                if words[2] != 'H' and words[-1] != 'H':
                    write_file.write(line)
                    continue

            #same terminal statements
            if words[0] == 'END' or words[0] == 'TER' or words[0] == 'ENDMDL':
                write_file.write(line)
                continue

        if test_MSE:
            print(">> mutated MSE to MET")

        if test_KCX:
            print(">> mutated removed a lysine carboxylation")

        write_file.close()
        read_file.close()

        #self.rewrite_pdb(write_file_path)
        #remove the original pdb file as it is not needed anymore.
        os.remove(read_file_path)

        return replacement_dict


    def split_struc_nmr(self, pdb):
        '''
        write a new file for each alternate NMR structure.
        '''

        #define the name format which will be followed for each file.
        number = 1
        name = pdb + '-alt-' + str(number) + '.pdb'

        path = os.path.join(self.raw_dir, f"{pdb}-clean.pdb")

        #Next it opens the file produced from the cleaning script and opens a new file to write in.
        f = open(path)
        endmdls = list()

        #Next it searches to see if there are any 'ENDMDL' statements in the folder (i.e. if there are multiple models).
        for line in f:
            if re.search('ENDMDL', line):
                endmdls.append(line)

        f.close()

        #If there are no 'ENDMDL' statements the clean file is renamed to suit the new format.
        if len(endmdls) == 0:

            path_rename = os.path.join(self.raw_dir, name)
            if os.path.exists(path_rename):
                os.remove(path_rename)

            os.rename(path, path_rename)

        #If there are 'ENDMDL' statements a new file is written for each model.
        elif len(endmdls) != 0:

            print("> Alternate model(s) found. Splitting...")

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

                    if re.search('END ', line):
                        f.close()
                        f_write.close()
                        os.remove(path)
                        os.remove(path_rename)
                        break

                    else:

                        if re.search('ENDMDL', line):
                            f_write.write(line)
                            f_write.close()
                            number = number + 1
                            name = pdb + '-alt-' + str(number) + '.pdb'

                            path_rename = os.path.join(self.raw_dir, name)
                            f_write = open(path_rename, 'w')

                        else:
                            f_write.write(line)
                except:
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

            ABC_list = ['A', 'B', 'C', 'D']
            ABC_dict = {"A": 0, "B": 0, "C": 0}

            #Next it checks if there are any alternate amino acid conformations present (i.e. if line[16 == A, B or C]).
            for i, ABC in enumerate(ABC_list):

                read_file = open(f)
                for line in read_file:
                    if (line[:4] == 'ATOM') and (line[15] == ABC):  # TODO NEED TO ENSURE THIS WORKS FOR ALL FILE WITH NEW ROUTINES
                        ABC_dict[ABC] = 1

            list_of_values = ABC_dict.values()

            #If there are none for this file it moves on to the next.
            if 1 not in list_of_values:
                continue
            else:
                print("> alternate amino acid conformations found. Splitting...")

            #If there are it rewrites a new file for all the As, Bs and Cs.
            for i, ABC in enumerate(ABC_list):
                try:
                    if ABC_dict.get(ABC) == 1:
                        write_path = f[:-6] + f[-5] + ABC + '.pdb'
                        f_write = open(write_path, 'w')

                        target_letter = ABC
                        non_target_letters = []
                        for letter in ABC_list:
                            if letter != target_letter:
                                non_target_letters.append(letter)

                        read = open(f)
                        for line in read:

                            if (line[:4] == 'ATOM') and (line[15] == target_letter):
                                newline = line[:15] + ' ' + line[16:]
                                f_write.write(newline)
                                continue

                            if (line[:4] == 'ATOM') and (line[15] in non_target_letters):
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

                    raise Exception(f"{e}")

            #The original is then removed if it has been replaced.
            read_file.close()
            os.remove(f)

        return


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
            raise Exception(f"Could not download FASTA file. {e}")

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
            print(f'FASTA file parsing failed. {e}')
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
                print(f'Failed renaming chains. {e}')
                continue

        #dictionary used to replace the auth chain names with the RCSB chain names
        return replacement_dict


    def replace_chains(self, path, replacement_dict):
        '''
        
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
                                    #line = line[:21] + replacement_chain_name + ' ' + line[22:]
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
            raise Exception(f'Failed replacing chains. {e}')

        return


    def new_replace_chains(self, pdb, path, replacement_dict):
        '''
        New version of the replace chains which correctly produces pdb files afterwards
        for the patching as the original left some proteins in space. 

        Parameters
        ----------
        pdb : string
            The name of the pdb file that the replacement is being done on

        path : string
            The path of the pdb file that the work is being done on

        replacement_dict : dict
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
            raise Exception(f'Failed replacing chains. {e}')

        return

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
        raise Exception(f'Failed rewriting pdb file. {e}')



    def replace_selenocysteine(self, path, pdb):
        '''
        Identify all the selenocysteine residues within the structure with cyteine in order
        for the model to process this protein structure.

        Parameters
        ----------
        path -> string
            the path to the PDB file to sort

        pdb -> string
            the name of the PDB file
        
        Example
        -------
        >> self.replace_selenocysteine(path, pdb)
        '''
        # replace the U residues in the pdb file
        try:
            M = bb.Molecule(path)
            indices = M.atomselect('*', '*', '*', True, False)[1]

            for index, row in M.data.iterrows():
                if row['resname'] == 'U':
                    replacement_chain_name = 'C'
                    M.data.at[index, 'resname'] = replacement_chain_name


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

        except Exception as e:
            raise Exception(f'Failed to replace selenocysteine in pdb file: {e}') from e


        # replace the U residues in the fasta file
        try:

            fasta = os.path.join(self.raw_dir, f"{pdb}.fasta")

            fin = open(fasta, "r")
            headers = [] # fasta headers
            sequences = [] #collection of sequences
            for line in fin:
                if ">" in line:
                    headers.append(line)
                else:
                    temp_sequence = line
                    for i, res in enumerate(temp_sequence):
                        if res == 'U':
                            temp_sequence = temp_sequence[:i] + 'C' + temp_sequence[(i+1):]
                    sequences.append(temp_sequence)

            new_fasta = open(fasta, 'w')
            for header, sequence in zip(headers, sequences):
                new_fasta.write(header)
                new_fasta.write(sequence)
            new_fasta.close()

        except Exception as e:
            raise Exception(f'Failed to replace selenocysteine in fasta file: {e}') from e




if __name__ == "__main__":

    PDB = PDB()
    if True:
        PDB.clean_and_split_pdb('6XZ7') # test MSE to MET mutation
        #PDB.clean_and_split_pdb('2MBH') # test splitting of models
        #PDB.clean_and_split_pdb('1U8F') # test splitting rotamers

    if False:

        from uniprot import Uniprot
        UP = Uniprot()
        UP.get_protein_data("P09167") # load strucutres for a single UNIPROT
        UP.from_csv_file("inputs\\input_codes_4.csv") # add structures from a .csv file
        print(UP.df)

        PDB.gather_proteins(UP.df)
        print(PDB.df)
