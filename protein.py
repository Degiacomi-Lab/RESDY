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

import pandas as pd
import numpy as np

import alphafold as af # to load alphafold data
import patcher # to patch PDB structures with missing regions
from helper import get_download_tool, ShutUp


class PDB(object):
    
    def __init__(self, outdir="result", gap=10, PDB_only=False):
        
        self._setup(outdir, gap, PDB_only)
        

    def _setup(self, outdir, gap, PDB_only):

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
            columns = ['PDB Code']
            self.df = pd.DataFrame(columns=columns)
        else:
            columns = ['Uniprot_Entry', 'PDB Code', 'Method', 'Resolution', 'Chains']
            self.df = pd.DataFrame(columns=columns)

        #gap to consider as small enough to justify patching
        self.gap = gap


    def save_state(self, outname="proteins.csv"):
        '''
        Save a csv file in output directory
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
        except Exception as e:
            print("Could not load csv file. %s"%e)


    def gather_proteins(self, uniprot_df, skip_if_found=True, gap=10):
        '''
        iterate over lines of a DataFrame containing PDB structure information.
        download every stucture, and curate it if necessary
        '''
        if self.PDB_only:
            try:
                # when the inputs are PDB codes only, convert them to a dataframe
                dic = {'PDB Code':uniprot_df}
                uniprot_df = pd.DataFrame(dic)
            except Exception as e:
                return e
        
        for index, row in uniprot_df.iterrows():
    
            PDBCODE = row["PDB Code"]
            
            if not self.PDB_only:
                uniprot_code = row["Uniprot_Entry"]
                print("\nUNIPROT: %s, PDB: %s"%(uniprot_code, PDBCODE))
            else:
                print("\nPDB: %s"%(PDBCODE))
    
            if PDBCODE[:2] == "AF": 
    
                if skip_if_found:
                    files=[os.path.basename(c).split(".")[0] for c in glob.glob(os.path.join(self.curated_dir, "*pdb"))]
                    if PDBCODE in files:
                        print(">> curated %s PDB found, continuing..."%PDBCODE)
                        if not self.PDB_only:
                            data = ({'Uniprot_Entry': uniprot_code, 'PDB Code': PDBCODE, 'Method': 'Predicted', 'Resolution': np.nan, 'Chains': "A"})
                        else:
                            data = ({'PDB Code': PDBCODE})
                        
                        self.df = pd.concat([self.df, pd.DataFrame.from_records(data, index=[0])], ignore_index=True)
                        continue

                try:
                    """possible bug fixed"""
                    af.download_AF_struc(PDBCODE,outfolder=self.outdir)
                except Exception as e:
                    print(">> FAILED: %s"%e)
                    continue
                
                # check if the AlphaFold file contains ATOM statements
                af_filename = os.path.join(self.curated_dir, "%s.pdb"%PDBCODE)
                fin = open(af_filename, "r")
                test = False
                for line in fin:
                    if line.split()[0] == "ATOM":
                        test = True
                        break
                fin.close()
                
                if test:
                    if not self.PDB_only:
                        data = ({'Uniprot_Entry': uniprot_code, 'PDB Code': PDBCODE, 'Method': 'Predicted', 'Resolution': np.nan, 'Chains': "A"})
                    else:
                        data = ({'PDB Code': PDBCODE})
                    self.df = pd.concat([self.df, pd.DataFrame.from_records(data, index=[0])], ignore_index=True)
                
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
                    if PDBCODE in files:
                        print(">> curated %s PDB found, continuing..."%PDBCODE)
                        if not self.PDB_only:
                            data = ({'Uniprot_Entry': uniprot_code, 'PDB Code': PDBCODE, 'Method': method_obtained, 'Resolution': resolution, 'Chains': chains})
                        else:
                            data = ({'PDB Code': PDBCODE})
                        self.df = pd.concat([self.df, pd.DataFrame.from_records(data, index=[0])], ignore_index=True)
                        continue
    
                # load, clean, and split it in alternate conformations
                try:
                    self.clean_and_split_pdb(PDBCODE)
                    if not self.PDB_only:
                        data = ({'Uniprot_Entry': uniprot_code, 'PDB Code': PDBCODE, 'Method': method_obtained, 'Resolution': resolution, 'Chains': chains})
                    else:
                        data = ({'PDB Code': PDBCODE})
                    self.df = pd.concat([self.df, pd.DataFrame.from_records(data, index=[0])], ignore_index=True)
                
                except Exception as e:
                    print(">> FAILED: %s"%(e))
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
            self.split_struc_NMR(pdb)
            self.split_struc_alt_aa(pdb)
    
        except Exception as e:
            raise Exception('Error cleaning %s: %s'%(pdb, e))
            
        files = glob.glob(os.path.join(self.raw_dir, "*%s*pdb"%pdb))
        test = False
        for cnt, f in enumerate(files):
            mypath = os.path.split(f)[0]
            fasta = os.path.join(mypath, "%s.fasta"%pdb)
            
            try: 
                
                fname = patcher.curate(f, fasta, outdir=self.curated_dir, gap=self.gap)
                
                if len(replacement_dict)>0:
                    reverse_replacement_dict = dict((v,k) for k,v in replacement_dict.items())            
                    self.replace_chains(fname, reverse_replacement_dict)
                    
                test = True
                    
            except Exception as e:
                print(">> Patching failed for conformer %s. %s"%(cnt, e))
                continue

        if not test:
            raise Exception("Patching failed for all conformers")

        return
    

    def download_pdb(self, pdb):
        
        cwd = os.getcwd()

        #go into [[outfolder]/conformations and downloads the .pdb file.
        print("> downloading PDB %s"%pdb)
        os.chdir(self.raw_dir)

        tool = get_download_tool()
        try:
            if tool == "curl":
                line = "curl -s -o %s.pdb https://files.rcsb.org/download/%s.pdb"%(pdb, pdb)
            elif tool == "wget":
                line = "wget https://files.rcsb.org/download/" + pdb + '.pdb'

            else:
                raise RuntimeError("You don't have a commandline tool for downloading files")

            subprocess.check_call(line, shell=True)
            os.chdir(cwd)
            
        except Exception as e:
            print('Error downloading file. %s'%e)
            os.chdir(cwd)
            
    
    def download_fasta(self, pdb):
            
        cwd = os.getcwd()

        tool = get_download_tool()
        
        try:
    
            print("> downloading FASTA for %s"%pdb)
            os.chdir(self.raw_dir)
    
            web_url = "https://www.rcsb.org/fasta/entry/" + pdb + '/download'
            name = pdb + '.fasta'
    
            if tool == "curl":
                line = "curl -s -o " + name + " " + web_url
            elif tool == "wget":
                line = "wget -O " + name + " " + web_url
            else:
                raise RuntimeError("You don't have a commandline tool for downloading files")
                
            subprocess.check_call(line, shell=True)
            os.chdir(cwd)
    
        except Exception as e:
            os.chdir(cwd)
            raise Exception('Failed downloading FASTA sequence for chain name comparison.%s'%e)
    
    
    
    def clean(self, pdb):
        '''
        Rename the protein's chains during the cleaning process so that they match the chain names given in the FASTA file.
        This is required as pdb files name their chains using the 'auth' name and fasta with the normal chain name.
        Therefore to avoid confusion we rename them all to what is used in the fasta file.
        '''
        try:
            replacement_dict = self.get_chain_replacement(pdb)
            if len(replacement_dict) > 0:
                
                path = os.path.join(self.raw_dir, "%s.pdb"%pdb)
                self.replace_chains(path, replacement_dict)
    
        except Exception as e:
            raise Exception('Error renaming chains. %s'%e)
                
        #Next it opens and starts reading the .pdb file and starts writing a new file with the ending '-clean.pdb'.
        list_of_metals = ['ZN', 'NI', 'CU', 'FE', 'MG', 'MN', 'NA', 'K', 'CA', 'CO', 'CL', 'MO']
    
        read_file_path = os.path.join(self.raw_dir, "%s.pdb"%pdb)
        read_file = open(read_file_path)
        
        write_file_path = os.path.join(self.raw_dir, "%s-clean.pdb"%pdb)
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
                    line = line.replace('KCX', 'MET')
    
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
    
        #remove the original pdb file as it is not needed anymore.
        os.remove(read_file_path)
    
        return replacement_dict
    
    
    def split_struc_NMR(self, pdb):
        '''
        write a new file for each alternate NMR structure.
        '''
        
        #define the name format which will be followed for each file.
        number = 1
        name = pdb + '-alt-' + str(number) + '.pdb'
        
        path = os.path.join(self.raw_dir, "%s-clean.pdb"%pdb)
            
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
        list_of_files = glob.glob(os.path.join(self.raw_dir, "*%s*.pdb"%pdb))
        
        for f in list_of_files:
    
                ABC_list = ['A', 'B', 'C', 'D']
                ABC_dict = {"A": 0, "B": 0, "C": 0}
    
                #Next it checks if there are any alternate amino acid conformations present (i.e. if line[16 == A, B or C]).
                for i in range(len(ABC_list)):
                    
                    read_file = open(f)
    
                    for line in read_file:
                        if (line[:4] == 'ATOM') and (line[16] == ABC_list[i]):
                            ABC_dict[ABC_list[i]] = 1
    
                list_of_values = ABC_dict.values()
    
                #If there are none for this file it moves on to the next.
                if 1 not in list_of_values:
                    continue
                else:
                    print("> alternate amino acid conformations found. Splitting...")
    
                #If there are it rewrites a new file for all the As, Bs and Cs.
                for i in range(len(ABC_list)):
                    try:
                        if ABC_dict.get(ABC_list[i]) == 1:
                            write_path = f[:-6] + f[-5] + ABC_list[i] + '.pdb'                
                            f_write = open(write_path, 'w')
    
                            target_letter = ABC_list[i]
                            non_target_letters = []
                            for letter in ABC_list:
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
                        
                        raise Exception("%s"%e)
                        
                #The original is then removed if it has been replaced.
                read_file.close()
                os.remove(f)
                
        return
    
    
    def get_chain_replacement(self, pdb_code):
        
        try:
            self.download_fasta(pdb_code)
        except Exception as e:
            raise Exception("Could not download FASTA file. %s"%e)
    
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
                m = re.findall('Chain[ a-z , 0-9, A-Z \[\]]*|', line)
                for entry in m:
                    if (len(entry) > 0):
                        chains_raw.append(entry)
    
            #The fasta file is then removed
            f.close()
            
        except Exception as e:
            print('FASTA file parsing failed. %s'%e)
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
                print('Failed renaming chains. %s'%e)
                continue
            
        #dictionary used to replace the auth chain names with the RCSB chain names
        return replacement_dict
    
    
    def replace_chains(self, path, replacement_dict):
    
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
                                    line = line[:21] + replacement_chain_name + line[22:]
                                    print(line, end ='')
                                else:
                                    print(line, end = '')
    
                        else:
                            print(line, end='')
                    except:
                        print(line, end='')
    
        #If the protein fails the replacement process it is removed from the conformations folder.
        except Exception as e:
            os.remove(path)
            raise Exception('Failed replacing chains. %s'%e)
    
        return


if __name__ == "__main__":

    PDB = PDB()    
    if True:
        PDB.clean_and_split_pdb('1CI4') # test MSE to MET mutation
        PDB.clean_and_split_pdb('2MBH') # test splitting of models
        PDB.clean_and_split_pdb('1U8F') # test splitting rotamers
      
    if False:

        from uniprot import Uniprot
        UP = Uniprot()
        UP.get_protein_data("P09167") # load strucutres for a single UNIPROT
        UP.from_csv_file("inputs\\input_codes_4.csv") # add structures from a .csv file
        print(UP.df)

        PDB.gather_proteins(UP.df)
        print(PDB.df)
        
        
    
