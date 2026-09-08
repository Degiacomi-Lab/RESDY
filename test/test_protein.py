import unittest
import sys
import os
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(sys.path[0]), "src"))
import coolpackage as CPN

class Test_Uniprot(unittest.TestCase):
    def setUp(self):
        self.outdir = 'test'
        self.gap = 10
        self.parallel = False
        self.PDB_only = False
        self.include_hetatm = False
        self.resnames_of_interest = ['LYS']
        self.PDB = CPN.PDB(outdir=self.outdir,
                           gap=self.gap,
                           parallel=self.parallel,
                           PDB_only=self.PDB_only,
                           include_hetatm=self.include_hetatm,
                           resnames_of_interest=self.resnames_of_interest)

    def test_downloads(self):
        self.PDB.download_pdb('1PAE')
        self.assertTrue(os.path.isfile(f'{self.outdir}{os.sep}conformations{os.sep}1PAE.pdb'))

        self.PDB.download_fasta('1PAE')
        self.assertTrue(os.path.isfile(f'{self.outdir}{os.sep}conformations{os.sep}1PAE.fasta'))


    def test_clean_split_KCXmut(self):
        # test as many difficult cases here as possible
        self.PDB.clean_and_split_pdb('13LD', 'P10724') # test KCX to LYS mutation
        self.assertTrue(os.path.isfile(f'{self.outdir}{os.sep}curated{os.sep}13LD-alt1A.pdb'))
    
    def test_clean_split_SECmut(self):
        self.PDB.clean_and_split_pdb('1PAE', 'P22887') # test SEC to CYS mutation
        self.assertTrue(os.path.isfile(f'{self.outdir}{os.sep}curated{os.sep}1PAE-alt-1.pdb'))
    
    def test_clean_split_MSEmut(self):
        self.PDB.clean_and_split_pdb('6XZ7', 'P60422') # test MSE to MET mutation
        self.assertTrue(os.path.isfile(f'{self.outdir}{os.sep}curated{os.sep}6XZ7-alt-1.pdb'))
    
    def test_clean_split_modelsplit(self):
        self.PDB.clean_and_split_pdb('1A6M', 'Q13351') # test splitting of models
        self.assertTrue(os.path.isfile(f'{self.outdir}{os.sep}curated{os.sep}1A6M-alt1A.pdb'))
    
    def test_clean_split_rotamersplit(self):
        self.PDB.clean_and_split_pdb('1A6M', 'P02185') # test splitting rotamers
        self.assertTrue(os.path.isfile(f'{self.outdir}{os.sep}curated{os.sep}1A6M-alt1A.pdb'))
    
    def test_clean_split_renumbering(self):
        self.PDB.clean_and_split_pdb('3DBJ', 'P50030', chains=['A', 'C', 'E', 'G']) # test renumbering residues with canonical uniprot sequence
        self.assertTrue(os.path.isfile(f'{self.outdir}{os.sep}curated{os.sep}3DBJ-alt-1.pdb'))
    
    def test_clean_split_generalmod(self):
        self.PDB.clean_and_split_pdb('2MWS', 'P0CG48') #  test removal of modified residue
        self.assertTrue(os.path.isfile(f'{self.outdir}{os.sep}curated{os.sep}2MWS-alt-1.pdb'))

    def test_auxiliary(self):
        self.PDB.rewrite_pdb(path=f'{self.outdir}{os.sep}curated{os.sep}13LD-alt1A.pdb')
        self.assertTrue(os.path.isfile(f'{self.outdir}{os.sep}curated{os.sep}13LD-alt1A.pdb'))

    def test_gathering(self):
        # runs on separate instance of PDB to work with demo
        UP = CPN.Uniprot()
        UP.from_csv_file(f'demo{os.sep}demo_input.csv')
        
        outdir = 'test'
        gap = 10
        parallel = False
        PDB_only = False
        include_hetatm = False
        resnames_of_interest = ['LYS']
        PDB = CPN.PDB(outdir=outdir,
                            gap=gap,
                            parallel=parallel,
                            PDB_only=PDB_only,
                            include_hetatm=include_hetatm,
                            resnames_of_interest=resnames_of_interest)
        PDB.gather_proteins(UP.df, skip_if_found=False)


if __name__ == "__main__":
    unittest.main()