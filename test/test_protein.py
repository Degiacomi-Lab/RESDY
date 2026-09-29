import unittest
import sys
import os
import shutil
import tempfile
import pandas as pd
import biobox as bb

sys.path.insert(0, os.path.join(os.path.dirname(sys.path[0]), "src"))
import resdy as RD
from resdy import alphafold as af

class Test_Protein(unittest.TestCase):
    def setUp(self):
        self.outdir = 'resdy_test_protein'
        os.makedirs(self.outdir, exist_ok=True)
        shutil.copyfile(os.path.join('demo', 'demo_input.csv'), os.path.join(self.outdir, 'demo_input.csv'))
        self.gap = 10
        self.parallel = False
        self.PDB_only = False
        self.include_hetatm = False
        self.resnames_of_interest = ['LYS']
        self.PDB = RD.PDB(outdir=self.outdir,
                           gap=self.gap,
                           parallel=self.parallel,
                           PDB_only=self.PDB_only,
                           include_hetatm=self.include_hetatm,
                           resnames_of_interest=self.resnames_of_interest,
                           max_nmr_conformers=3)

    def tearDown(self):
        shutil.rmtree(self.outdir, ignore_errors=True)

    def test_downloads(self):
        self.PDB.download_pdb('1PAE')
        self.assertTrue(os.path.isfile(f'{self.outdir}{os.sep}conformations{os.sep}1PAE.pdb'))

        self.PDB.download_fasta('1PAE')
        self.assertTrue(os.path.isfile(f'{self.outdir}{os.sep}conformations{os.sep}1PAE.fasta'))

    def test_clean_split_KCXmut(self):
        # test KCX to LYS mutation
        self.PDB.clean_and_split_pdb('3Q7V', 'Q7WU28')
        self.assertTrue(os.path.isfile(f'{self.outdir}{os.sep}curated{os.sep}3Q7V-alt1A.pdb'))
        M_13ld = bb.Molecule(f'{self.outdir}{os.sep}curated{os.sep}3Q7V-alt1A.pdb')
        self.assertFalse(any(a in ['UNK', 'KCX'] for a in list(M_13ld.data['resname'].unique())))

    def test_clean_split_SECmut(self):
        # test SEC to CYS mutation
        self.PDB.clean_and_split_pdb('1PAE', 'P22887')
        self.assertTrue(os.path.isfile(f'{self.outdir}{os.sep}curated{os.sep}1PAE-alt-1.pdb'))
        M_1pae = bb.Molecule(f'{self.outdir}{os.sep}curated{os.sep}1PAE-alt-1.pdb')
        self.assertFalse(any(a in ['UNK', 'SEC'] for a in list(M_1pae.data['resname'].unique())))

    def test_clean_split_MSEmut(self):
        # test MSE to MET mutation
        self.PDB.clean_and_split_pdb('1A8O', 'P12497')
        self.assertTrue(os.path.isfile(f'{self.outdir}{os.sep}curated{os.sep}1A8O-alt-1.pdb'))
        M_1a8o = bb.Molecule(f'{self.outdir}{os.sep}curated{os.sep}1A8O-alt-1.pdb')
        self.assertFalse(any(a in ['UNK', 'MSE'] for a in list(M_1a8o.data['resname'].unique())))

    def test_clean_split_modelsplit(self):
        # test splitting of models
        self.PDB.clean_and_split_pdb('1A6M', 'Q13351')
        self.assertTrue(os.path.isfile(f'{self.outdir}{os.sep}curated{os.sep}1A6M-alt1A.pdb'))

    def test_clean_split_rotamersplit(self):
        # test splitting rotamers
        self.PDB.clean_and_split_pdb('1A6M', 'P02185')
        self.assertTrue(os.path.isfile(f'{self.outdir}{os.sep}curated{os.sep}1A6M-alt1A.pdb'))

    def test_clean_split_renumbering(self):
        # test renumbering residues with canonical uniprot sequence
        self.PDB.clean_and_split_pdb('3DBJ', 'P50030', chains=['A', 'C', 'E', 'G'])
        self.assertTrue(os.path.isfile(f'{self.outdir}{os.sep}curated{os.sep}3DBJ-alt-1.pdb'))

    def test_clean_split_generalmod(self):
        #  test removal of modified residue
        self.PDB.clean_and_split_pdb('2MWS', 'P0CG48')
        self.assertTrue(os.path.isfile(f'{self.outdir}{os.sep}curated{os.sep}2MWS-alt-1.pdb'))
        M_2mws = bb.Molecule(f'demo{os.sep}curated{os.sep}2MWS-alt-1.pdb')
        self.assertFalse(any(a in ['UNK', '3X9'] for a in list(M_2mws.data['resname'].unique())))

    def test_auxiliary(self):
        # test all other random functions from protein class
        if not os.path.exists(f'{self.outdir}{os.sep}curated{os.sep}13LD-alt1A.pdb'):
            self.PDB.clean_and_split_pdb('13LD', 'P10724')
        self.PDB.rewrite_pdb(path=f'{self.outdir}{os.sep}curated{os.sep}13LD-alt1A.pdb')
        self.assertTrue(os.path.isfile(f'{self.outdir}{os.sep}curated{os.sep}13LD-alt1A.pdb'))

        self.PDB.load_state(fname=f'demo{os.sep}proteins.csv',
                            outdir=self.outdir,
                            PDB_only=self.PDB_only)
        self.PDB.df = pd.DataFrame(columns=['Uniprot_Entry', 'PDB_Code', 'Method', 'Resolution', 'Chains'])

    def test_minimisation(self):
        # test the minimisation of the protein structure with openmm
        code = 'AF-P40616-F1-model_v6'
        af.download_AF_struc(code, outfolder=self.outdir)
        self.PDB.apply_minimisation(code, max_iterations=10)
        self.assertTrue(os.path.exists(f'{self.outdir}{os.sep}curated{os.sep}{code}_relaxed.pdb'))

    def test_gathering(self):
        UP = RD.Uniprot()
        UP.from_csv_file(f'{self.outdir}{os.sep}demo_input.csv')

        gap = 10
        parallel = False
        PDB_only = False
        include_hetatm = False
        resnames_of_interest = ['LYS']
        minimise_strucs='AF'
        PDB = RD.PDB(outdir=self.outdir,
                    gap=gap,
                    parallel=parallel,
                    PDB_only=PDB_only,
                    include_hetatm=include_hetatm,
                    resnames_of_interest=resnames_of_interest,
                    minimise_strucs=minimise_strucs,
                    max_nmr_conformers=2)
        PDB.gather_proteins(UP.df.head(1), skip_if_found=False)
        self.assertEqual(1, len(PDB.df))


if __name__ == "__main__":
    unittest.main()
