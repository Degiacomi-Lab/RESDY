import unittest
import sys
import os
import shutil
import tempfile
import pandas as pd
import biobox as bb

sys.path.insert(0, os.path.join(os.path.dirname(sys.path[0]), "src"))
import resdy as RD

class Test_Protein(unittest.TestCase):
    def setUp(self):
        self.tmpdir = tempfile.mkdtemp(prefix='resdy_test_')
        self.outdir = self.tmpdir
        shutil.copyfile(os.path.join('demo', 'demo_input.csv'), os.path.join(self.tmpdir, 'demo_input.csv'))
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
                           minimise_strucs=None,
                           max_nmr_conformers=3)

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

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
        self._assert_sulfur_named(f'{self.outdir}{os.sep}curated{os.sep}1PAE-alt-1.pdb', 'CYS', 'SG')

    def test_clean_split_MSEmut(self):
        # test MSE to MET mutation
        self.PDB.clean_and_split_pdb('1A8O', 'P12497')
        self.assertTrue(os.path.isfile(f'{self.outdir}{os.sep}curated{os.sep}1A8O-alt-1.pdb'))
        M_1a8o = bb.Molecule(f'{self.outdir}{os.sep}curated{os.sep}1A8O-alt-1.pdb')
        self.assertFalse(any(a in ['UNK', 'MSE'] for a in list(M_1a8o.data['resname'].unique())))
        self._assert_sulfur_named(f'{self.outdir}{os.sep}curated{os.sep}1A8O-alt-1.pdb', 'MET', 'SD')

    def _assert_sulfur_named(self, path, resname, sulfur):
        '''
        every residue resname in path has its sulfur named sulfur, none is named S, and the
        sulfur keeps the deposited coordinates (occupancy above 0)
        '''
        M = bb.Molecule(path)
        d = M.data[M.data['resname'] == resname]
        self.assertGreater(len(d), 0)
        self.assertNotIn('S', set(d['name']))
        for _, residue in d.groupby(['chain', 'resid']):
            self.assertIn(sulfur, set(residue['name']))
            self.assertTrue((residue.loc[residue['name'] == sulfur, 'occupancy'] > 0).all())

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

    # ---- heteroatom retention policy -------------------------------------------

    @staticmethod
    def _hetatm_counts(path):
        import collections
        return collections.Counter(l[17:20].strip() for l in open(path)
                                   if l.startswith('HETATM'))

    def _curate_1a6m(self, **kwargs):
        '''Curate 1A6M (HEM, waters, SO4, OXY; one chain, no gaps) under a policy.'''
        out = tempfile.mkdtemp(prefix='resdy_het_')
        P = RD.PDB(outdir=out, gap=10, parallel=False, minimise_strucs=None, **kwargs)
        P.clean_and_split_pdb('1A6M', 'P02185', chains=[])
        f = sorted(a for a in os.listdir(os.path.join(out, 'curated')) if '1A6M' in a and a.endswith('.pdb'))[0]
        return P, out, os.path.join(out, 'curated', f)

    def test_hetatm_default_keeps_nothing(self):
        P, out, f = self._curate_1a6m()
        try:
            self.assertFalse(P.include_hetatm)
            self.assertEqual(self._hetatm_counts(f), {})
        finally:
            shutil.rmtree(out, ignore_errors=True)

    def test_hetatm_keep_waters(self):
        P, out, f = self._curate_1a6m(keep_waters=True)
        try:
            self.assertTrue(P.include_hetatm)
            counts = self._hetatm_counts(f)
            self.assertGreater(counts.get('HOH', 0), 100)
            self.assertNotIn('HEM', counts)
        finally:
            shutil.rmtree(out, ignore_errors=True)

    def test_hetatm_keep_named_ligand(self):
        P, out, f = self._curate_1a6m(keep_ligands=('HEM',))
        try:
            counts = self._hetatm_counts(f)
            self.assertEqual(counts.get('HEM', 0), 43)
            self.assertNotIn('HOH', counts)
            self.assertNotIn('SO4', counts)
        finally:
            shutil.rmtree(out, ignore_errors=True)

    def test_hetatm_keep_all(self):
        P, out, f = self._curate_1a6m(keep_ligands='all', keep_waters=True)
        try:
            counts = self._hetatm_counts(f)
            for resname in ('HEM', 'HOH', 'SO4', 'OXY'):
                self.assertIn(resname, counts, f'{resname} should have been kept')
        finally:
            shutil.rmtree(out, ignore_errors=True)

    def test_include_hetatm_is_backwards_compatible(self):
        '''include_hetatm=True must mean what it always documented: keep metal ions.'''
        from resdy.protein import DEFAULT_KEPT_IONS
        P = RD.PDB(outdir=self.outdir, gap=10, include_hetatm=True, minimise_strucs=None)
        self.assertEqual(P.keep_ions, DEFAULT_KEPT_IONS)
        self.assertTrue(P.include_hetatm)
        self.assertFalse(P.keep_waters)

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

    def test_gathering(self):
        UP = RD.Uniprot()
        UP.from_csv_file(f'{self.tmpdir}{os.sep}demo_input.csv')

        gap = 10
        parallel = False
        PDB_only = False
        include_hetatm = False
        resnames_of_interest = ['LYS']
        # minimisation is tested separately, in test_minimisation.py
        minimise_strucs=None
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


class Test_SeleniumToSulfur(unittest.TestCase):
    '''selenomethionine and selenocysteine take the sulfur names of their parents'''

    MSE_SE = 'HETATM  465 SE   MSE A  62      19.231  25.744  31.405  1.00 20.00          SE  \n'
    MSE_CA = 'HETATM  462  CA  MSE A  62      21.331  22.844  28.505  1.00 20.00           C  \n'
    SEC_SE = 'HETATM 1027 SE   SEC A 135      10.100  11.200  12.300  1.00 15.00          SE  \n'

    def test_selenium_to_sulfur(self):
        from resdy.protein import selenium_to_sulfur
        met = selenium_to_sulfur(self.MSE_SE, 'MET', 'SD')
        self.assertEqual(met, self.MSE_SE.replace('SE   MSE', ' SD  MET').replace('SE  \n', ' S  \n'))
        cys = selenium_to_sulfur(self.SEC_SE, 'CYS', 'SG')
        self.assertEqual((cys[17:20], cys[12:16], cys[76:78]), ('CYS', ' SG ', ' S'))
        # every other atom only changes residue name
        self.assertEqual(selenium_to_sulfur(self.MSE_CA, 'MET', 'SD'),
                         self.MSE_CA.replace('MSE', 'MET'))


if __name__ == "__main__":
    unittest.main()
