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
                           minimise_strucs=None,
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
        f = sorted(a for a in os.listdir(os.path.join(out, 'curated')) if '1A6M' in a)[0]
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
        UP.from_csv_file(f'{self.outdir}{os.sep}demo_input.csv')

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


class Test_Old_Curation(unittest.TestCase):
    """
    A curated file with an occupancy outside [0, 1] was written before RESDY required
    biobox 1.1.5: it is deleted, with the rest of its structure, rather than reused.
    """

    def setUp(self):
        self.outdir = tempfile.mkdtemp(prefix='resdy_test_')
        self.curated = os.path.join(self.outdir, 'curated')
        os.makedirs(self.curated)

    def tearDown(self):
        shutil.rmtree(self.outdir, ignore_errors=True)

    def _write(self, stem, occupancy):
        path = os.path.join(self.curated, f'{stem}.pdb')
        with open(path, 'w') as fh:
            fh.write(f'ATOM      1  CA  LYS A   1       1.000   2.000   3.000{occupancy:6.2f} 25.00           C\nEND\n')
        return path

    def test_old_layout_is_detected_and_discarded(self):
        from resdy.protein import has_old_column_layout
        new = self._write('1ABC-alt1A', 1.0)
        built = self._write('1ABC-alt1B', 0.0)
        self.assertFalse(has_old_column_layout(new))
        self.assertFalse(has_old_column_layout(built))

        old = self._write('9XYZ-alt1A', 25.0)
        relaxed = self._write('9XYZ-alt1A_relaxed', 1.0)
        record = os.path.join(self.curated, '9XYZ-alt1A.minimisation.json')
        open(record, 'w').write('{}')
        self.assertTrue(has_old_column_layout(old))

        P = RD.PDB(outdir=self.outdir, minimise_strucs=None)
        self.assertFalse(P._discard_old_curation([new, built]))
        self.assertTrue(os.path.exists(new) and os.path.exists(built))
        self.assertTrue(P._discard_old_curation([old, relaxed]))
        for f in (old, relaxed, record):
            self.assertFalse(os.path.exists(f))



if __name__ == "__main__":
    unittest.main()
