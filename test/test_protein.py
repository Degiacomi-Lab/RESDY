import unittest
import sys
import os
import re
import shutil
import tempfile
import pandas as pd
import biobox as bb

sys.path.insert(0, os.path.join(os.path.dirname(sys.path[0]), "src"))
import resdy as RD

class Test_Protein(unittest.TestCase):
    def setUp(self):
        self.outdir = tempfile.mkdtemp(prefix='resdy_test_protein_')
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

    def _check_single_conformer(self, path):
        '''
        A conformation file written by the rotamer split holds one conformer: no alternate
        location letter is left, no atom is written twice, and every chain is one contiguous
        block followed by one TER record.
        '''
        lines = open(path).read().splitlines()
        atoms = [l for l in lines if l.startswith(('ATOM', 'HETATM'))]
        self.assertTrue(len(atoms) > 0)
        self.assertEqual({l[16] for l in atoms}, {' '})
        keys = [(l[21], l[22:27], l[12:16]) for l in atoms]
        self.assertEqual(len(keys), len(set(keys)))
        chains = [l[21] for l in atoms]
        n_chains = len(set(chains))
        self.assertEqual(sum(a != b for a, b in zip(chains, chains[1:])), n_chains - 1)
        self.assertEqual(sum(l.startswith('TER') for l in lines), n_chains)

    def test_downloads(self):
        self.PDB.download_pdb('1PAE')
        self.assertTrue(os.path.isfile(f'{self.outdir}{os.sep}conformations{os.sep}1PAE.pdb'))

        self.PDB.download_fasta('1PAE')
        self.assertTrue(os.path.isfile(f'{self.outdir}{os.sep}conformations{os.sep}1PAE.fasta'))

    def test_clean_split_rmv_KCXmut(self):
        # test KCX to LYS mutation
        self.PDB.clean_and_split_pdb('3Q7V', 'Q7WU28')
        self.assertTrue(os.path.isfile(f'{self.outdir}{os.sep}curated{os.sep}3Q7V-alt1A.pdb'))
        M_13ld = bb.Molecule(f'{self.outdir}{os.sep}curated{os.sep}3Q7V-alt1A.pdb')
        self.assertFalse(any(a in ['UNK', 'KCX'] for a in list(M_13ld.data['resname'].unique())))

    '''
    def test_clean_split_keep_KCXmut(self):
        # test keeping the KCX mutation
        PDB = RD.PDB(outdir=self.outdir,
                    gap=self.gap,
                    parallel=self.parallel,
                    PDB_only=self.PDB_only,
                    include_hetatm=self.include_hetatm,
                    resnames_of_interest=self.resnames_of_interest,
                    minimise_strucs=None,
                    max_nmr_conformers=3,
                    keep_modifications='KCX')
        PDB.clean_and_split_pdb('3Q7V', 'Q7WU28')
        self.assertTrue(os.path.isfile(f'{self.outdir}{os.sep}conformations{os.sep}3Q7V-alt1A.pdb'))
        M_13ld = bb.Molecule(f'{self.outdir}{os.sep}conformations{os.sep}3Q7V-alt1A.pdb')
        self.assertTrue(any(a in ['UNK', 'KCX'] for a in list(M_13ld.data['resname'].unique())))
    '''

    def test_clean_split_rmv_SECmut(self):
        # test SEC to CYS mutation
        self.PDB.clean_and_split_pdb('1PAE', 'P22887')
        self.assertTrue(os.path.isfile(f'{self.outdir}{os.sep}curated{os.sep}1PAE-alt-1.pdb'))
        M_1pae = bb.Molecule(f'{self.outdir}{os.sep}curated{os.sep}1PAE-alt-1.pdb')
        self.assertFalse(any(a in ['UNK', 'SEC'] for a in list(M_1pae.data['resname'].unique())))
        self._assert_sulfur_named(f'{self.outdir}{os.sep}curated{os.sep}1PAE-alt-1.pdb', 'CYS', 'SG')

    '''
    def test_clean_split_keep_SECmut(self):
        # test keeping the SEC mutation
        PDB = RD.PDB(outdir=self.outdir,
                    gap=self.gap,
                    parallel=self.parallel,
                    PDB_only=self.PDB_only,
                    include_hetatm=self.include_hetatm,
                    resnames_of_interest=self.resnames_of_interest,
                    minimise_strucs=None,
                    max_nmr_conformers=3,
                    keep_modifications='SEC')
        PDB.clean_and_split_pdb('1PAE', 'P22887')
        self.assertTrue(os.path.isfile(f'{self.outdir}{os.sep}conformations{os.sep}1PAE-alt-1.pdb'))
        M_1pae = bb.Molecule(f'{self.outdir}{os.sep}conformations{os.sep}1PAE-alt-1.pdb')
        self.assertTrue(any(a in ['UNK', 'SEC'] for a in list(M_1pae.data['resname'].unique())))
        self._assert_sulfur_named(f'{self.outdir}{os.sep}conformations{os.sep}1PAE-alt-1.pdb', 'CYS', 'SG')
    '''

    def test_clean_split_rmv_MSEmut(self):
        # test MSE to MET mutation
        self.PDB.clean_and_split_pdb('1A8O', 'P12497')
        self.assertTrue(os.path.isfile(f'{self.outdir}{os.sep}curated{os.sep}1A8O-alt-1.pdb'))
        M_1a8o = bb.Molecule(f'{self.outdir}{os.sep}curated{os.sep}1A8O-alt-1.pdb')
        self.assertFalse(any(a in ['UNK', 'MSE'] for a in list(M_1a8o.data['resname'].unique())))
        self._assert_sulfur_named(f'{self.outdir}{os.sep}curated{os.sep}1A8O-alt-1.pdb', 'MET', 'SD')

    def test_clean_split_keep_MSEmut(self):
        # test keeping the MSE mutation in the structure
        PDB = RD.PDB(outdir=self.outdir,
            gap=self.gap,
            parallel=self.parallel,
            PDB_only=self.PDB_only,
            include_hetatm=self.include_hetatm,
            resnames_of_interest=self.resnames_of_interest,
            minimise_strucs=None,
            max_nmr_conformers=3,
            keep_modifications='MSE')
        PDB.clean_and_split_pdb('1A8O', 'P12497')
        self.assertTrue(os.path.isfile(f'{self.outdir}{os.sep}curated{os.sep}1A8O-alt-1.pdb'))
        M_1a8o = bb.Molecule(f'{self.outdir}{os.sep}curated{os.sep}1A8O-alt-1.pdb')
        self.assertTrue(any(a in ['UNK', 'MSE'] for a in list(M_1a8o.data['resname'].unique())))

    def test_clean_split_keepallmod(self):
        # test the 'all' flag in keep_modifications
        PDB = RD.PDB(outdir=self.outdir,
            gap=self.gap,
            parallel=self.parallel,
            PDB_only=self.PDB_only,
            include_hetatm=self.include_hetatm,
            resnames_of_interest=self.resnames_of_interest,
            minimise_strucs=None,
            max_nmr_conformers=3,
            keep_modifications='all')
        PDB.clean_and_split_pdb('1A8O', 'P12497')
        self.assertTrue(os.path.isfile(f'{self.outdir}{os.sep}curated{os.sep}1A8O-alt-1.pdb'))
        M_1a8o = bb.Molecule(f'{self.outdir}{os.sep}curated{os.sep}1A8O-alt-1.pdb')
        self.assertTrue(any(a in ['UNK', 'MSE'] for a in list(M_1a8o.data['resname'].unique())))

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

    '''
    def test_nmr_conformer_selection_all(self):
        PDB = RD.PDB(outdir=self.outdir,
            gap=self.gap,
            parallel=self.parallel,
            PDB_only=self.PDB_only,
            include_hetatm=self.include_hetatm,
            resnames_of_interest=self.resnames_of_interest,
            minimise_strucs=None,
            max_nmr_conformers=5,
            nmr_conformer_selection='all')
        PDB.clean_and_split_pdb('2MWS', 'P0CG48')
        self.assertTrue(all(os.path.isfile(f'{self.outdir}{os.sep}curated{os.sep}2MWS-alt-{str(a)}.pdb')
                            for a in range(1,9)))
    '''

    def test_nmr_conformer_selection_first(self):
        PDB = RD.PDB(outdir=self.outdir,
            gap=self.gap,
            parallel=self.parallel,
            PDB_only=self.PDB_only,
            include_hetatm=self.include_hetatm,
            resnames_of_interest=self.resnames_of_interest,
            minimise_strucs=None,
            max_nmr_conformers=2,
            nmr_conformer_selection='first')
        PDB.clean_and_split_pdb('2MWS', 'P0CG48')
        self.assertTrue(all(os.path.isfile(f'{self.outdir}{os.sep}curated{os.sep}2MWS-alt-{str(a)}.pdb') for a in range(1,3)))
        self.assertEqual(len([a for a in os.listdir(f'{self.outdir}{os.sep}curated') if re.match('2MWS-alt-*..pdb', a)]), 2)

    def test_nmr_conformer_selection_representative(self):
        PDB = RD.PDB(outdir=self.outdir,
            gap=self.gap,
            parallel=self.parallel,
            PDB_only=self.PDB_only,
            include_hetatm=self.include_hetatm,
            resnames_of_interest=self.resnames_of_interest,
            minimise_strucs=None,
            max_nmr_conformers=10,
            nmr_conformer_selection='representative')
        PDB.clean_and_split_pdb('2MWS', 'P0CG48')
        self.assertTrue(os.path.isfile(f'{self.outdir}{os.sep}curated{os.sep}2MWS-alt-1.pdb'))
        self.assertEqual(len([a for a in os.listdir(f'{self.outdir}{os.sep}curated')
                              if re.match('2MWS-alt-*..pdb', a)]), 1)

    def test_nmr_conformer_selection_list(self):
        PDB = RD.PDB(outdir=self.outdir,
            gap=self.gap,
            parallel=self.parallel,
            PDB_only=self.PDB_only,
            include_hetatm=self.include_hetatm,
            resnames_of_interest=self.resnames_of_interest,
            minimise_strucs=None,
            max_nmr_conformers=10,
            nmr_conformer_selection=[1,3])
        PDB.clean_and_split_pdb('2MWS', 'P0CG48')
        self.assertTrue(all(os.path.isfile(f'{self.outdir}{os.sep}curated{os.sep}2MWS-alt-{str(a)}.pdb') for a in [1,3]))

    def test_rotamer_selection_all(self):
        # test keeping all available rotamers from the pdb structure
        PDB = RD.PDB(outdir=self.outdir,
                    gap=self.gap,
                    parallel=self.parallel,
                    PDB_only=self.PDB_only,
                    include_hetatm=self.include_hetatm,
                    resnames_of_interest=self.resnames_of_interest,
                    minimise_strucs=None,
                    max_nmr_conformers=3,
                    rotamer_selection='all')
        PDB.clean_and_split_pdb('3Q7V', 'Q7WU28')
        self.assertTrue(os.path.isfile(f'{self.outdir}{os.sep}conformations{os.sep}3Q7V-alt1A.pdb'))
        self.assertTrue(os.path.isfile(f'{self.outdir}{os.sep}conformations{os.sep}3Q7V-alt1B.pdb'))
        self._check_single_conformer(f'{self.outdir}{os.sep}conformations{os.sep}3Q7V-alt1A.pdb')
        self._check_single_conformer(f'{self.outdir}{os.sep}conformations{os.sep}3Q7V-alt1B.pdb')

    def test_rotamer_selection_first(self):
        # test keeping only the first rotamer from the pdb structure
        PDB = RD.PDB(outdir=self.outdir,
                    gap=self.gap,
                    parallel=self.parallel,
                    PDB_only=self.PDB_only,
                    include_hetatm=self.include_hetatm,
                    resnames_of_interest=self.resnames_of_interest,
                    minimise_strucs=None,
                    max_nmr_conformers=3,
                    rotamer_selection='first')
        PDB.clean_and_split_pdb('3Q7V', 'Q7WU28')
        self.assertTrue(os.path.isfile(f'{self.outdir}{os.sep}conformations{os.sep}3Q7V-alt1A.pdb'))
        self.assertFalse(os.path.isfile(f'{self.outdir}{os.sep}conformations{os.sep}3Q7V-alt1B.pdb'))
        self._check_single_conformer(f'{self.outdir}{os.sep}conformations{os.sep}3Q7V-alt1A.pdb')

    def test_rotamer_selection_highest_occupancy(self):
        # test keeping only the first rotamer from the pdb structure
        PDB = RD.PDB(outdir=self.outdir,
                    gap=self.gap,
                    parallel=self.parallel,
                    PDB_only=self.PDB_only,
                    include_hetatm=self.include_hetatm,
                    resnames_of_interest=self.resnames_of_interest,
                    minimise_strucs=None,
                    max_nmr_conformers=3,
                    rotamer_selection='highest_occupancy')
        PDB.clean_and_split_pdb('3Q7V', 'Q7WU28')
        print(os.listdir(f'{self.outdir}{os.sep}curated'))
        self.assertTrue(os.path.isfile(f'{self.outdir}{os.sep}conformations{os.sep}3Q7V-alt1HiOcc.pdb'))
        self._check_single_conformer(f'{self.outdir}{os.sep}conformations{os.sep}3Q7V-alt1HiOcc.pdb')
        # the structure that was curated from it is free of clashes
        df_geometry = pd.read_csv(f'{self.outdir}{os.sep}curation_geometry.csv')
        self.assertEqual(int(df_geometry['n_clashes'].max()), 0)

    def test_rotamer_selection_is_per_residue(self):
        # a file with two chains: A1 has alternates A (0.30) and B (0.70), A2 has none, B1 has
        # alternates B (0.60) and C (0.40) only, B2 has A (0.50) and B (0.50)
        def atom(serial, name, alt, chain, resid, x, occ):
            return (f'ATOM  {serial:5d}  {name:<3s}{alt}ALA {chain}{resid:4d}    '
                    f'{x:8.3f}{0.0:8.3f}{0.0:8.3f}{occ:6.2f}{10.0:6.2f}           C\n')
        text = ''.join([
            atom(1, 'N', ' ', 'A', 1, 0.0, 1.0),
            atom(2, 'CA', 'A', 'A', 1, 1.0, 0.30), atom(3, 'CA', 'B', 'A', 1, 2.0, 0.70),
            atom(4, 'CA', ' ', 'A', 2, 3.0, 1.0),
            'TER\n',
            atom(5, 'CA', 'B', 'B', 1, 4.0, 0.60), atom(6, 'CA', 'C', 'B', 1, 5.0, 0.40),
            atom(7, 'CA', 'A', 'B', 2, 6.0, 0.50), atom(8, 'CA', 'B', 'B', 2, 7.0, 0.50),
            'TER\n'])
        # x coordinate of the CA kept for A1, A2, B1 and B2
        expected = {'first': ('TEST-alt1A.pdb', [1.0, 3.0, 4.0, 6.0]),
                    'highest_occupancy': ('TEST-alt1HiOcc.pdb', [2.0, 3.0, 4.0, 6.0])}
        for selection, (name, xs) in expected.items():
            outdir = os.path.join(self.outdir, selection)
            PDB = RD.PDB(outdir=outdir, minimise_strucs=None, rotamer_selection=selection)
            with open(os.path.join(PDB.raw_dir, 'TEST-alt-1.pdb'), 'w') as fh:
                fh.write(text)
            PDB.split_struc_alt_aa('TEST')
            self.assertEqual(sorted(os.listdir(PDB.raw_dir)), [name])
            path = os.path.join(PDB.raw_dir, name)
            self._check_single_conformer(path)
            ca = [float(l[30:38]) for l in open(path) if l.startswith('ATOM') and l[12:16].strip() == 'CA']
            self.assertEqual(ca, xs)

        # 'all' writes one file per letter of the file, each with the atoms of its own letter
        outdir = os.path.join(self.outdir, 'all')
        PDB = RD.PDB(outdir=outdir, minimise_strucs=None, rotamer_selection='all')
        with open(os.path.join(PDB.raw_dir, 'TEST-alt-1.pdb'), 'w') as fh:
            fh.write(text)
        PDB.split_struc_alt_aa('TEST')
        self.assertEqual(sorted(os.listdir(PDB.raw_dir)),
                         ['TEST-alt1A.pdb', 'TEST-alt1B.pdb', 'TEST-alt1C.pdb'])
        ca = {name: [float(l[30:38]) for l in open(os.path.join(PDB.raw_dir, name))
                     if l.startswith('ATOM') and l[12:16].strip() == 'CA']
              for name in os.listdir(PDB.raw_dir)}
        self.assertEqual(ca, {'TEST-alt1A.pdb': [1.0, 3.0, 6.0],
                              'TEST-alt1B.pdb': [2.0, 3.0, 4.0, 7.0],
                              'TEST-alt1C.pdb': [3.0, 5.0]})

    def test_get_prot_org(self):
        shutil.copyfile(os.path.join('demo', 'curated', 'AF-P40616-F1-model_v6.pdb'),
                        os.path.join(self.outdir, 'curated', 'AF-P40616-F1-model_v6.pdb'))
        f = f'{self.outdir}{os.sep}curated{os.sep}AF-P40616-F1-model_v6.pdb'
        org_dets = self.PDB._get_prot_org(f)
        expected_dets = {'organism_name': 'HOMO SAPIENS', 'organism_taxid': '9606'}
        self.assertTrue(org_dets == expected_dets)

    def test_capped_vs_charged(self):
        for meth in ['capped', 'charged']:
            PDB = RD.PDB(outdir=self.outdir,
                        gap=self.gap,
                        parallel=self.parallel,
                        PDB_only=self.PDB_only,
                        include_hetatm=self.include_hetatm,
                        resnames_of_interest=self.resnames_of_interest,
                        minimise_strucs=None,
                        max_nmr_conformers=1,
                        minimise_termini=meth)
            PDB.clean_and_split_pdb('1PAE', 'P22887')

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
        f = sorted(a for a in os.listdir(os.path.join(out, 'curated'))
                   if '1A6M' in a and a.endswith('.pdb'))[0]
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
