import unittest
import sys
import os
import shutil
import tempfile
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(sys.path[0]), "src"))
import resdy as RD

class Test_Measure(unittest.TestCase):
    def setUp(self):
        print('-> Setting up measures tests')
        self.outdir = tempfile.mkdtemp(prefix='resdy_test_')
        shutil.copytree(os.path.join('demo', 'curated'), os.path.join(self.outdir, 'curated'))
        shutil.copyfile(os.path.join('demo', 'demo_input.csv'), os.path.join(self.outdir, 'demo_input.csv'))
        shutil.copyfile(os.path.join('demo', 'proteins.csv'), os.path.join(self.outdir, 'proteins.csv'))
        self.df_input = pd.read_csv(f'{self.outdir}{os.sep}demo_input.csv')
        self.df_prot = pd.read_csv(f'{self.outdir}{os.sep}proteins.csv')
        self.all_features = ['propka', 'sasa', 'depth', 'aev',
                             'das', 'seqcharge', 'melodia', 'frustration']
        self.features_dict = {k: {} for k in self.all_features}
        self.df_pdb_only_test = pd.DataFrame({'PDB_Code': ['1A6M-alt1A.pdb', '1A6M-alt1a_relaxed.pdb']})

    def tearDown(self):
        shutil.rmtree(self.outdir, ignore_errors=True)

    def test_full_setup(self):
        # test full setup of measures class, no run, test runs individually
        # standard track
        print('-> Testing full setup')
        M = RD.Measure(df_input=self.df_prot,
                        outdir=self.outdir,
                        features_dict=self.features_dict)
        reg_feats = {m[0] for m in M.measures}
        dereg_feats = set(self.features_dict) - reg_feats - set(M.features)
        self.assertEqual(dereg_feats, set(), f'> Some features vanished without being removed: {dereg_feats}')

        # pdb only
        M = RD.Measure(df_input=self.df_pdb_only_test,
                        outdir=self.outdir,
                        features_dict=self.features_dict)
        reg_feats = {m[0] for m in M.measures}
        dereg_feats = set(self.features_dict) - reg_feats - set(M.features)
        self.assertEqual(dereg_feats, set(), f'> Some features vanished without being removed: {dereg_feats}')

    def test_propka_standardtrack(self):
        print('-> Test measuring PROPKA - Standard Track')
        M_propka = RD.Measure(df_input=self.df_prot[:2],
                        outdir=self.outdir,
                        features_dict={'propka': {}},
                        parallel=True,
                        include_modified=False)
        M_propka.measure_data()
        print(M_propka.df)
        print('len propka', len(M_propka.df))
        self.assertTrue(len(M_propka.df) > 0)

    def test_propka_pdbonlytrack(self):
        print('-> Test measuring PROPKA - PDB Only')
        M_propka = RD.Measure(df_input=self.df_pdb_only_test,
                        outdir=self.outdir,
                        features_dict={'propka': {}},
                        parallel=True,
                        include_modified=False)
        M_propka.measure_data()
        print(M_propka.df)
        print('len propka', len(M_propka.df))
        self.assertTrue(len(M_propka.df) > 0)

    def test_sasa_standardtrack(self):
        print('-> Test measuring SASA')
        M_sasa = RD.Measure(df_input=self.df_prot[:2],
                        outdir=self.outdir,
                        features_dict={'sasa': {}},
                        parallel=True,
                        include_modified=False)
        M_sasa.measure_data()
        print(M_sasa.df)
        print('len sasa', len(M_sasa.df))
        self.assertTrue(len(M_sasa.df) > 0)

    def test_sasa_pdbonlytrack(self):
        print('-> Test measuring SASA')
        M_sasa = RD.Measure(df_input=self.df_pdb_only_test,
                        outdir=self.outdir,
                        features_dict={'sasa': {}},
                        parallel=True,
                        include_modified=False)
        M_sasa.measure_data()
        print(M_sasa.df)
        print('len sasa', len(M_sasa.df))
        self.assertTrue(len(M_sasa.df) > 0)

    def test_depth_standardtrack(self):
        print('-> Test measuring DEPTH')
        M_depth = RD.Measure(df_input=self.df_prot[:2],
                        outdir=self.outdir,
                        features_dict={'depth': {}},
                        parallel=True,
                        include_modified=False)
        M_depth.measure_data()
        print(M_depth.df)
        print('len depth', len(M_depth.df))
        self.assertTrue(len(M_depth.df) > 0)

    def test_depth_pdbonlytrack(self):
        print('-> Test measuring DEPTH')
        M_depth = RD.Measure(df_input=self.df_pdb_only_test,
                        outdir=self.outdir,
                        features_dict={'depth': {}},
                        parallel=True,
                        include_modified=False)
        M_depth.measure_data()
        print(M_depth.df)
        print('len depth', len(M_depth.df))
        self.assertTrue(len(M_depth.df) > 0)

    def test_aev_standardtrack(self):
        print('-> Test measuring AEV')
        M_aev = RD.Measure(df_input=self.df_prot[:2],
                        outdir=self.outdir,
                        features_dict={'aev': {}},
                        parallel=True,
                        include_modified=False)
        M_aev.measure_data()
        print(M_aev.df)
        print('len aev', len(M_aev.df))
        self.assertTrue(len(M_aev.df) > 0)

    def test_aev_pdbonlytrack(self):
        print('-> Test measuring AEV')
        M_aev = RD.Measure(df_input=self.df_pdb_only_test,
                        outdir=self.outdir,
                        features_dict={'aev': {}},
                        parallel=True,
                        include_modified=False)
        M_aev.measure_data()
        print(M_aev.df)
        print('len aev', len(M_aev.df))
        self.assertTrue(len(M_aev.df) > 0)

    def test_das(self):
        print('-> Test measuring DAS')
        M_das = RD.Measure(df_input=self.df_prot[:2],
                        outdir=self.outdir,
                        features_dict={'das': {}},
                        parallel=True,
                        include_modified=False)
        M_das.measure_data()
        print(M_das.df)
        print('len das', len(M_das.df))
        self.assertTrue(len(M_das.df) > 0)

    def test_charge(self):
        print('-> Test measuring CHARGE')
        M_seqcharge = RD.Measure(df_input=self.df_prot[:2],
                        outdir=self.outdir,
                        features_dict={'seqcharge': {}},
                        parallel=True,
                        include_modified=False)
        M_seqcharge.measure_data()
        print(M_seqcharge.df)
        print('len seqcharge', len(M_seqcharge.df))
        self.assertTrue(len(M_seqcharge.df) > 0)

    def test_flexibility(self):
        print('-> Test measuring FLEXIBILITY')
        M_flex = RD.Measure(df_input=self.df_prot[:2],
                        outdir=self.outdir,
                        features_dict={'flexibility': {}},
                        parallel=True,
                        include_modified=False)
        M_flex.measure_data()
        print(M_flex.df)
        print('len flexibility', len(M_flex.df))

    def test_legolas(self):
        print('-> Test measuring LEGOLAS')
        M_legolas = RD.Measure(df_input=self.df_prot[:2],
                        outdir=self.outdir,
                        features_dict={'legolas': {}},
                        parallel=False,
                        include_modified=False)
        M_legolas.measure_data()
        print(M_legolas.df)
        print('len legolas:', len(M_legolas.df))
    '''
    def test_frustration(self):
        print('-> Test measuring FRUSTRATION')
        M_frustration = RD.Measure(df_input=self.df_prot[:2],
                        outdir=self.outdir,
                        features_dict={'frustration': {}, 'density': {}},
                        parallel=False,
                        include_modified=False)
        M_frustration.measure_data()
        print(M_frustration.df)
        print('len frustration', len(M_frustration.df))
    '''
    def test_melodia(self):
        print('-> Test measuring MELODIA')
        M_melodia = RD.Measure(df_input=self.df_prot[:2],
                        outdir=self.outdir,
                        features_dict={'melodia': {}},
                        parallel=False,
                        include_modified=False)
        M_melodia.measure_data()
        print(M_melodia.df)
        print('len melodia', len(M_melodia.df))
        self.assertTrue(len(M_melodia.df) > 0)

    def test_phi(self):
        print('-> Test measuring PHI')
        M_phi = RD.Measure(df_input=self.df_prot[:2],
                        outdir=self.outdir,
                        features_dict={'phi': {}},
                        parallel=False,
                        include_modified=False)
        M_phi.measure_data()
        print(M_phi.df)
        print('len phi', len(M_phi.df))
        self.assertTrue(len(M_phi.df) > 0)

    def test_evolution(self):
        print('-> Test measuring EVOLUTION')
        M_evolution = RD.Measure(df_input=self.df_prot[:2],
                        outdir=self.outdir,
                        features_dict={'evolution': {}},
                        parallel=True,
                        include_modified=False)
        M_evolution.measure_data()
        print(M_evolution.df)
        print('len evolution', len(M_evolution.df))
        self.assertTrue(len(M_evolution.df) > 0)

    def test_rmsf(self):
        print('-> Test measuring RMSF')
        M_rmsf = RD.Measure(df_input=self.df_prot[:2],
                        outdir=self.outdir,
                        features_dict={'rmsf': {}},
                        parallel=True,
                        include_modified=False)
        M_rmsf.measure_data()
        print(M_rmsf.df)
        print('len rmsf', len(M_rmsf.df))
        self.assertTrue(len(M_rmsf.df) > 0)

    def test_pkaani(self):
        print('-> Test measuring PKAANI')
        M_pkaani = RD.Measure(df_input=self.df_prot[:2],
                        outdir=self.outdir,
                        features_dict={'pkaani': {}},
                        parallel=True,
                        include_modified=False)
        M_pkaani.measure_data()
        print(M_pkaani.df)
        print('len pkaani', len(M_pkaani.df))

    def test_secondarystructure(self):
        print('-> Test measuring SECONDARYSTRUCTURE')
        M_ss = RD.Measure(df_input=self.df_prot[:2],
                        outdir=self.outdir,
                        features_dict={'secondarystructure': {}},
                        parallel=True,
                        include_modified=False)
        M_ss.measure_data()
        print(M_ss.df)
        print('len secondarystructure', len(M_ss.df))

    '''
    def test_restart_measure_data(self):
        print('-> Test restarting measures')
        M = RD.Measure(df_input=self.df_prot,
                        outdir=self.outdir,
                        features_dict={'das': {}},
                        parallel=True,
                        include_modified=False)
        M.restart_measure_data()

    def test_recovery(self):
        print('-> Test recovering data')
        M = RD.Measure(df_input=self.df_prot,
                        outdir=self.outdir,
                        features_dict={'das': {}},
                        parallel=True,
                        include_modified=False)
        M.recover_from_log('measure_log.txt')
    '''

    def test_min_dist_other_chain(self):
        # metadata column: NaN for one chain, a distance otherwise, no row without an anchor
        print('-> Testing Min_Dist_Other_Chain')
        import numpy as np
        import biobox as bb
        rows = self.df_prot[self.df_prot['PDB_Code'].isin(['AF-P40616-F1-model_v6', '1UPT', '8Y8Y'])]
        M = RD.Measure(df_input=rows, outdir=self.outdir, features_dict={'sasa': {}},
                       parallel=False)
        M.measure_data()
        df = M.df
        self.assertIn('Min_Dist_Other_Chain', df.columns)
        self.assertNotIn('Min_Dist_Other_Chain', M.features)

        single = df[df['PDB_Code'].str.startswith(('AF-', '8Y8Y'))]
        self.assertGreater(len(single), 0)
        self.assertTrue(single['Min_Dist_Other_Chain'].isna().all())
        multi = df[df['PDB_Code'].str.startswith('1UPT')]
        self.assertGreater(len(multi), 0)
        self.assertTrue((multi['Min_Dist_Other_Chain'] > 0).all())

        # the five 8Y8Y lysines deposited without NZ cannot be featurised
        y8 = df[df['PDB_Code'].str.startswith('8Y8Y')]
        self.assertFalse(set(y8['Resid'].astype(int)) & {86, 166, 244, 246, 276})

        # brute force on one 1UPT lysine
        r = multi.iloc[0]
        S = bb.Molecule(os.path.join(self.outdir, 'curated', f"{r['PDB_Code']}.pdb"))
        d = S.data
        nz = S.points[((d['chain'] == r['Chain']) & (d['resid'] == int(r['Resid']))
                       & (d['name'] == 'NZ')).to_numpy()][0]
        other = S.points[(d['chain'] != r['Chain']).to_numpy()]
        self.assertAlmostEqual(r['Min_Dist_Other_Chain'],
                               np.linalg.norm(other - nz, axis=1).min(), places=6)

    def test_residue_preset(self):
        # a three-letter code of any standard residue resolves to its preset
        print('-> Testing a residue preset other than lysine')
        M = RD.Measure(df_input=self.df_prot,
                       outdir=self.outdir,
                       features_dict={'sasa': {}},
                       residue_of_interest='tyr')
        self.assertEqual(M.aa_properties['non_modified_codes'], ['TYR'])
        self.assertEqual(M.aa_properties['atom_select_names_nonmod'], ['OH'])

    def test_lysine_preset_unchanged(self):
        # the shipped lysine behaviour has to survive the move into residues.py
        print('-> Testing that the lysine default is unchanged')
        M = RD.Measure(df_input=self.df_prot,
                       outdir=self.outdir,
                       features_dict={'sasa': {}})
        self.assertEqual(M.aa_properties['atom_select_names_nonmod'], ['NZ'])
        self.assertEqual(M.aa_properties['atom_select_names_modified'], ['NZ', 'N07'])
        self.assertIn('LYS', M.aa_properties['non_modified_codes'])
        self.assertIn('KCX', M.aa_properties['modified_codes'])
        self.assertIn('LYE', M.aa_properties['modified_codes'])

    def test_unsupported_features_dropped(self):
        # alanine has no pKa from either backend and no half-sphere shells
        print('-> Testing that features that cannot act on the residue are dropped')
        M = RD.Measure(df_input=self.df_prot,
                       outdir=self.outdir,
                       features_dict={'propka': {}, 'pkaANI': {}, 'das': {}, 'sasa': {}},
                       residue_of_interest='ALA')
        self.assertEqual(set(M.features_dict), {'sasa'})
        self.assertEqual({m[0] for m in M.measures}, {'sasa'})

    def test_supported_features_kept(self):
        print('-> Testing that the pKa backends survive for a residue they support')
        M = RD.Measure(df_input=self.df_prot,
                       outdir=self.outdir,
                       features_dict={'propka': {}, 'pkaANI': {}, 'das': {}},
                       residue_of_interest='TYR')
        self.assertEqual(set(M.features_dict), {'propka', 'pkaANI', 'das'})

    def test_das_radii_injected(self):
        # the residue's own shells reach the feature, and a user setting wins over them
        print('-> Testing that the DAS shells of the residue are passed through')
        features_dict = {'das': {}}
        M = RD.Measure(df_input=self.df_prot,
                       outdir=self.outdir,
                       features_dict=features_dict,
                       residue_of_interest='TRP')
        self.assertEqual(M.features_dict['das']['radii'], [4.7, 4.4, 4.0, 3.6])
        # the caller's own dictionary must not have been written into
        self.assertEqual(features_dict['das'], {})

        M_user = RD.Measure(df_input=self.df_prot,
                            outdir=self.outdir,
                            features_dict={'das': {'radii': [5.0, 4.0, 3.0]}},
                            residue_of_interest='TRP')
        self.assertEqual(M_user.features_dict['das']['radii'], [5.0, 4.0, 3.0])

    def test_custom_dict_any_key_order(self):
        print('-> Testing a custom residue dictionary given in a different key order')
        M = RD.Measure(df_input=self.df_prot,
                       outdir=self.outdir,
                       features_dict={'sasa': {}},
                       residue_of_interest={'atom_select_names_nonmod': ['SG'],
                                            'modified_codes': '',
                                            'non_modified_codes': ['CYS'],
                                            'atom_select_names_modified': ['SG']})
        self.assertEqual(M.aa_properties['modified_codes'], [])
        self.assertEqual(M.aa_properties['atom_select_names_nonmod'], ['SG'])

    def test_plugin_feature_declarations_are_honoured(self):
        '''
        A feature the core has never heard of should be able to declare which residues it can
        act on, and its per-residue settings, without an edit to measure.py or residues.py.
        '''
        print('-> Testing that a new feature can declare its own residue support')
        import resdy.measure as measure_module

        class PLUGINPROBE():
            SUPPORTED_RESIDUES = {'TYR', 'TRP'}
            RESIDUE_KWARGS = {'TYR': {'window': 11}, 'TRP': {'window': 7}}
            UNSUPPORTED_REASON = {'LYS': 'it was written for aromatics'}

            def __init__(self, include_modified=False, aa_properties=None,
                         error_filename='measure_errors.txt', window=0):
                self.window = window

            def calculate(self, path):
                return pd.DataFrame(columns=['Chain', 'Resid', 'pluginprobe'])

        measure_module.PLUGINPROBE = PLUGINPROBE
        try:
            # the residue it supports: kept, and given the settings of that residue
            M = RD.Measure(df_input=self.df_prot, outdir=self.outdir,
                           features_dict={'pluginprobe': {}}, residue_of_interest='TYR')
            self.assertEqual(M.features_dict['pluginprobe'], {'window': 11})
            self.assertIn('pluginprobe', {m[0] for m in M.measures})

            # a residue it excludes: dropped, quoting its own reason
            M_out = RD.Measure(df_input=self.df_prot, outdir=self.outdir,
                               features_dict={'pluginprobe': {}, 'sasa': {}},
                               residue_of_interest='LYS')
            self.assertEqual(set(M_out.features_dict), {'sasa'})

            # what the user asks for still wins over the settings of the residue
            M_user = RD.Measure(df_input=self.df_prot, outdir=self.outdir,
                                features_dict={'pluginprobe': {'window': 99}},
                                residue_of_interest='TRP')
            self.assertEqual(M_user.features_dict['pluginprobe'], {'window': 99})
        finally:
            del measure_module.PLUGINPROBE

    def test_unknown_residue_raises(self):
        print('-> Testing that an unknown residue of interest raises')
        with self.assertRaises(KeyError):
            RD.Measure(df_input=self.df_prot,
                       outdir=self.outdir,
                       features_dict={'sasa': {}},
                       residue_of_interest='ZZZ')


class Test_Flexibility_Built_Atoms(unittest.TestCase):
    """
    Atoms with occupancy 0 were built by Modeller and have no measured B-factor: they are
    left out of the normalisation and of each residue's average.
    """

    def setUp(self):
        self.outdir = tempfile.mkdtemp(prefix='resdy_test_')
        self.path = os.path.join(self.outdir, 'built.pdb')
        names = ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'CE', 'NZ']
        lines, serial = [], 1
        # LYS 1 measured; LYS 2 with a built side chain; LYS 3 built entirely; GLY 4-6 measured
        for resid, resname in [(1, 'LYS'), (2, 'LYS'), (3, 'LYS'), (4, 'GLY'), (5, 'GLY'), (6, 'GLY')]:
            for k, name in enumerate(names if resname == 'LYS' else names[:4]):
                built = resid == 3 or (resid == 2 and k >= 4)
                occ, b = (0.0, 0.0) if built else (1.0, 10.0 * resid + k)
                lines.append(f'ATOM  {serial:5d}  {name:<3s} {resname} A{resid:4d}    '
                             f'{float(serial):8.3f}{0.0:8.3f}{0.0:8.3f}{occ:6.2f}{b:6.2f}'
                             f'           {name[0]}\n')
                serial += 1
        with open(self.path, 'w') as fh:
            fh.write(''.join(lines) + 'END\n')

    def tearDown(self):
        shutil.rmtree(self.outdir, ignore_errors=True)

    def test_built_atoms_are_ignored(self):
        from resdy.features.FLEXIBILITY import FLEXIBILITY
        import numpy as np
        df = FLEXIBILITY(error_filename='no_record').calculate(self.path).set_index('Resid')

        b = {}
        for line in open(self.path):
            if line.startswith('ATOM') and float(line[54:60]) > 0:
                b.setdefault(int(line[22:26]), []).append(float(line[60:66]))
        measured = np.concatenate([np.array(v) for v in b.values()])
        z = lambda v: (np.array(v) - measured.mean()) / measured.std(ddof=1)

        self.assertAlmostEqual(df.loc[1, 'flexibility'], z(b[1]).mean(), places=10)
        self.assertAlmostEqual(df.loc[2, 'flexibility'], z(b[2]).mean(), places=10)
        self.assertTrue(np.isnan(df.loc[3, 'flexibility']))


class Test_Feature_Settings(unittest.TestCase):
    """
    The geometric constants of the features are settable, validated, and recorded with the
    measurement in measure_settings.json.
    """
    SOURCE = os.path.join('demo', 'conformations', '1A6M-alt1A.pdb')

    def setUp(self):
        self.outdir = tempfile.mkdtemp(prefix='resdy_test_')

    def tearDown(self):
        shutil.rmtree(self.outdir, ignore_errors=True)

    def test_settings_are_recorded_and_checked_on_restart(self):
        import json
        from resdy.measure import MEASURE_SETTINGS_FILE
        df = pd.DataFrame({'PDB_Code': ['1A6M']})
        M = RD.Measure(df_input=df, outdir=self.outdir,
                       features_dict={'sasa': {'probe': 1.2}, 'das': {}})
        settings = M.measurement_settings()
        self.assertEqual(settings['features']['sasa'],
                         {'probe': 1.2, 'n_sphere_point': 960, 'threshold': 0})
        # the radii DAS is given for the residue of interest are recorded, not its None default
        self.assertEqual(settings['features']['das']['radii'], [6.3, 5.9, 5.4, 4.8])
        self.assertEqual(settings['residue_of_interest']['non_modified_codes'][0], 'LYS')

        M._write_settings()
        with open(os.path.join(self.outdir, MEASURE_SETTINGS_FILE)) as fh:
            self.assertEqual(json.load(fh), settings)
        M._check_settings()
        changed = RD.Measure(df_input=df, outdir=self.outdir,
                             features_dict={'sasa': {'probe': 1.4}, 'das': {}})
        with self.assertRaises(ValueError):
            changed._check_settings()

    def test_sasa_cut_out_matches_the_whole_structure(self):
        """At the default probe, the cut-out gives the area computed on the whole structure."""
        import biobox as bb
        import numpy as np
        from resdy.features.SASA import SASA
        df = SASA(error_filename='no_record').calculate(self.SOURCE)
        M = bb.Molecule()
        M.import_pdb(self.SOURCE, include_hetatm=True)
        for _, row in df.iterrows():
            sel = ((M.data['chain'] == row['Chain']) & (M.data['resid'] == row['Resid'])
                   & ~M.data['name'].isin(['CA', 'C', 'N', 'O'])).to_numpy()
            whole = bb.sasa(M, indices=np.where(sel)[0], threshold=0)[0]
            self.assertAlmostEqual(row['sasa'], whole, places=6)

    def test_sasa_cut_out_grows_with_the_probe(self):
        """
        An atom 16 A from NZ, beyond a fixed 15 A cut-out, still occludes part of the side
        chain when the probe is 7 A: the cut-out has to grow with the probe to include it.
        """
        import biobox as bb
        import numpy as np
        from resdy.features.SASA import SASA
        names = ['N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'CE', 'NZ']
        x = [-6.0, -4.5, -4.5, -4.5, -3.0, -1.5, 0.0, 1.5, 3.0]
        y = [0.0, 0.0, 1.5, 2.7, 0.0, 0.0, 0.0, 0.0, 0.0]
        lines = [f'ATOM  {k + 1:5d}  {n:<3s} LYS A   1    {xx:8.3f}{yy:8.3f}{0.0:8.3f}'
                 f'{1.0:6.2f}{20.0:6.2f}           {n[0]}\n' for k, (n, xx, yy) in enumerate(zip(names, x, y))]
        lines.append(f'ATOM     10  CA  GLY A   2    {3.0 + 16.0:8.3f}{0.0:8.3f}{0.0:8.3f}'
                     f'{1.0:6.2f}{20.0:6.2f}           C\n')
        path = os.path.join(self.outdir, 'far_occluder.pdb')
        with open(path, 'w') as fh:
            fh.write(''.join(lines) + 'END\n')

        probe = 7.0
        got = SASA(error_filename='no_record', probe=probe).calculate(path)['sasa'].iloc[0]
        M = bb.Molecule()
        M.import_pdb(path)
        side_chain = np.arange(4, 9)
        whole = bb.sasa(M, indices=side_chain, probe=probe, threshold=0)[0]
        without = bb.sasa(M.get_subset(indices=np.arange(9)), indices=side_chain,
                          probe=probe, threshold=0)[0]
        self.assertLess(whole, without)              # the far atom does occlude
        self.assertAlmostEqual(got, whole, places=6)

    def test_feature_settings_are_validated(self):
        from resdy.features.SASA import SASA
        from resdy.features.DAS import DAS
        from resdy.features.FLEXIBILITY import FLEXIBILITY
        for kwargs in ({'probe': 0}, {'n_sphere_point': 0}, {'threshold': 1.5}):
            with self.assertRaises(ValueError):
                SASA(**kwargs)
        for bad in ({'radii': [5.0, 4.0]}, {'i': 3}, {'pts': 2.0}):
            with self.assertRaises(ValueError):
                DAS(half_sphere_kwargs=bad)
        with self.assertRaises(ValueError):
            FLEXIBILITY(outlier_z=0)

    def test_das_half_sphere_kwargs_are_used(self):
        from resdy.features.DAS import DAS
        coarse = DAS(error_filename='no_record').calculate(self.SOURCE)['das']
        fine = DAS(error_filename='no_record',
                   half_sphere_kwargs={'pts_surf': 2.0}).calculate(self.SOURCE)['das']
        self.assertGreater(fine.sum(), coarse.sum())

    def test_flexibility_raw_b_factors_and_flat_column(self):
        import numpy as np
        from resdy.features.FLEXIBILITY import FLEXIBILITY
        raw = FLEXIBILITY(error_filename='no_record', normalise=False).calculate(self.SOURCE)
        b = {}
        for line in open(self.SOURCE):
            if line.startswith('ATOM') and float(line[54:60]) > 0:
                b.setdefault((line[21], int(line[22:26])), []).append(float(line[60:66]))
        for _, row in raw.head(5).iterrows():
            self.assertAlmostEqual(row['flexibility'], np.mean(b[(row['Chain'], row['Resid'])]), places=6)

        # more than half of the B-factors equal: the median absolute deviation is zero
        flat = os.path.join(self.outdir, 'flat.pdb')
        lines = [l for l in open(self.SOURCE) if l.startswith('ATOM')]
        with open(flat, 'w') as fh:
            for k, l in enumerate(lines):
                fh.write(f'{l[:54]}{1.0:6.2f}{(20.0 if k % 3 else 20.0 + k):6.2f}{l[66:]}')
        out = FLEXIBILITY(error_filename='no_record', remove_outliers=True).calculate(flat)
        self.assertGreater(out['flexibility'].notna().sum(), 0)

    def test_aev_cutoff_beyond_model_radius_changes_nothing(self):
        try:
            from resdy.features.AEV import AEV
            default = AEV(error_filename='no_record')
        except ImportError:
            self.skipTest('torchani is not installed')
        wider = AEV(error_filename='no_record', cutoff=8.0)
        a = default.calculate(self.SOURCE)['aev']
        b = wider.calculate(self.SOURCE)['aev']
        import numpy as np
        from ast import literal_eval
        for x, y in zip(a, b):
            np.testing.assert_allclose(literal_eval(x), literal_eval(y), atol=1e-5)
        with self.assertRaises(ValueError):
            AEV(cutoff=0)


if __name__ == "__main__":
    unittest.main()
