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

    def tearDown(self):
        shutil.rmtree(self.outdir, ignore_errors=True)

    def test_full_setup(self):
        # test full setup of measures class, no run, test runs individually
        print('-> Testing full setup')
        M = RD.Measure(df_input=self.df_prot,
                        outdir=self.outdir,
                        features_dict=self.features_dict)
        reg_feats = {m[0] for m in M.measures}
        dereg_feats = set(self.features_dict) - reg_feats - set(M.features)
        self.assertEqual(dereg_feats, set(), f'> Some features vanished without being removed: {dereg_feats}')

    def test_propka(self):
        print('-> Test measuring PROPKA')
        M_propka = RD.Measure(df_input=self.df_prot[:2],
                        outdir=self.outdir,
                        features_dict={'propka': {}},
                        parallel=True,
                        include_modified=False)
        M_propka.measure_data()
        print(M_propka.df)
        print('len propka', len(M_propka.df))
        self.assertTrue(len(M_propka.df) > 0)

    def test_sasa(self):
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

    def test_depth(self):
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

    def test_aev(self):
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
        import pandas as pd
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


if __name__ == "__main__":
    unittest.main()
