import unittest
from unittest.mock import Mock, patch
import sys
import os
import shutil
import tempfile
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(sys.path[0]), "src"))
import resdy as RD


class Test_Aggregation(unittest.TestCase):
    def setUp(self):
        self.outdir = tempfile.mkdtemp(prefix='resdy_test_')
        shutil.copyfile(os.path.join('demo', 'measures.csv'), os.path.join(self.outdir, 'measures.csv'))
        self.demo_measure = pd.read_csv(f'{self.outdir}{os.sep}measures.csv')

    def tearDown(self):
        shutil.rmtree(self.outdir, ignore_errors=True)

    # TEST SCALAR AGGREGATION

    def test_min(self):
        Agg = RD.Aggregation(df_measurements=self.demo_measure,
                             outdir=self.outdir,
                            aggregation_method='min',
                            features_to_include=['all'],
                            aev_red_method='pca')
        df_agg = Agg.aggregate_data()
        self.assertTrue(len(df_agg) > 0)

    def test_max(self):
        Agg = RD.Aggregation(df_measurements=self.demo_measure,
                             outdir=self.outdir,
                                    aggregation_method='max',
                                    features_to_include=['all'],
                                    aev_red_method='pca')
        df_agg = Agg.aggregate_data()
        self.assertTrue(len(df_agg) > 0)

    def test_avg(self):
        Agg = RD.Aggregation(df_measurements=self.demo_measure,
                             outdir=self.outdir,
                                    aggregation_method='avg',
                                    features_to_include=['all'],
                                    aev_red_method='pca')
        df_agg = Agg.aggregate_data()
        self.assertTrue(len(df_agg) > 0)

    def test_random(self):
        Agg = RD.Aggregation(df_measurements=self.demo_measure,
                             outdir=self.outdir,
                                    aggregation_method='random',
                                    features_to_include=['all'],
                                    aev_red_method='pca')
        df_agg = Agg.aggregate_data()
        self.assertTrue(len(df_agg) > 0)

    def test_median(self):
        Agg = RD.Aggregation(df_measurements=self.demo_measure,
                             outdir=self.outdir,
                                    aggregation_method='median',
                                    features_to_include=['all'],
                                    aev_red_method='pca')
        df_agg = Agg.aggregate_data()
        self.assertTrue(len(df_agg) > 0)

    def test_mixmatch(self):
        Agg = RD.Aggregation(df_measurements=self.demo_measure,
                             outdir=self.outdir,
                                    aggregation_method='mixmatch',
                                    features_to_include=['all'],
                                    aev_red_method='pca')
        df_agg = Agg.aggregate_data()
        self.assertTrue(len(df_agg) > 0)

    def test_minmax(self):
        Agg = RD.Aggregation(df_measurements=self.demo_measure,
                             outdir=self.outdir,
                                    aggregation_method='minmax',
                                    features_to_include=['all'],
                                    aev_red_method='pca')
        df_agg = Agg.aggregate_data()
        self.assertTrue(len(df_agg) > 0)

    def test_minmaxavg(self):
        Agg = RD.Aggregation(df_measurements=self.demo_measure,
                             outdir=self.outdir,
                                    aggregation_method='minmaxavg',
                                    features_to_include=['all'],
                                    aev_red_method='pca')
        df_agg = Agg.aggregate_data()
        self.assertTrue(len(df_agg) > 0)
        Agg.save_state(outname=f'measures_aggregated.csv')

    def test_all(self):
        Agg = RD.Aggregation(df_measurements=self.demo_measure,
                             outdir=self.outdir,
                                    aggregation_method='all',
                                    features_to_include=['all'],
                                    aev_red_method='pca')
        df_agg = Agg.aggregate_data()
        self.assertTrue(len(df_agg) > 0)

    def test_choose(self):
        resp_base = '>> Choose statistics for feature: '
        responses = {f'{resp_base}propka': '1',
                     f'{resp_base}sasa': '2',
                     f'{resp_base}depth': '3',
                     f'{resp_base}seqcharge': '4',
                     f'{resp_base}curvature': '5',
                     f'{resp_base}writhing': '6',
                     f'{resp_base}torsion': '7',
                     f'{resp_base}arc_length': '8',
                     f'{resp_base}phi': '1;2',
                     f'{resp_base}psi': '9',
                     f'{resp_base}das': 'all',
                     f'{resp_base}flexibility': 'max',
                     f'{resp_base}rmsf': 'med'}

        fake_input = Mock(side_effect=responses.get)
        expected_columns = ['Uniprot_Entry', 'Resid', 'class', 'propka_min', 'propka_max',
                            'propka_avg', 'propka_sd', 'propka_range', 'propka_rand', 'sasa_max',
                            'depth_min', 'seqcharge_med', 'curvature_avg', 'writhing_sd',
                            'torsion_range', 'phi_min', 'phi_max', 'phi_avg', 'phi_sd', 'phi_range',
                            'phi_rand', 'psi_min', 'psi_max', 'psi_med', 'psi_avg', 'psi_sd',
                            'psi_range', 'psi_rand', 'das_min', 'das_max', 'das_avg', 'das_sd',
                            'das_range', 'das_rand', 'flexibility_max', 'rmsf_med']

        # keep every feature whatever its NaN fraction (rmsf is NaN for every structure with a
        # single model), so that each of the statistics asked for above is exercised
        with patch('builtins.input', fake_input):
            Agg = RD.Aggregation(df_measurements=self.demo_measure,
                                 outdir=self.outdir,
                                        aggregation_method='choose',
                                        features_to_include=['all'],
                                        aev_red_method='pca',
                                        max_feature_nan_fraction=1.0)
            df_agg = Agg.aggregate_data()
            self.assertTrue(len(df_agg) > 0)
            self.assertListEqual(expected_columns, list(df_agg.columns))

    def test_unknown(self):
        Agg = RD.Aggregation(df_measurements=self.demo_measure,
                             outdir=self.outdir,
                                        aggregation_method='x',
                                        features_to_include=['all'],
                                        aev_red_method='pca')
        df_agg = Agg.aggregate_data()
        self.assertTrue(len(df_agg) > 0)
        self.assertEqual(Agg.aggregation_method, 'minmax')


    # TEST VECTOR FEATURES AGGREGATION
    def test_pca(self):
        Agg = RD.Aggregation(df_measurements=self.demo_measure,
                             outdir=self.outdir,
                            aggregation_method='minmax',
                            features_to_include=['all'],
                            aev_red_method='pca')
        df_agg = Agg.aggregate_data()
        self.assertTrue(len(df_agg) > 0)

    def test_sd(self):
        Agg = RD.Aggregation(df_measurements=self.demo_measure,
                             outdir=self.outdir,
                            aggregation_method='minmax',
                            features_to_include=['all'],
                            aev_red_method='sd')
        df_agg = Agg.aggregate_data()
        self.assertTrue(len(df_agg) > 0)

    def test_null(self):
        Agg = RD.Aggregation(df_measurements=self.demo_measure,
                             outdir=self.outdir,
                            aggregation_method='minmax',
                            features_to_include=['all'],
                            aev_red_method='null')
        df_agg = Agg.aggregate_data()
        self.assertTrue(len(df_agg) > 0)


class Test_Aggregation_Settings(unittest.TestCase):
    """aev_pca_variance and vif_threshold, which used to be fixed at 0.99 and 5."""

    def setUp(self):
        self.outdir = tempfile.mkdtemp(prefix='resdy_test_')
        self.demo_measure = pd.read_csv(os.path.join('demo', 'measures.csv'))

    def tearDown(self):
        shutil.rmtree(self.outdir, ignore_errors=True)

    def _agg(self, **kwargs):
        return RD.Aggregation(df_measurements=self.demo_measure.copy(), outdir=self.outdir,
                              aggregation_method='min', features_to_include=['all'],
                              aev_red_method='pca', **kwargs)

    def test_pca_variance_fraction_and_component_count(self):
        import numpy as np
        n_aev = lambda agg: len([c for c in agg.df_measurements.columns if c.startswith('AEV_')])
        few = self._agg(aev_pca_variance=0.5)
        few.aggregate_data()
        many = self._agg(aev_pca_variance=0.999)
        many.aggregate_data()
        self.assertLess(n_aev(few), n_aev(many))
        fixed = self._agg(aev_pca_variance=np.int64(3))
        fixed.aggregate_data()
        self.assertEqual(n_aev(fixed), 3)

    def test_invalid_settings_fail_at_construction(self):
        for bad in (0.0, 1.0, 1.5, 0, True, 'all'):
            with self.assertRaises(ValueError):
                self._agg(aev_pca_variance=bad)
        with self.assertRaises(ValueError):
            self._agg(vif_threshold=1.0)

    def test_vif_threshold_reaches_preprocessing(self):
        with patch('resdy.aggregation.Preprocessing') as P:
            P.return_value.calculate_diff_features.return_value = []
            agg = RD.Aggregation(df_measurements=self.demo_measure.copy(), outdir=self.outdir,
                                 aggregation_method='min', features_to_include=['all'],
                                 aev_red_method='vif', vif_threshold=7.5)
            agg.aggregate_data()
        self.assertEqual(P.return_value.calculate_diff_features.call_args.kwargs['vif_threshold'], 7.5)


if __name__ == "__main__":
    unittest.main()
