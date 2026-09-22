import unittest
from unittest.mock import Mock, patch
import sys
import os
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(sys.path[0]), "src"))
import resdy as RD

class Test_Uniprot(unittest.TestCase):
    def setUp(self):
        self.df_measures = f'demo{os.sep}measures.csv'
        self.outdir = 'demo'
        self.features_to_analyse = ['depth', 'sasa', 'propka', 'das', 'curvature']
        self.A = RD.Analysis(df=self.df_measures,
                              outdir=self.outdir,
                              features_to_analyse=self.features_to_analyse)

    def test_auxiliary(self):
        df_spec = self.A.get_data(uniprot_entry='P02185',
                                  resid=17)
        self.assertEqual(len(df_spec), 2)

        df_af = self.A.get_data_alphafold()
        self.assertEqual(list(set(list(df_af['Method']))), ['Predicted'])

        df_outlier = self.A.get_outliers(uniprot_entry='P02185',
                                         resid=17,
                                         feature='propka')

        df_extreme = self.A.get_extreme_values(feature='propka',
                                               lower=1,
                                               upper=14)

    def test_remove_unimportant_resids(self):
        df_req_res = pd.read_csv(f'demo{os.sep}demo_input.csv')
        df_cut = self.A.remove_not_important_residues(df_req_res).drop_duplicates(subset=['Uniprot_Entry', 'Resid'])
        self.assertEqual(len(df_cut), 4)

    def test_add_extra_measures(self):
        responses = {'Enter old or new for data to keep: ': 'old'}
        fake_input = Mock(side_effect=responses.get)

        with patch('builtins.input', fake_input):
            self.A.add_extra_measures(extra_measures_filename=self.df_measures,
                                    write_new_file=False)
            self.assertEqual(len(self.A.df.drop_duplicates(subset=['Uniprot_Entry'])), 4)

    def test_GO_analysis(self):
        self.A.GO_get_data()
        df_term = self.A.GO_search_term(df=pd.read_csv(self.df_measures),
                                        code='0031982',
                                        name='')
        go_list = self.A.GO_search_protein(uniprot_entry='P02185')
        df_enrichment = self.A.enrichment_analysis(feature_one = ['propka', 7, 11],
                                                   feature_two = ['sasa', 0, 10],
                                                   uniprot_cnt_cutoff=1)

    def test_relative_best(self):
        resp_base = 'Enter bias towards min or max for '
        resp_end = ' (enter "min" or "max"): '
        responses = {f'{resp_base}curvature{resp_end}': 'min'}
        fake_input = Mock(side_effect=responses.get)

        with patch('builtins.input', fake_input):
            self.A.relative_best(weights=0.5,
                                features=['depth', 'sasa', 'propka', 'das', 'curvature'])
            self.A.relative_best(weights=[0.5, 0.3, 0.1, 0.4],
                                features=['depth', 'sasa', 'propka', 'das', 'curvature'])

    def test_plotting(self):
        self.A.plot_feature_histogram(plot_type='single', feature='depth', save_name='hist_depth.svg')
        self.A.plot_feature_histogram(plot_type='all', agg_type='avg', save_name='hist_all.svg')
        self.A.plot_feature_histogram(plot_type='agg', feature='depth', save_name='hist_agg.svg')
        self.A.plot_feature_violins(features='depth', save_name='violin_depth_allagg.svg')
        self.A.plot_feature_violins(agg_type='all', save_name='violin_all_features_all.svg')



if __name__ == "__main__":
    unittest.main()
