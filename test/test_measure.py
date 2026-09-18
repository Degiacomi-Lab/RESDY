import unittest
import sys
import os
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(sys.path[0]), "src"))
import resdy as RD

class Test_Uniprot(unittest.TestCase):
    def setUp(self):
        self.outdir = 'demo'
        self.df_input = pd.read_csv(f'{self.outdir}{os.sep}demo_input.csv')
        self.df_prot = pd.read_csv(f'{self.outdir}{os.sep}proteins.csv')
        self.all_features = ['propka', 'pkaANI', 'sasa', 'depth', 'aev',
                             'das', 'seqcharge', 'melodia', 'frustration']
        self.features_dict = {k: {} for k in self.all_features}

    def test_full_setup(self):
        # test full setup of measures class, no run, test runs individually
        M = RD.Measure(df_input=self.df_prot,
                        outdir=self.outdir,
                        features_dict=self.features_dict)
        self.assertTrue(set(M.features) <= set(self.all_features))
        self.assertTrue(len(M.measures) > 0)

    def test_propka(self):
        M_propka = RD.Measure(df_input=self.df_prot,
                        outdir=self.outdir,
                        features_dict={'propka': {}},
                        parallel=True,
                        include_modified=False)
        M_propka.measure_data()
        print(M_propka.df)
        print('len propka', len(M_propka.df))
        self.assertTrue(len(M_propka.df) > 0)

    def test_sasa(self):
        M_sasa = RD.Measure(df_input=self.df_prot,
                        outdir=self.outdir,
                        features_dict={'sasa': {}},
                        parallel=True,
                        include_modified=False)
        M_sasa.measure_data()
        print(M_sasa.df)
        print('len sasa', len(M_sasa.df))
        self.assertTrue(len(M_sasa.df) > 0)

    def test_depth(self):
        M_depth = RD.Measure(df_input=self.df_prot,
                        outdir=self.outdir,
                        features_dict={'depth': {}},
                        parallel=True,
                        include_modified=False)
        M_depth.measure_data()
        print(M_depth.df)
        print('len depth', len(M_depth.df))
        self.assertTrue(len(M_depth.df) > 0)

    def test_aev(self):
        M_aev = RD.Measure(df_input=self.df_prot,
                        outdir=self.outdir,
                        features_dict={'aev': {}},
                        parallel=True,
                        include_modified=False)
        M_aev.measure_data()
        print(M_aev.df)
        print('len aev', len(M_aev.df))
        self.assertTrue(len(M_aev.df) > 0)

    def test_das(self):
        M_das = RD.Measure(df_input=self.df_prot,
                        outdir=self.outdir,
                        features_dict={'das': {}},
                        parallel=True,
                        include_modified=False)
        M_das.measure_data()
        print(M_das.df)
        print('len das', len(M_das.df))
        self.assertTrue(len(M_das.df) > 0)

    def test_charge(self):
        M_seqcharge = RD.Measure(df_input=self.df_prot,
                        outdir=self.outdir,
                        features_dict={'seqcharge': {}},
                        parallel=True,
                        include_modified=False)
        M_seqcharge.measure_data()
        print(M_seqcharge.df)
        print('len seqcharge', len(M_seqcharge.df))
        self.assertTrue(len(M_seqcharge.df) > 0)

    def test_flexibility(self):
        M_flex = RD.Measure(df_input=self.df_prot,
                        outdir=self.outdir,
                        features_dict={'flexibility': {}},
                        parallel=True,
                        include_modified=False)
        M_flex.measure_data()
        print(M_flex.df)
        print('len flexibility', len(M_flex.df))
        self.assertTrue(len(M_flex.df) > 0)

    def test_legolas(self):
        M_legolas = RD.Measure(df_input=self.df_prot,
                        outdir=self.outdir,
                        features_dict={'legolas': {}},
                        parallel=True,
                        include_modified=False)
        M_legolas.measure_data()
        print(M_legolas.df)
        print('len legolas:', len(M_legolas.df))

    def test_frustration(self):
        M_frustration = RD.Measure(df_input=self.df_prot,
                        outdir=self.outdir,
                        features_dict={'frustration': {}, 'density': {}},
                        parallel=False,
                        include_modified=False)
        M_frustration.measure_data()
        print(M_frustration.df)
        print('len frustration', len(M_frustration.df))

    def test_melodia(self):
        M_melodia = RD.Measure(df_input=self.df_prot,
                        outdir=self.outdir,
                        features_dict={'melodia': {}},
                        parallel=True,
                        include_modified=False)
        M_melodia.measure_data()
        print(M_melodia.df)
        print('len propka', len(M_melodia.df))
        self.assertTrue(len(M_melodia.df) > 0)

    def test_phi(self):
        M_phi = RD.Measure(df_input=self.df_prot,
                        outdir=self.outdir,
                        features_dict={'phi': {}},
                        parallel=True,
                        include_modified=False)
        M_phi.measure_data()
        print(M_phi.df)
        print('len phi', len(M_phi.df))
        self.assertTrue(len(M_phi.df) > 0)

    def test_evolution(self):
        M_evolution = RD.Measure(df_input=self.df_prot,
                        outdir=self.outdir,
                        features_dict={'evolution': {}},
                        parallel=True,
                        include_modified=False)
        M_evolution.measure_data()
        print(M_evolution.df)
        print('len evolution', len(M_evolution.df))
        self.assertTrue(len(M_evolution.df) > 0)

    def test_rmsf(self):
        M_rmsf = RD.Measure(df_input=self.df_prot,
                        outdir=self.outdir,
                        features_dict={'rmsf': {}},
                        parallel=True,
                        include_modified=False)
        M_rmsf.measure_data()
        print(M_rmsf.df)
        print('len rmsf', len(M_rmsf.df))
        self.assertTrue(len(M_rmsf.df) > 0)

    def test_pkaani(self):
        M_pkaani = RD.Measure(df_input=self.df_prot,
                        outdir=self.outdir,
                        features_dict={'pkaani': {}},
                        parallel=True,
                        include_modified=False)
        M_pkaani.measure_data()
        print(M_pkaani.df)
        print('len pkaani', len(M_pkaani.df))

    def test_secondarystructure(self):
        M_ss = RD.Measure(df_input=self.df_prot,
                        outdir=self.outdir,
                        features_dict={'secondarystructure': {}},
                        parallel=True,
                        include_modified=False)
        M_ss.measure_data()
        print(M_ss.df)
        print('len secondarystructure', len(M_ss.df))

    '''
    def test_restart_measure_data(self):
        M = RD.Measure(df_input=self.df_prot,
                        outdir=self.outdir,
                        features_dict={'das': {}},
                        parallel=True,
                        include_modified=False)
        M.restart_measure_data()

    
    def test_recovery(self):
        M = RD.Measure(df_input=self.df_prot,
                        outdir=self.outdir,
                        features_dict={'das': {}},
                        parallel=True,
                        include_modified=False)
        M.recover_from_log('measure_log.txt')
    '''



if __name__ == "__main__":
    unittest.main()