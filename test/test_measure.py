import unittest
import sys
import os
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(sys.path[0]), "src"))
import coolpackage as CPN

class Test_Uniprot(unittest.TestCase):
    def setUp(self):
        self.outdir = 'demo'
        self.df_input = pd.read_csv(f'{self.outdir}{os.sep}demo_input.csv')
        self.df_prot = pd.read_csv(f'{self.outdir}{os.sep}proteins.csv')
        self.all_features = ['propka', 'pkaANI', 'sasa', 'depth', 'aev', 'das', 'seqcharge', 'melodia', 'frustration']

    def test_full_setup(self):
        # test full setup of measures class, no run, test runs individually
        M = CPN.Measure(df_input=self.df_prot,
                        outdir=self.outdir,
                        features=self.all_features)

    def test_propka(self):
        M_propka = CPN.Measure(df_input=self.df_prot,
                        outdir=self.outdir,
                        features=['propka'],
                        parallel=True,
                        include_modified=False)
        M_propka.measure_data()
        print(M_propka.df)
        print('len propka', len(M_propka.df))
        self.assertTrue(len(M_propka.df) > 0)

    def test_sasa(self):
        M_sasa = CPN.Measure(df_input=self.df_prot,
                        outdir=self.outdir,
                        features=['sasa'],
                        parallel=True,
                        include_modified=False)
        M_sasa.measure_data()
        print(M_sasa.df)
        print('len sasa', len(M_sasa.df))
        self.assertTrue(len(M_sasa.df) > 0)

    def test_depth(self):
        M_depth = CPN.Measure(df_input=self.df_prot,
                        outdir=self.outdir,
                        features=['depth'],
                        parallel=True,
                        include_modified=False)
        M_depth.measure_data()
        print(M_depth.df)
        print('len depth', len(M_depth.df))
        self.assertTrue(len(M_depth.df) > 0)

    def test_aev(self):
        M_aev = CPN.Measure(df_input=self.df_prot,
                        outdir=self.outdir,
                        features=['aev'],
                        parallel=True,
                        include_modified=False)
        M_aev.measure_data()
        print(M_aev.df)
        print('len aev', len(M_aev.df))
        self.assertTrue(len(M_aev.df) > 0)

    def test_das(self):
        M_das = CPN.Measure(df_input=self.df_prot,
                        outdir=self.outdir,
                        features=['das'],
                        parallel=True,
                        include_modified=False)
        M_das.measure_data()
        print(M_das.df)
        print('len das', len(M_das.df))
        self.assertTrue(len(M_das.df) > 0)

    def test_seqcharge(self):
        M_seqcharge = CPN.Measure(df_input=self.df_prot,
                        outdir=self.outdir,
                        features=['seqcharge'],
                        parallel=True,
                        include_modified=False)
        M_seqcharge.measure_data()
        print(M_seqcharge.df)
        print('len seqcharge', len(M_seqcharge.df))
        self.assertTrue(len(M_seqcharge.df) > 0)

    def test_flexibility(self):
        M_flex = CPN.Measure(df_input=self.df_prot,
                        outdir=self.outdir,
                        features=['flexibility'],
                        parallel=True,
                        include_modified=False)
        M_flex.measure_data()
        print(M_flex.df)
        print('len flexibility', len(M_flex.df))
        self.assertTrue(len(M_flex.df) > 0)

    def test_legolas(self):
        M_nmr = CPN.Measure(df_input=self.df_prot,
                        outdir=self.outdir,
                        features=['legolas'],
                        parallel=True,
                        include_modified=False)
        M_nmr.measure_data()
        print(M_nmr.df)
        print('len legolas:', len(M_nmr.df))

    def test_frustration(self):
        M_frustration = CPN.Measure(df_input=self.df_prot,
                        outdir=self.outdir,
                        features=['frustration', 'density'],
                        parallel=False,
                        include_modified=False)
        M_frustration.measure_data()
        print(M_frustration.df)
        print('len frustration', len(M_frustration.df))

    def test_melodia(self):
        M_melodia = CPN.Measure(df_input=self.df_prot,
                        outdir=self.outdir,
                        features=['melodia'],
                        parallel=True,
                        include_modified=False)
        M_melodia.measure_data()
        print(M_melodia.df)
        print('len propka', len(M_melodia.df))
        self.assertTrue(len(M_melodia.df) > 0)

    def test_phi(self):
        M_phi = CPN.Measure(df_input=self.df_prot,
                        outdir=self.outdir,
                        features=['phi'],
                        parallel=True,
                        include_modified=False)
        M_phi.measure_data()
        print(M_phi.df)
        print('len phi', len(M_phi.df))
        self.assertTrue(len(M_phi.df) > 0)

    def test_evolution(self):
        M_evolution = CPN.Measure(df_input=self.df_prot,
                        outdir=self.outdir,
                        features=['evolution'],
                        parallel=True,
                        include_modified=False)
        M_evolution.measure_data()
        print(M_evolution.df)
        print('len evolution', len(M_evolution.df))
        self.assertTrue(len(M_evolution.df) > 0)

    def test_rmsf(self):
        M_rmsf = CPN.Measure(df_input=self.df_prot,
                        outdir=self.outdir,
                        features=['rmsf'],
                        parallel=True,
                        include_modified=False)
        M_rmsf.measure_data()
        print(M_rmsf.df)
        print('len rmsf', len(M_rmsf.df))
        self.assertTrue(len(M_rmsf.df) > 0)

    def test_pkaani(self):
        M_pkaani = CPN.Measure(df_input=self.df_prot,
                        outdir=self.outdir,
                        features=['pkaANI'],
                        parallel=True,
                        include_modified=False)
        M_pkaani.measure_data()
        print(M_pkaani.df)
        print('len pkaani', len(M_pkaani.df))

    def test_secondarystructure(self):
        M_ss = CPN.Measure(df_input=self.df_prot,
                        outdir=self.outdir,
                        features=['SECONDARYSTRUCTURE'],
                        parallel=True,
                        include_modified=False)
        M_ss.measure_data()
        print(M_ss.df)
        print('len secondarystructure', len(M_ss.df))

    '''
    def test_restart_measure_data(self):
        M = CPN.Measure(df_input=self.df_prot,
                        outdir=self.outdir,
                        features=['das'],
                        parallel=True,
                        include_modified=False)
        M.restart_measure_data()

    
    def test_recovery(self):
        M = CPN.Measure(df_input=self.df_prot,
                        outdir=self.outdir,
                        features=['das'],
                        parallel=True,
                        include_modified=False)
        M.recover_from_log('measure_log.txt')
    '''



if __name__ == "__main__":
    unittest.main()