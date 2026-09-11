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
        M = CPN.Measure(df_input=self.df_prot,
                        outdir=self.outdir,
                        features=['propka'],
                        parallel=True,
                        include_modified=False)
        M.measure_data()
        print('len propka', len(M.df))
        self.assertTrue(len(M.df) > 0)

    def test_sasa(self):
        M = CPN.Measure(df_input=self.df_prot,
                        outdir=self.outdir,
                        features=['sasa'],
                        parallel=True,
                        include_modified=False)
        M.measure_data()
        print('len sasa', len(M.df))
        self.assertTrue(len(M.df) > 0)

    def test_depth(self):
        M = CPN.Measure(df_input=self.df_prot,
                        outdir=self.outdir,
                        features=['depth'],
                        parallel=True,
                        include_modified=False)
        M.measure_data()
        print('len depth', len(M.df))
        self.assertTrue(len(M.df) > 0)

    def test_aev(self):
        M = CPN.Measure(df_input=self.df_prot,
                        outdir=self.outdir,
                        features=['aev'],
                        parallel=True,
                        include_modified=False)
        M.measure_data()
        print('len aev', len(M.df))
        self.assertTrue(len(M.df) > 0)

    def test_das(self):
        M = CPN.Measure(df_input=self.df_prot,
                        outdir=self.outdir,
                        features=['das'],
                        parallel=True,
                        include_modified=False)
        M.measure_data()
        print('len das', len(M.df))
        self.assertTrue(len(M.df) > 0)

    def test_seqcharge(self):
        M = CPN.Measure(df_input=self.df_prot,
                        outdir=self.outdir,
                        features=['seqcharge'],
                        parallel=True,
                        include_modified=False)
        M.measure_data()
        print('len seqcharge', len(M.df))
        self.assertTrue(len(M.df) > 0)

    def test_flexibility(self):
        M = CPN.Measure(df_input=self.df_prot,
                        outdir=self.outdir,
                        features=['flexibility'],
                        parallel=True,
                        include_modified=False)
        M.measure_data()
        print('len flexibility', len(M.df))
        self.assertTrue(len(M.df) > 0)

    def test_legolas(self):
        M = CPN.Measure(df_input=self.df_prot,
                        outdir=self.outdir,
                        features=['legolas'],
                        parallel=True,
                        include_modified=False)
        M.measure_data()
        print('len legolas', len(M.df))

    def test_frustration(self):
        M = CPN.Measure(df_input=self.df_prot,
                        outdir=self.outdir,
                        features=['frustration', 'density'],
                        parallel=True,
                        include_modified=False)
        M.measure_data()
        print('len frustration', len(M.df))

    def test_melodia(self):
        M = CPN.Measure(df_input=self.df_prot,
                        outdir=self.outdir,
                        features=['melodia'],
                        parallel=True,
                        include_modified=False)
        M.measure_data()
        print('len propka', len(M.df))
        self.assertTrue(len(M.df) > 0)

    def test_phi(self):
        M = CPN.Measure(df_input=self.df_prot,
                        outdir=self.outdir,
                        features=['phi'],
                        parallel=True,
                        include_modified=False)
        M.measure_data()
        print('len phi', len(M.df))
        self.assertTrue(len(M.df) > 0)

    def test_evolution(self):
        M = CPN.Measure(df_input=self.df_prot,
                        outdir=self.outdir,
                        features=['evolution'],
                        parallel=True,
                        include_modified=False)
        M.measure_data()
        print('len evolution', len(M.df))
        self.assertTrue(len(M.df) > 0)

    def test_rmsf(self):
        M = CPN.Measure(df_input=self.df_prot,
                        outdir=self.outdir,
                        features=['rmsf'],
                        parallel=True,
                        include_modified=False)
        M.measure_data()
        print('len rmsf', len(M.df))
        self.assertTrue(len(M.df) > 0)

    def test_pkaani(self):
        M = CPN.Measure(df_input=self.df_prot,
                        outdir=self.outdir,
                        features=['pkaANI'],
                        parallel=True,
                        include_modified=False)
        M.measure_data()
        print('len pkaani', len(M.df))

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

# ADD IN PDB ONLY TESTS


if __name__ == "__main__":
    unittest.main()