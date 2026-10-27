import unittest
import sys
import os
import tempfile
import shutil
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(sys.path[0]), "src"))
import resdy as RD

class Test_Uniprot(unittest.TestCase):
    def setUp(self):
        self.outdir = tempfile.mkdtemp(prefix='resdy_test_')
        shutil.copyfile(os.path.join('demo', 'measures.csv'), os.path.join(self.outdir, 'measures.csv'))
        self.UP = RD.Uniprot(outdir=self.outdir)
        self.demo_csv = f'demo{os.sep}demo_input.csv'

    def test_organism_funcs(self):
        tot_review = self.UP.count_organism_proteins(code='UP000007445', reviewed_only=True)
        tot_all = self.UP.count_organism_proteins(code='UP000007445', reviewed_only=False)
        self.assertEqual(int(tot_review), 4)
        self.assertEqual(int(tot_all), 5)

        codes_review = self.UP.get_organism_proteins(code='UP000007445', reviewed_only=True)
        codes_all = self.UP.get_organism_proteins(code='UP000007445', reviewed_only=False)
        self.assertEqual(len(codes_review), 4)
        self.assertEqual(len(codes_all), 5)

        self.UP.from_organism(code='UP000007445', reviewed_only=True)
        self.assertEqual(len(self.UP.df), 4)
        self.UP.from_organism(code='UP000007445', reviewed_only=False)
        self.assertEqual(len(self.UP.df), 5)

    def test_protein(self):
        self.UP.get_protein_data('P09167')
        self.UP.get_protein_data('P0CG48', '2MWS')

        self.UP.from_csv_file(csv_file=self.demo_csv)

    def test_auxiliary(self):
        self.UP.from_csv_file(csv_file=self.demo_csv)

        self.UP.filter_by_technique(list_of_techniques=['X-ray'])
        self.UP.filter_by_technique(list_of_techniques=['NMR'])
        self.UP.filter_by_technique(list_of_techniques=['Predicted'])

        #self.UP.save_state()

    def test_GO_data(self):
        test_uniprots = ['P40616', 'P02185', 'Q16825', 'P0CG48', 'P10724']
        self.UP.GO_get_data(uniprot=test_uniprots)
        test_df = pd.read_csv(f'{self.outdir}{os.sep}measures.csv')
        df_term = self.UP.GO_search_term(df=test_df,
                                        code='0031982',
                                        name='')
        self.assertTrue(len(df_term) > 0)
        go_list = self.UP.GO_search_protein(uniprot_entry='P02185')
        self.assertTrue(len(go_list) > 0)
        self.UP.save_GO_metadata()
        self.assertTrue(os.path.isfile(os.path.join(self.outdir, 'GO_metadata.json')))
        self.UP.load_GO_metadata()
        reloaded_go_list = self.UP.GO_search_protein(uniprot_entry='P02185')
        self.assertEqual(set(go_list), set(reloaded_go_list))

if __name__ == "__main__":
    unittest.main()
