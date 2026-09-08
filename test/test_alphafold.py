import unittest
import sys
import os
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(sys.path[0]), "src"))
from coolpackage import alphafold as af

class Test_Uniprot(unittest.TestCase):
    def test_download_af(self):
        af.download_AF_struc('AF-P0CG48-F1-model_v6', outfolder='test')
        self.assertTrue(os.path.isfile(f'test{os.sep}curated{os.sep}AF-P0CG48-F1-model_v6.pdb'))

    def test_plddt(self):
        af.find_af_plddt('AF-P0CG48-F1-model_v6', outfolder='test')
        self.assertTrue(os.path.isfile(f'test{os.sep}curated{os.sep}AF_PLDDT_Output.csv'))
        df_plddt = pd.read_csv(f'test{os.sep}curated{os.sep}AF_PLDDT_Output.csv')
        self.assertTrue(len(df_plddt['PDB_Code']) > 0)


if __name__ == "__main__":
    unittest.main()