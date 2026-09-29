import unittest
import sys
import os
import shutil
import tempfile
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(sys.path[0]), "src"))
from resdy import alphafold as af

class Test_Alphafold(unittest.TestCase):
    def setUp(self):
        self.outdir = tempfile.mkdtemp(prefix='resdy_test_')

    def tearDown(self):
        shutil.rmtree(self.outdir, ignore_errors=True)

    def test_af(self):
        af.download_AF_struc('AF-P0CG48-F1-model_v6', outfolder=self.outdir)
        self.assertTrue(os.path.isfile(f'{self.outdir}{os.sep}curated{os.sep}AF-P0CG48-F1-model_v6.pdb'))

        af.find_af_plddt('AF-P0CG48-F1-model_v6', outfolder=self.outdir)
        self.assertTrue(os.path.isfile(f'{self.outdir}{os.sep}curated{os.sep}AF_PLDDT_Output.csv'))
        df_plddt = pd.read_csv(f'{self.outdir}{os.sep}curated{os.sep}AF_PLDDT_Output.csv')
        self.assertTrue(len(df_plddt['PDB_Code']) > 0)

if __name__ == "__main__":
    unittest.main()
