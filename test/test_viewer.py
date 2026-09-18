import unittest
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(sys.path[0]), "src"))
import resdy as RD

class Test_Uniprot(unittest.TestCase):
    def setUp(self):
        self.outdir = 'demo'
        self.df_input = f'{self.outdir}{os.sep}demo_input.csv'
        self.df_measures = f'{self.outdir}{os.sep}measures.csv'
        self.A = RD.Analysis(df=self.df_measures,
                             outdir=self.outdir)
        self.V = RD.Viewer(outdir=self.outdir,
                           analysis=self.A)


if __name__ == "__main__":
    unittest.main()