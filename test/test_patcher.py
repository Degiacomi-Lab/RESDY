import unittest
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(sys.path[0]), "src"))
from coolpackage import patcher

class Test_Uniprot(unittest.TestCase):
    def test_patcher_pipeline(self):
        pdb = f"demo{os.sep}conformations{os.sep}2MWS-alt-1.pdb"
        fasta = f"demo{os.sep}conformations{os.sep}2MWS.fasta"
        outdir = f'test{os.sep}curated'
        gap = 10
        fname = patcher.curate(pdb=pdb, fasta=fasta, outdir=outdir, gap=gap)
        self.assertTrue(os.path.isfile(f'test{os.sep}curated{os.sep}2MWS-alt-1.pdb'))


if __name__ == "__main__":
    unittest.main()