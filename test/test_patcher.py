import unittest
import sys
import os
import shutil
import tempfile

sys.path.insert(0, os.path.join(os.path.dirname(sys.path[0]), "src"))
from resdy import patcher

class Test_Patcher(unittest.TestCase):
    def setUp(self):
        self.outdir = tempfile.mkdtemp(prefix='resdy_test_')
        os.makedirs(os.path.join(self.outdir, 'conformations'), exist_ok=True)
        os.makedirs(os.path.join(self.outdir, 'curated'), exist_ok=True)
        shutil.copyfile(os.path.join('demo', 'conformations', '2MWS-alt-1.pdb'), os.path.join(self.outdir, 'conformations', '2MWS-alt-1.pdb'))
        shutil.copyfile(os.path.join('demo', 'conformations', '2MWS.fasta'), os.path.join(self.outdir, 'conformations', '2MWS.fasta'))

    def tearDown(self):
        shutil.rmtree(self.outdir, ignore_errors=True)

    def test_patcher_pipeline(self):
        pdb = f"{self.outdir}{os.sep}conformations{os.sep}2MWS-alt-1.pdb"
        fasta = f"{self.outdir}{os.sep}conformations{os.sep}2MWS.fasta"
        gap = 10
        fname = patcher.curate(pdb=pdb, fasta=fasta, outdir=self.outdir, gap=gap)
        self.assertTrue(os.path.isfile(f'{self.outdir}{os.sep}curated{os.sep}2MWS-alt-1.pdb'))


if __name__ == "__main__":
    unittest.main()
