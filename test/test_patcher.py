import unittest
import sys
import os
import shutil
import tempfile
import itertools

import numpy as np

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

    @staticmethod
    def _chain_centroids(path):
        '''
        :returns: dict mapping chain name to the centroid of its CA atoms.
        '''
        pts = {}
        for line in open(path):
            if line.startswith('ATOM') and line[12:16].strip() == 'CA':
                pts.setdefault(line[21], []).append(
                    (float(line[30:38]), float(line[38:46]), float(line[46:54])))
        return {c: np.mean(v, axis=0) for c, v in pts.items()}

    def test_patcher_pipeline(self):
        pdb = f"{self.outdir}{os.sep}conformations{os.sep}2MWS-alt-1.pdb"
        fasta = f"{self.outdir}{os.sep}conformations{os.sep}2MWS.fasta"
        gap = 10
        outname, largest = patcher.curate(pdb=pdb, fasta=fasta, outdir=self.outdir, gap=gap)
        self.assertTrue(os.path.isfile(outname))
        self.assertEqual(largest, 0)

    def test_chains_keep_their_relative_placement(self):
        '''
        Curation must not move the chains of a complex with respect to one another.
        Modeller returns each patched chain in its own frame, so reassembling the chains
        without putting each one back destroys the quaternary structure.
        '''
        pdb = f"{self.outdir}{os.sep}conformations{os.sep}2MWS-alt-1.pdb"
        fasta = f"{self.outdir}{os.sep}conformations{os.sep}2MWS.fasta"
        before = self._chain_centroids(pdb)
        outname, _ = patcher.curate(pdb=pdb, fasta=fasta, outdir=self.outdir, gap=10)
        after = self._chain_centroids(outname)

        self.assertEqual(sorted(before), sorted(after))
        chains = sorted(before)
        for a, b in itertools.combinations(chains, 2):
            sep_before = np.linalg.norm(before[a] - before[b])
            sep_after = np.linalg.norm(after[a] - after[b])
            self.assertAlmostEqual(
                sep_before, sep_after, delta=1.0,
                msg=f'chains {a} and {b} moved from {sep_before:.2f} A apart to '
                    f'{sep_after:.2f} A during curation')

    def test_superpose_onto_restores_a_displaced_chain(self):
        '''
        superpose_onto must undo a rigid displacement exactly.
        '''
        src = os.path.join(self.outdir, 'conformations', '2MWS-alt-1.pdb')
        moved = os.path.join(self.outdir, 'moved.pdb')
        shift = np.array([12.0, -7.5, 3.25])
        with open(src) as fin, open(moved, 'w') as fout:
            for line in fin:
                if line.startswith(('ATOM', 'HETATM')):
                    v = np.array([float(line[30:38]), float(line[38:46]),
                                  float(line[46:54])]) + shift
                    line = f'{line[:30]}{v[0]:8.3f}{v[1]:8.3f}{v[2]:8.3f}{line[54:]}'
                fout.write(line)

        n_fit, rmsd = patcher.superpose_onto(moved, src)
        self.assertGreater(n_fit, 100)
        self.assertLess(rmsd, 1e-3)
        recovered = self._chain_centroids(moved)
        original = self._chain_centroids(src)
        for c in original:
            self.assertLess(float(np.linalg.norm(recovered[c] - original[c])), 1e-3)


if __name__ == "__main__":
    unittest.main()
