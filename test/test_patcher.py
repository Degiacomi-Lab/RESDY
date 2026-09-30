import unittest
import sys
import os
import shutil
import tempfile
import itertools

import numpy as np
import biobox as bb

sys.path.insert(0, os.path.join(os.path.dirname(sys.path[0]), "src"))
from resdy import patcher
from resdy.geometry import check_geometry

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
        outname, largest, geometry = patcher.curate(pdb=pdb, fasta=fasta, outdir=self.outdir, gap=gap)
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
        outname, _, _ = patcher.curate(pdb=pdb, fasta=fasta, outdir=self.outdir, gap=10)
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

    def test_curate_reports_geometry(self):
        '''curate must hand back a geometry report for the assembled structure.'''
        pdb = f"{self.outdir}{os.sep}conformations{os.sep}2MWS-alt-1.pdb"
        fasta = f"{self.outdir}{os.sep}conformations{os.sep}2MWS.fasta"
        outname, largest, geometry = patcher.curate(pdb=pdb, fasta=fasta,
                                                    outdir=self.outdir, gap=10)
        for key in ('n_atoms', 'n_clashes', 'min_contact', 'n_modelled_residues',
                    'involves_modelled', 'inter_chain'):
            self.assertIn(key, geometry)
        # 2MWS has no gaps, so nothing is rebuilt and nothing should clash
        self.assertEqual(geometry['n_modelled_residues'], 0)
        self.assertEqual(geometry['n_clashes'], 0)
        self.assertGreater(geometry['min_contact'], 2.0)

    def test_geometry_check_finds_a_planted_clash(self):
        '''
        The check has to notice an atom driven into a neighbouring chain, which is the
        failure a rebuilt loop would produce: each chain is modelled on its own and knows
        nothing about the chains packed against it.
        '''
        src = os.path.join(self.outdir, 'conformations', '2MWS-alt-1.pdb')
        clean, _ = check_geometry(src)
        self.assertEqual(clean['n_clashes'], 0)

        # move one atom of chain B onto an atom of chain A
        target = None
        for line in open(src):
            if line.startswith('ATOM') and line[21] == 'A' and line[12:16].strip() == 'CA':
                target = line[30:54]
                break
        planted = os.path.join(self.outdir, 'planted.pdb')
        done = False
        with open(src) as fin, open(planted, 'w') as fout:
            for line in fin:
                if (not done and line.startswith('ATOM') and line[21] == 'B'
                        and line[12:16].strip() == 'CA'):
                    line = line[:30] + target + line[54:]
                    done = True
                fout.write(line)

        summary, offending = check_geometry(planted)
        self.assertGreater(summary['n_clashes'], 0)
        self.assertTrue(summary['inter_chain'])
        self.assertLess(summary['min_contact'], 0.1)

    def test_analyze_protein_ignores_heteroatoms(self):
        '''
        Gap detection must not count heteroatoms. Waters and ions are numbered in their
        own range past the end of the chain, so counting them reports the whole span
        between the last residue and the first water as missing, and the structure is
        then rejected as having too large a gap.
        '''
        src = os.path.join(self.outdir, 'conformations', '2MWS-alt-1.pdb')
        M = bb.Molecule()
        M.import_pdb(src, include_hetatm=True)
        before = patcher.analyze_protein(M)

        last = int(M.data['resid'].max())
        with_water = os.path.join(self.outdir, 'with_water.pdb')
        with open(src) as fin, open(with_water, 'w') as fout:
            for line in fin:
                if not line.startswith('END'):
                    fout.write(line)
            # waters numbered well past the chain, as a deposited entry numbers them
            for n in range(1, 21):
                fout.write(f'HETATM{9000 + n:5d}  O   HOH A{last + 30 + n:4d}    '
                           f'{0.0:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00           O\n')
            fout.write('END\n')

        N = bb.Molecule()
        N.import_pdb(with_water, include_hetatm=True)
        self.assertGreater(len(N.data), len(M.data))
        self.assertEqual(patcher.analyze_protein(N), before)

    def test_analyze_protein_with_no_polymer(self):
        '''A chain holding only heteroatoms has no sequence, so it has no gaps.'''
        only_water = os.path.join(self.outdir, 'only_water.pdb')
        with open(only_water, 'w') as fout:
            for n in range(1, 11):
                fout.write(f'HETATM{n:5d}  O   HOH A{n:4d}    '
                           f'{0.0:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00           O\n')
            fout.write('END\n')
        M = bb.Molecule()
        M.import_pdb(only_water, include_hetatm=True)
        self.assertEqual(patcher.analyze_protein(M), [0, 0, 0])

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
