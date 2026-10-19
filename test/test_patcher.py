import unittest
import sys
import os
from pathlib import Path
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
        self.outdir = 'resdy_test_patcher'
        os.makedirs(self.outdir, exist_ok=True)
        os.makedirs(os.path.join(self.outdir, 'conformations'), exist_ok=True)

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


    def test_restore_template_columns(self):
        """
        Every atom the template provides takes the template's occupancy and B-factor, and
        every atom Modeller built is marked with occupancy 0. Modeller itself puts the
        template B-factor in the occupancy column and a value of its own in the B column.
        """
        def atom(serial, name, resname, resid, occ, b):
            return (f'ATOM  {serial:5d}  {name:<3s} {resname} A{resid:4d}    '
                    f'{1.0:8.3f}{2.0:8.3f}{3.0:8.3f}{occ:6.2f}{b:6.2f}           {name[0]}\n')

        template = os.path.join(self.outdir, 'template.pdb')
        with open(template, 'w') as fh:
            fh.write(atom(1, 'N', 'GLY', 10, 1.00, 11.11) + atom(2, 'CA', 'GLY', 10, 1.00, 12.22)
                     + atom(3, 'N', 'LYS', 11, 0.50, 21.11) + atom(4, 'CA', 'LYS', 11, 0.50, 22.22))
        # Modeller layout: a residue inserted between the two, and an atom (CB of LYS)
        # that the template residue lacks
        model = os.path.join(self.outdir, 'model.pdb')
        with open(model, 'w') as fh:
            fh.write(atom(1, 'N', 'GLY', 1, 11.11, 37.0) + atom(2, 'CA', 'GLY', 1, 12.22, 37.0)
                     + atom(3, 'N', 'ALA', 2, 0.0, 41.0) + atom(4, 'CA', 'ALA', 2, 0.0, 41.0)
                     + atom(5, 'N', 'LYS', 3, 21.11, 52.0) + atom(6, 'CA', 'LYS', 3, 22.22, 52.0)
                     + atom(7, 'CB', 'LYS', 3, 0.0, 52.0))

        n_template, n_built = patcher.restore_template_columns(
            model, template, pairs=[(0, 0), (2, 1)])
        self.assertEqual((n_template, n_built), (4, 3))
        columns = {(int(l[22:26]), l[12:16].strip()): (float(l[54:60]), float(l[60:66]))
                   for l in open(model) if l.startswith('ATOM')}
        self.assertEqual(columns[(1, 'N')], (1.00, 11.11))
        self.assertEqual(columns[(1, 'CA')], (1.00, 12.22))
        self.assertEqual(columns[(3, 'N')], (0.50, 21.11))
        self.assertEqual(columns[(3, 'CA')], (0.50, 22.22))
        for built in [(2, 'N'), (2, 'CA'), (3, 'CB')]:
            self.assertEqual(columns[built], (0.0, 0.0))

    def test_geometry_counts_built_atoms_as_modelled(self):
        """
        An atom with occupancy 0 was built by Modeller: a clash involving it is reported,
        and reported as involving a modelled atom.
        """
        src = os.path.join(self.outdir, 'conformations', '2MWS-alt-1.pdb')
        target = next(line[30:54] for line in open(src)
                      if line.startswith('ATOM') and line[21] == 'A' and line[12:16].strip() == 'CA')
        planted = os.path.join(self.outdir, 'planted_built.pdb')
        done = False
        with open(src) as fin, open(planted, 'w') as fout:
            for line in fin:
                if (not done and line.startswith('ATOM') and line[21] == 'B'
                        and line[12:16].strip() == 'CA'):
                    line = line[:30] + target + f'{0.0:6.2f}{0.0:6.2f}' + line[66:]
                    done = True
                fout.write(line)
        summary, _ = check_geometry(planted)
        self.assertGreater(summary['n_clashes'], 0)
        self.assertTrue(summary['involves_modelled'])

    @unittest.skipUnless(patcher.modeller_available, 'Modeller is not installed')
    def test_curated_gap_keeps_template_columns(self):
        """
        A residue cut out of a chain is rebuilt by Modeller. In the curated file the
        rebuilt atoms carry occupancy 0, and every other atom keeps the occupancy and
        B-factor it had in the structure that was curated.
        """
        src = os.path.join('demo', 'conformations', '1A6M-alt1A.pdb')
        gapped = os.path.join(self.outdir, 'conformations', '1A6M-alt1A.pdb')
        removed = 50
        with open(src) as fin, open(gapped, 'w') as fout:
            for line in fin:
                if line.startswith(('ATOM', 'HETATM')) and int(line[22:26]) == removed:
                    continue
                fout.write(line)
        fasta = os.path.join(self.outdir, 'conformations', '1A6M.fasta')
        shutil.copyfile(os.path.join('demo', 'conformations', '1A6M.fasta'), fasta)

        outname, largest, _ = patcher.curate(pdb=gapped, fasta=fasta, outdir=self.outdir, gap=10)
        self.assertEqual(largest, 1)

        def columns(path):
            return {(int(l[22:26]), l[12:16].strip()): (l[54:60], l[60:66])
                    for l in open(path) if l.startswith('ATOM')}
        before, after = columns(src), columns(outname)
        occupancies = [float(o) for o, _ in after.values()]
        self.assertTrue(all(0.0 <= o <= 1.0 for o in occupancies))

        built = {k for k, (o, _) in after.items() if float(o) == 0.0}
        self.assertEqual({r for r, _ in built if r != removed}, set(),
                         msg=f'atoms outside residue {removed} marked as built: {sorted(built)}')
        self.assertEqual(built, {k for k in before if k[0] == removed})
        kept = [k for k in before if k[0] != removed]
        self.assertTrue(all(after[k] == before[k] for k in kept))


if __name__ == "__main__":
    unittest.main()
