import unittest
import sys
import os
import json
import shutil
import tempfile
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(sys.path[0]), "src"))
import biobox as bb
from resdy import assembly

# a 180 degree rotation about z, and a translation along x
ROT_Z = [[-1, 0, 0, 0.0], [0, -1, 0, 0.0], [0, 0, 1, 0.0]]
SHIFT_X = [[1, 0, 0, 30.0], [0, 1, 0, 0.0], [0, 0, 1, 0.0]]
IDENTITY = [[1, 0, 0, 0.0], [0, 1, 0, 0.0], [0, 0, 1, 0.0]]

ATOMS = [('A', 1, 'LYS', 'N', 'N', (5.0, 1.0, 0.0)),
         ('A', 1, 'LYS', 'CA', 'C', (6.0, 1.5, 0.5)),
         ('B', 1, 'GLY', 'CA', 'C', (12.0, -3.0, 1.0))]


def _biomolecule(number, chains, matrices, author=False, software=False):
    lines = [f'REMARK 350 BIOMOLECULE: {number}']
    if author:
        lines.append('REMARK 350 AUTHOR DETERMINED BIOLOGICAL UNIT: DIMERIC')
    if software:
        lines.append('REMARK 350 SOFTWARE DETERMINED QUATERNARY STRUCTURE: DIMERIC')
    lines.append(f'REMARK 350 APPLY THE FOLLOWING TO CHAINS: {", ".join(chains)}')
    for k, m in enumerate(matrices, start=1):
        for row in range(3):
            r = m[row]
            lines.append(f'REMARK 350   BIOMT{row + 1} {k:3d} {r[0]:9.6f} {r[1]:9.6f} '
                         f'{r[2]:9.6f} {r[3]:14.5f}')
    return lines


def _write_pdb(path, remarks):
    lines = list(remarks)
    for serial, (chain, resid, resname, name, element, (x, y, z)) in enumerate(ATOMS, start=1):
        lines.append(f'ATOM  {serial:5d}  {name:<3s} {resname:3s} {chain}{resid:4d}    '
                     f'{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2s}')
    lines.append('END')
    with open(path, 'w') as fout:
        fout.write('\n'.join(lines) + '\n')


def _read(path):
    with open(path) as fin:
        return fin.read()


def _coordinates(path):
    M = bb.Molecule()
    M.import_pdb(path, include_hetatm=True)
    return M.data['chain'].values, M.points


class Test_Assembly(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.raw = os.path.join(self.tmp, 'raw.pdb')
        self.curated = os.path.join(self.tmp, '1ABC.pdb')

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def _build(self, remarks, **kwargs):
        _write_pdb(self.raw, remarks)
        _write_pdb(self.curated, [])
        biomatrix, determined = assembly.read_biomolecules(self.raw)
        return assembly.build_assembly(self.curated, biomatrix, determined, **kwargs)

    def test_author_determined_is_chosen_over_an_earlier_software_one(self):
        _write_pdb(self.raw, _biomolecule(1, ['A'], [IDENTITY, SHIFT_X], software=True)
                   + _biomolecule(2, ['A', 'B'], [IDENTITY, ROT_Z], author=True, software=True))
        biomatrix, determined = assembly.read_biomolecules(self.raw)
        self.assertEqual(determined, {1: 'software', 2: 'author'})
        self.assertEqual(assembly.choose_biomolecule(biomatrix, determined), (2, 'author'))

    def test_first_listed_without_an_author_determined_one(self):
        _write_pdb(self.raw, _biomolecule(1, ['A'], [IDENTITY, SHIFT_X], software=True)
                   + _biomolecule(2, ['A', 'B'], [IDENTITY, ROT_Z]))
        biomatrix, determined = assembly.read_biomolecules(self.raw)
        self.assertEqual(determined, {1: 'software', 2: 'unspecified'})
        self.assertEqual(assembly.choose_biomolecule(biomatrix, determined), (1, 'software'))

    def test_copies_are_transformed_and_renamed(self):
        record = self._build(_biomolecule(1, ['A', 'B'], [IDENTITY, ROT_Z], author=True))
        self.assertEqual(record['status'], 'built')
        self.assertEqual(record['n_operators'], 2)
        chain, xyz = _coordinates(self.curated)
        self.assertEqual(len(xyz), 6)
        original = np.array([a[5] for a in ATOMS])
        np.testing.assert_allclose(xyz[:3], original, atol=1e-3)
        np.testing.assert_allclose(xyz[3:], original * [-1, -1, 1], atol=1e-3)
        self.assertEqual(list(chain[:3]), ['A', 'A', 'B'])
        self.assertTrue(set(chain[3:]).isdisjoint({'A', 'B'}))
        self.assertEqual(record['measured_chains'], ['A', 'B'])

    def test_measured_chains_are_read_back_for_the_relaxed_copy_too(self):
        self._build(_biomolecule(1, ['A', 'B'], [IDENTITY, ROT_Z], author=True))
        relaxed = os.path.join(self.tmp, '1ABC_relaxed.pdb')
        self.assertEqual(assembly.measured_chains(self.curated), {'A', 'B'})
        self.assertEqual(assembly.measured_chains(relaxed), {'A', 'B'})
        self.assertIsNone(assembly.measured_chains(os.path.join(self.tmp, 'other.pdb')))

    def test_chain_left_out_of_the_biomolecule_is_dropped(self):
        record = self._build(_biomolecule(1, ['A'], [IDENTITY, SHIFT_X], author=True))
        chain, xyz = _coordinates(self.curated)
        self.assertNotIn('B', set(chain))
        self.assertEqual(len(xyz), 4)
        self.assertEqual(record['measured_chains'], ['A'])

    def test_identity_over_every_chain_leaves_the_file_untouched(self):
        _write_pdb(self.raw, _biomolecule(1, ['A', 'B'], [IDENTITY], author=True))
        _write_pdb(self.curated, [])
        before = _read(self.curated)
        biomatrix, determined = assembly.read_biomolecules(self.raw)
        record = assembly.build_assembly(self.curated, biomatrix, determined)
        self.assertEqual(record['status'], 'identity')
        self.assertEqual(_read(self.curated), before)

    def test_no_remark_350(self):
        record = self._build([])
        self.assertEqual(record['status'], 'none')
        self.assertEqual(assembly.measured_chains(self.curated), {'A', 'B'})

    def test_max_atoms_keeps_the_asymmetric_unit(self):
        record = self._build(_biomolecule(1, ['A', 'B'], [IDENTITY, ROT_Z], author=True),
                             max_atoms=5)
        self.assertEqual(record['status'], 'kept')
        self.assertEqual(len(_coordinates(self.curated)[1]), 3)
        self.assertEqual(record['measured_chains'], ['A', 'B'])

    def test_two_character_chain_names_keep_the_asymmetric_unit(self):
        # 70 copies of chain A need more names than single characters provide
        shifts = [[[1, 0, 0, 30.0 * k], [0, 1, 0, 0.0], [0, 0, 1, 0.0]] for k in range(70)]
        record = self._build(_biomolecule(1, ['A'], shifts, author=True))
        self.assertEqual(record['status'], 'kept')
        self.assertIn('single-character', record['reason'])
        self.assertEqual(len(_coordinates(self.curated)[1]), 3)

    def test_record_is_written_beside_the_structure(self):
        self._build(_biomolecule(1, ['A', 'B'], [IDENTITY, ROT_Z], author=True))
        with open(os.path.join(self.tmp, f'1ABC{assembly.RECORD_SUFFIX}')) as fin:
            record = json.load(fin)
        self.assertEqual(record['biomolecule'], 1)
        self.assertEqual(record['determined'], 'author')


if __name__ == '__main__':
    unittest.main()
