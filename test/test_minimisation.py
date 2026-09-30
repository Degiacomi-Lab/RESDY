import unittest
import sys
import os
import json
import shutil
import tempfile
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(sys.path[0]), "src"))
import resdy as RD
from resdy import protein

STRUCTURE = 'frag'
SOURCE = os.path.join('demo', 'curated', '1A6M-alt1A.pdb')
N_RESIDUES = 25
WITH_TIP3P = ('amber14-all.xml', 'implicit/gbn2.xml', 'amber14/tip3p.xml')


def _write_fragment(path, zinc=False, water=False):
    '''
    The first N_RESIDUES residues of a curated demo structure, capped with an OXT so the
    last residue matches a C-terminal template. Small enough to minimise in seconds.
    Optionally with a zinc ion and a water placed 4 A from the C-terminus.
    '''
    atoms = [l for l in open(SOURCE) if l.startswith('ATOM')]
    resids = sorted({int(l[22:26]) for l in atoms})[:N_RESIDUES]
    frag = [l for l in atoms if int(l[22:26]) in resids]
    last = [l for l in frag if int(l[22:26]) == resids[-1]]
    xyz = lambda l: np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
    c = [l for l in last if l[12:16].strip() == 'C'][0]
    o = [l for l in last if l[12:16].strip() == 'O'][0]
    ca = [l for l in last if l[12:16].strip() == 'CA'][0]
    # trigonal placement, pointing away from both CA and O. Placing OXT opposite O
    # (O-C-OXT at 180 degrees) gives a singular force that stalls the minimiser.
    unit = lambda v: v / np.linalg.norm(v)
    ox = xyz(c) + 1.25 * unit(unit(xyz(c) - xyz(ca)) + unit(xyz(c) - xyz(o)))
    frag.append(o[:12] + ' OXT' + o[16:30] + f'{ox[0]:8.3f}{ox[1]:8.3f}{ox[2]:8.3f}' + o[54:])
    x, y, z = xyz(c)
    het = []
    if zinc:
        het.append(f'HETATM 9001 ZN    ZN A 901    {x + 4:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00          ZN\n')
    if water:
        het.append(f'HETATM 9002  O   HOH A 902    {x:8.3f}{y + 4:8.3f}{z:8.3f}  1.00 20.00           O\n')
    with open(path, 'w') as fh:
        fh.write(''.join(frag) + 'TER\n' + ''.join(het) + 'END\n')


@unittest.skipUnless(protein.openmm_available, 'openmm is not installed')
class Test_Minimisation(unittest.TestCase):
    '''
    Minimisation settings, the per-structure record, heteroatom retention, and reuse of
    curated structures. Works on a small curated fragment directly, so no download or
    patching is involved.
    '''

    def setUp(self):
        self.outdir = tempfile.mkdtemp(prefix='resdy_min_')
        os.makedirs(os.path.join(self.outdir, 'curated'))
        _write_fragment(self._path(STRUCTURE))

    def tearDown(self):
        shutil.rmtree(self.outdir, ignore_errors=True)

    def _path(self, stem):
        return os.path.join(self.outdir, 'curated', f'{stem}.pdb')

    def _pdb(self, **kwargs):
        kwargs.setdefault('minimise_strucs', 'ALL')
        kwargs.setdefault('minimisation_platform', 'CPU')
        kwargs.setdefault('minimise_max_iterations', 200)
        return RD.PDB(outdir=self.outdir, **kwargs)

    def _record(self, stem=STRUCTURE):
        with open(os.path.join(self.outdir, 'curated', f'{stem}.minimisation.json')) as fh:
            return json.load(fh)

    def _relaxed(self, stem=STRUCTURE):
        return os.path.join(self.outdir, 'curated', f'{stem}_relaxed.pdb')

    def _hetatm_resnames(self, stem=STRUCTURE):
        return {l[17:20].strip() for l in open(self._relaxed(stem)) if l.startswith('HETATM')}

    def test_defaults_are_recorded(self):
        defaults = RD.PDB(outdir=self.outdir, minimise_strucs='ALL').minimisation_settings
        self.assertEqual(defaults['forcefield'], list(protein.DEFAULT_FORCEFIELD))
        self.assertEqual(defaults['tolerance_kj_mol_nm'], 10.0)
        self.assertEqual(defaults['max_iterations'], 1000)
        self.assertFalse(defaults['restrain_heavy_atoms'])
        self.assertIsNone(defaults['platform'])
        P = self._pdb()
        P.apply_minimisation(STRUCTURE)
        self.assertTrue(os.path.exists(self._relaxed()))
        record = self._record()
        self.assertEqual(record['status'], 'relaxed')
        self.assertEqual(record['settings'], P.minimisation_settings)
        self.assertEqual(record['platform_used'], 'CPU')
        self.assertLess(record['energy_final_kj_mol'], record['energy_initial_kj_mol'])
        self.assertGreater(record['heavy_atom_rmsd_a'], 0.0)

    def test_restraints_keep_heavy_atoms_closer_to_input(self):
        self._pdb().apply_minimisation(STRUCTURE)
        rmsd_free = self._record()['heavy_atom_rmsd_a']
        self._pdb(restrain_heavy_atoms=True, restraint_k=10.0).apply_minimisation(STRUCTURE)
        record = self._record()
        self.assertTrue(record['settings']['restrain_heavy_atoms'])
        self.assertLess(record['heavy_atom_rmsd_a'], rmsd_free)

    def test_retained_ion_survives_minimisation(self):
        _write_fragment(self._path(STRUCTURE), zinc=True)
        self._pdb(keep_ions=('ZN',), forcefield=WITH_TIP3P).apply_minimisation(STRUCTURE)
        self.assertEqual(self._record()['status'], 'relaxed')
        self.assertIn('ZN', self._hetatm_resnames())

    def test_retained_water_survives_vacuum_minimisation(self):
        _write_fragment(self._path(STRUCTURE), water=True)
        self._pdb(keep_waters=True,
                  forcefield=('amber14-all.xml', 'amber14/tip3p.xml')).apply_minimisation(STRUCTURE)
        self.assertEqual(self._record()['status'], 'relaxed')
        self.assertIn('HOH', self._hetatm_resnames())

    def test_invalid_settings_fail_at_construction(self):
        for kwargs in ({'forcefield': 'no_such_forcefield.xml'},
                       {'minimise_tolerance': -1},
                       {'minimise_max_iterations': -5},
                       {'restraint_k': -1.0},
                       {'minimisation_platform': 'NoSuchPlatform'}):
            with self.subTest(**kwargs):
                with self.assertRaises(ValueError):
                    self._pdb(**kwargs)

    def test_kept_waters_skip_implicit_minimisation_and_remove_a_stale_relaxed_file(self):
        with open(self._relaxed(), 'w') as fh:
            fh.write('stale\n')
        P = self._pdb(keep_waters=True)
        P.apply_minimisation(STRUCTURE)
        self.assertFalse(os.path.exists(self._relaxed()))
        record = self._record()
        self.assertEqual(record['status'], 'skipped')
        self.assertIn('water', record['reason'])
        self.assertIn(STRUCTURE, P.minimisation_skipped)

    def test_failure_is_recorded(self):
        with open(self._path('broken'), 'w') as fh:
            fh.write('not a structure\n')
        P = self._pdb()
        P.apply_minimisation('broken')
        self.assertFalse(os.path.exists(self._relaxed('broken')))
        self.assertEqual(self._record('broken')['status'], 'failed')
        self.assertIn('broken', P.minimisation_skipped)

    def test_reuse_minimises_again_only_when_needed(self):
        P = self._pdb()
        P.apply_minimisation(STRUCTURE)
        stamp = self._record()['date'], os.path.getmtime(self._relaxed())

        # same settings: the existing relaxed file is kept as it is
        self._pdb()._ensure_minimised(STRUCTURE)
        self.assertEqual((self._record()['date'], os.path.getmtime(self._relaxed())), stamp)

        # different settings: minimised again, and the record says with what
        Q = self._pdb(minimise_tolerance=5.0)
        Q._ensure_minimised(STRUCTURE)
        self.assertEqual(self._record()['settings']['tolerance_kj_mol_nm'], 5.0)
        self.assertNotEqual(os.path.getmtime(self._relaxed()), stamp[1])

        # same settings but the relaxed file has gone: minimised again
        os.remove(self._relaxed())
        Q._ensure_minimised(STRUCTURE)
        self.assertTrue(os.path.exists(self._relaxed()))

    def test_per_call_override_is_recorded_and_not_mistaken_for_current(self):
        P = self._pdb()
        P.apply_minimisation(STRUCTURE, max_iterations=5)
        self.assertEqual(self._record()['settings']['max_iterations'], 5)
        P._ensure_minimised(STRUCTURE)
        self.assertEqual(self._record()['settings'], P.minimisation_settings)

    def test_reuse_minimises_a_structure_with_no_record(self):
        self._pdb()._ensure_minimised(STRUCTURE)
        self.assertTrue(os.path.exists(self._relaxed()))
        self.assertEqual(self._record()['status'], 'relaxed')

    def test_log_collects_records_and_refills_skipped(self):
        self._pdb(keep_waters=True).apply_minimisation(STRUCTURE)
        # a fresh instance stands in for the parent process of a parallel run, which
        # never sees the workers' minimisation_skipped
        P = self._pdb()
        self.assertEqual(P.minimisation_skipped, {})
        df = P.collect_minimisation_log()
        self.assertEqual(list(df['structure']), [STRUCTURE])
        self.assertEqual(df['status'].iloc[0], 'skipped')
        self.assertIn(STRUCTURE, P.minimisation_skipped)
        self.assertTrue(os.path.exists(os.path.join(self.outdir, 'minimisation_log.csv')))


if __name__ == '__main__':
    unittest.main()
