import unittest
import sys
import os
import glob

sys.path.insert(0, os.path.join(os.path.dirname(sys.path[0]), "src"))
from resdy.residues import (AA_PRESETS, AA_PROPERTY_KEYS, residue_key, resolve_residue)
from resdy.features.DAS import DAS
from resdy.features.PKAANI import PKAANI
from resdy.features.PROPKA import PROPKA

STANDARD = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
            'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']


def _residues_in_demo():
    '''
    Read the relaxed demo structures once and collect the atom names seen for each residue.

    :returns: Mapping of three-letter code to the set of atom names observed for it.
    :rtype: dict
    '''
    seen = {}
    for path in glob.glob(os.path.join('demo', 'curated', '*_relaxed.pdb')):
        with open(path) as fin:
            for line in fin:
                if line[:4] != 'ATOM':
                    continue
                seen.setdefault(line[17:20].strip(), set()).add(line[12:16].strip())
    return seen


class Test_Presets(unittest.TestCase):

    def test_all_twenty_present(self):
        print('-> Testing that all twenty standard residues have a preset')
        self.assertEqual(sorted(AA_PRESETS), STANDARD)

    def test_preset_shape(self):
        # every preset carries exactly the four keys, as lists of strings, and names its own
        # canonical code first so that the pKa parsers find it
        print('-> Testing preset shape')
        for code, preset in AA_PRESETS.items():
            self.assertEqual(set(preset), set(AA_PROPERTY_KEYS), code)
            for key in AA_PROPERTY_KEYS:
                self.assertIsInstance(preset[key], list, f'{code}: {key}')
                self.assertTrue(all(isinstance(v, str) for v in preset[key]), f'{code}: {key}')
            self.assertEqual(preset['non_modified_codes'][0], code)
            self.assertTrue(preset['atom_select_names_nonmod'], code)
            self.assertEqual(residue_key(preset), code)

    def test_no_code_serves_two_presets(self):
        # a residue code should resolve to one preset only, otherwise a structure would be
        # measured twice under two different anchors
        print('-> Testing that residue codes are not shared between presets')
        owners = {}
        for code, preset in AA_PRESETS.items():
            for res in preset['non_modified_codes'] + preset['modified_codes']:
                self.assertNotIn(res, owners, f'{res} claimed by {owners.get(res)} and {code}')
                owners[res] = code

    def test_modified_anchor_covers_parent(self):
        # the modified anchor list has to keep the parent atom name, since the modified branch
        # selects over the unmodified codes as well
        print('-> Testing that modified anchors keep the parent atom name')
        for code, preset in AA_PRESETS.items():
            if not preset['modified_codes']:
                continue
            for atom in preset['atom_select_names_nonmod']:
                self.assertIn(atom, preset['atom_select_names_modified'], code)


class Test_Feature_Declarations(unittest.TestCase):
    '''
    The per-residue facts a feature needs belong to the feature, not to residues.py. These
    check the contract a dropped-in feature is expected to follow.
    '''

    def test_residues_module_names_no_feature(self):
        # a new feature with per-residue settings must not need an edit to residues.py
        print('-> Testing that residues.py holds no feature-specific tables')
        import resdy.residues as residues
        for gone in ('PKA_SUPPORT', 'DAS_RADII', 'DAS_UNSUPPORTED', 'DAS_WEAK', 'das_note'):
            self.assertFalse(hasattr(residues, gone), gone)

    def test_declared_residues_are_known(self):
        print('-> Testing that declared residues are real presets')
        for cls in (DAS, PROPKA, PKAANI):
            for code in (getattr(cls, 'SUPPORTED_RESIDUES', None) or ()):
                self.assertIn(code, AA_PRESETS, f'{cls.__name__}: {code}')
            for code in (getattr(cls, 'RESIDUE_KWARGS', None) or {}):
                self.assertIn(code, AA_PRESETS, f'{cls.__name__}: {code}')
            reason = getattr(cls, 'UNSUPPORTED_REASON', None)
            if isinstance(reason, dict):
                for code in reason:
                    self.assertIn(code, AA_PRESETS, f'{cls.__name__}: {code}')

    def test_pka_support(self):
        print('-> Testing the pKa backends')
        # pKa-ANI ships models for a subset of what PROPKA3 reports
        self.assertTrue(PKAANI.SUPPORTED_RESIDUES.issubset(PROPKA.SUPPORTED_RESIDUES))

    def test_das_radii(self):
        print('-> Testing the DAS radii')
        for code, kwargs in DAS.RESIDUE_KWARGS.items():
            radii = kwargs['radii']
            # biobox takes the largest gap between consecutive shells and raises on a list of
            # fewer than two
            self.assertGreaterEqual(len(radii), 2, code)
            self.assertEqual(radii, sorted(radii, reverse=True), code)
            # every shell has to clear biobox's 2.0 A clash threshold to return anything
            self.assertGreaterEqual(radii[-1], 2.0, code)
        self.assertEqual(DAS.RESIDUE_KWARGS['LYS']['radii'], [6.3, 5.9, 5.4, 4.8])
        self.assertEqual(DAS.SUPPORTED_RESIDUES, set(DAS.RESIDUE_KWARGS))

    def test_every_residue_is_either_supported_or_explained(self):
        print('-> Testing that DAS is either available or explained for every residue')
        for code in AA_PRESETS:
            supported = code in DAS.SUPPORTED_RESIDUES
            explained = code in DAS.UNSUPPORTED_REASON
            self.assertTrue(supported != explained, code)

    def test_glycine_has_no_shells(self):
        # biobox raises "For flexible mode, a side chain atom must be provided!" on a CA anchor
        print('-> Testing that glycine is excluded from DAS')
        self.assertNotIn('GLY', DAS.SUPPORTED_RESIDUES)

    def test_kwargs_reach_the_constructor(self):
        # what the class declares has to be accepted by its own __init__
        print('-> Testing that declared kwargs are real constructor arguments')
        das = DAS(aa_properties=AA_PRESETS['TRP'], error_filename='no_record',
                  **DAS.RESIDUE_KWARGS['TRP'])
        self.assertEqual(das.radii, [4.7, 4.4, 4.0, 3.6])


class Test_Resolver(unittest.TestCase):

    def test_preset_by_name(self):
        print('-> Testing preset lookup by name')
        self.assertEqual(resolve_residue('tyr'), resolve_residue('TYR'))
        self.assertEqual(resolve_residue(' LYS ')['atom_select_names_nonmod'], ['NZ'])

    def test_result_is_a_copy(self):
        print('-> Testing that the resolver hands back a copy')
        got = resolve_residue('LYS')
        got['non_modified_codes'].append('NOPE')
        self.assertNotIn('NOPE', AA_PRESETS['LYS']['non_modified_codes'])

    def test_unknown_name(self):
        print('-> Testing that an unknown preset name raises')
        with self.assertRaises(KeyError):
            resolve_residue('XYZ')

    def test_dict_key_order_is_free(self):
        # the old implementation compared list(dict) against a fixed order and rejected
        # anything else
        print('-> Testing that dictionary key order does not matter')
        spec = {'atom_select_names_modified': ['NZ'],
                'non_modified_codes': ['LYS'],
                'atom_select_names_nonmod': ['NZ'],
                'modified_codes': ['KCX']}
        self.assertEqual(resolve_residue(spec)['non_modified_codes'], ['LYS'])

    def test_dict_normalisation(self):
        print('-> Testing that empty fields and bare strings are normalised')
        spec = {'non_modified_codes': 'LYS', 'modified_codes': '',
                'atom_select_names_nonmod': ['NZ'], 'atom_select_names_modified': None}
        got = resolve_residue(spec)
        self.assertEqual(got['non_modified_codes'], ['LYS'])
        self.assertEqual(got['modified_codes'], [])
        self.assertEqual(got['atom_select_names_modified'], [])

    def test_dict_missing_key(self):
        print('-> Testing that a missing key raises')
        with self.assertRaises(KeyError):
            resolve_residue({'non_modified_codes': ['LYS'], 'modified_codes': []})

    def test_dict_empty_required_fields(self):
        print('-> Testing that the required fields cannot be empty')
        for empty in ('non_modified_codes', 'atom_select_names_nonmod'):
            spec = {'non_modified_codes': ['LYS'], 'modified_codes': [],
                    'atom_select_names_nonmod': ['NZ'], 'atom_select_names_modified': []}
            spec[empty] = []
            with self.assertRaises(ValueError):
                resolve_residue(spec)

    def test_wrong_type(self):
        print('-> Testing that a non string, non dictionary argument raises')
        with self.assertRaises(ValueError):
            resolve_residue(7)


class Test_Against_Demo_Structures(unittest.TestCase):
    '''
    Offline check of the anchor atom names against real coordinates, so that a typo in the
    table is caught without reaching the CCD over the network.
    '''

    @classmethod
    def setUpClass(cls):
        cls.seen = _residues_in_demo()

    def test_demo_structures_were_found(self):
        print('-> Testing that the demo structures are readable')
        self.assertTrue(self.seen, 'run the tests from the repository root')

    def test_anchor_atoms_exist(self):
        print('-> Testing that every unmodified anchor atom exists in the demo structures')
        for code, preset in AA_PRESETS.items():
            if code not in self.seen:
                continue
            for atom in preset['atom_select_names_nonmod']:
                self.assertIn(atom, self.seen[code],
                              f'{code}: anchor {atom} never seen in demo/curated')


if __name__ == '__main__':
    unittest.main()
