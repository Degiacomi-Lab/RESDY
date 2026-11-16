import unittest
import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(sys.path[0]), "src"))
from resdy.geometry import distance_to_other_chains


def _atoms(*rows):
    '''(chain, resname, element, x) rows on the x axis -> the arrays the function takes'''
    chain, resname, element, x = zip(*rows)
    xyz = np.zeros((len(rows), 3))
    xyz[:, 0] = x
    return np.array(chain), np.array(resname), np.array(element), xyz


class Test_DistanceToOtherChains(unittest.TestCase):
    '''Min_Dist_Other_Chain: anchor to the closest protein heavy atom of any other chain'''

    def test_two_chains(self):
        chain, resname, element, xyz = _atoms(('A', 'LYS', 'N', 0.0), ('A', 'LYS', 'C', 1.0),
                                              ('B', 'ALA', 'C', 7.5))
        anchor = np.array([True, False, False])
        d = distance_to_other_chains(chain, resname, element, xyz, anchor)
        np.testing.assert_allclose(d, [7.5])

    def test_closest_of_several_chains(self):
        # one float, whichever chain is the closest
        chain, resname, element, xyz = _atoms(('A', 'LYS', 'N', 0.0), ('B', 'GLY', 'C', 9.0),
                                              ('C', 'GLY', 'C', -4.0), ('D', 'GLY', 'C', 20.0))
        anchor = np.array([True, False, False, False])
        np.testing.assert_allclose(distance_to_other_chains(chain, resname, element, xyz, anchor),
                                   [4.0])

    def test_each_anchor_measures_against_the_chains_other_than_its_own(self):
        chain, resname, element, xyz = _atoms(('A', 'LYS', 'N', 0.0), ('A', 'GLY', 'C', 0.5),
                                              ('B', 'LYS', 'N', 6.0), ('B', 'GLY', 'C', 5.5))
        anchor = np.array([True, False, True, False])
        np.testing.assert_allclose(distance_to_other_chains(chain, resname, element, xyz, anchor),
                                   [5.5, 5.5])

    def test_single_chain_is_nan(self):
        chain, resname, element, xyz = _atoms(('A', 'LYS', 'N', 0.0), ('A', 'GLY', 'C', 3.0))
        d = distance_to_other_chains(chain, resname, element, xyz, np.array([True, False]))
        self.assertTrue(np.isnan(d).all())

    def test_ions_waters_and_ligands_are_not_a_chain(self):
        chain, resname, element, xyz = _atoms(('A', 'LYS', 'N', 0.0), ('B', 'ZN', 'ZN', 2.0),
                                              ('C', 'HOH', 'O', 2.5), ('D', 'HEM', 'FE', 3.0))
        d = distance_to_other_chains(chain, resname, element, xyz,
                                     np.array([True, False, False, False]))
        self.assertTrue(np.isnan(d).all())

        chain, resname, element, xyz = _atoms(('A', 'LYS', 'N', 0.0), ('B', 'ZN', 'ZN', 2.0),
                                              ('B', 'GLY', 'C', 8.0))
        d = distance_to_other_chains(chain, resname, element, xyz,
                                     np.array([True, False, False]))
        np.testing.assert_allclose(d, [8.0])

    def test_hydrogens_are_ignored(self):
        chain, resname, element, xyz = _atoms(('A', 'LYS', 'N', 0.0), ('B', 'GLY', 'H', 1.5),
                                              ('B', 'GLY', 'C', 2.5))
        d = distance_to_other_chains(chain, resname, element, xyz,
                                     np.array([True, False, False]))
        np.testing.assert_allclose(d, [2.5])

    def test_modified_and_alias_codes_count_as_protein(self):
        # residue names from the presets: a protonation alias and a modified lysine
        for partner in ('HID', 'KCX'):
            chain, resname, element, xyz = _atoms(('A', 'LYS', 'N', 0.0), ('B', partner, 'C', 3.0))
            d = distance_to_other_chains(chain, resname, element, xyz, np.array([True, False]))
            np.testing.assert_allclose(d, [3.0])

    def test_no_anchor(self):
        chain, resname, element, xyz = _atoms(('A', 'LYS', 'N', 0.0), ('B', 'GLY', 'C', 3.0))
        d = distance_to_other_chains(chain, resname, element, xyz, np.array([False, False]))
        self.assertEqual(d.shape, (0,))


if __name__ == "__main__":
    unittest.main()
