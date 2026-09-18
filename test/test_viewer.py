import unittest
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(sys.path[0]), "src"))
import resdy as RD

class Test_Uniprot(unittest.TestCase):
    def setUp(self):
        # NEEDS FINISHING ONCE VIEWER UPDATED
        self.V = RD.Viewer()
        self.demo_csv = f'demo{os.sep}demo_input.csv'

    def test_organism_funcs(self):
        self.UP.count_organism_proteins(code='UP000007445', reviewed_only=True)
        self.UP.count_organism_proteins(code='UP000007445', reviewed_only=False)


if __name__ == "__main__":
    unittest.main()