import unittest
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(sys.path[0]), "src"))
import coolpackage as CPN

class Test_Uniprot(unittest.TestCase):
    def setUp(self):
        self.UP = CPN.Uniprot()
        self.demo_csv = f'demo{os.sep}demo_input.csv'

    def test_organism_funcs(self):
        self.UP.count_organism_proteins(code='UP000007445', reviewed_only=True)
        self.UP.count_organism_proteins(code='UP000007445', reviewed_only=False)

        self.UP.get_organism_proteins(code='UP000007445', reviewed_only=True)
        self.UP.get_organism_proteins(code='UP000007445', reviewed_only=False)

        self.UP.from_organism(code='UP000007445', reviewed_only=True)
        self.UP.from_organism(code='UP000007445', reviewed_only=False)

    def test_protein(self):
        self.UP.get_protein_data('P09167')
        self.UP.get_protein_data('P0CG48', '2MWS')

        self.UP.from_csv_file(csv_file=self.demo_csv)

    def test_auxiliary(self):
        self.UP.from_csv_file(csv_file=self.demo_csv)

        self.UP.filter_by_technique(list_of_techniques=['X-ray'])
        self.UP.filter_by_technique(list_of_techniques=['NMR'])
        self.UP.filter_by_technique(list_of_techniques=['Predicted'])

        #self.UP.save_state()

if __name__ == "__main__":
    unittest.main()