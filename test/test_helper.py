import unittest
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(sys.path[0]), "src"))
from resdy import helper

class Test_Helper(unittest.TestCase):
    def test_download_tool(self):
        tool = helper.get_download_tool()

    def test_require_biobox(self):
        # biobox before 1.1.5 swapped the occupancy and B-factor columns when writing
        module = lambda version: type('biobox', (), {'__version__': version})
        for old in ('1.1.4', '1.0.2', '0.9'):
            with self.assertRaises(ImportError):
                helper.require_biobox(module(old))
        for new in ('1.1.5', '1.2.0', '2.0'):
            helper.require_biobox(module(new))
        # a mocked module, as in the documentation build, has no string version
        helper.require_biobox(type('mock', (), {'__version__': object()}))

if __name__ == "__main__":
    unittest.main()
