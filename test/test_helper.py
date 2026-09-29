import unittest
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(sys.path[0]), "src"))
from resdy import helper

class Test_Helper(unittest.TestCase):
    def test_download_tool(self):
        tool = helper.get_download_tool()

if __name__ == "__main__":
    unittest.main()
