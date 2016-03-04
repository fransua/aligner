
import unittest
import sys
sys.path.append('/home/francisco/Box/aligner/align')

from aligner import Aligner

class TestAligner(unittest.TestCase):
    def test_check_tool(self):
        for tool in ['muscle', 'mafft', 'probcons', 'dialign', 'prout']:
            print tool
            ali = Aligner(tool)
            print ali.exe
            print ali.run
            print '-' * 80





            
if __name__ == '__main__':
    unittest.main()
