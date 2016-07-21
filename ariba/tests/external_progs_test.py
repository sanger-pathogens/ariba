import unittest
from ariba import external_progs

class TestExternalProgs(unittest.TestCase):
    def test_external_progs_ok(self):
        '''Test that external programs are found'''
        external_progs.ExternalProgs(verbose=True)

