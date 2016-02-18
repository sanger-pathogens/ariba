import unittest
import os
from ariba import external_progs

class TestExternalProgs(unittest.TestCase):
    def test_external_progs_ok(self):
        '''Test that external programs are found'''
        progs = external_progs.ExternalProgs(verbose=True)

