import unittest
import os
import copy
import shutil
import filecmp
from ariba import aln_to_metadata

modules_dir = os.path.dirname(os.path.abspath(aln_to_metadata.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


class TestAlnToMetadata(unittest.TestCase):
    def test_load_aln_file(self):
        '''test _load_aln_file'''
        pass


    def test_load_vars_file(self):
        '''test _load_vars_file'''
        pass
