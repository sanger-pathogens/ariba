import unittest
import sys
import os
import shutil
import filecmp
import pyfastaq
from ariba import ref_preparer

modules_dir = os.path.dirname(os.path.abspath(ref_preparer.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


class TestRefPreparer(unittest.TestCase):
    def test_get_ref_files(self):
        '''test _get_ref_files using ref_prefix'''
        #_get_ref_files(ref_prefix, presabs, varonly, noncoding, metadata, verbose):
        pass


