import unittest
import sys
import os
import shutil
import filecmp
from ariba import read_store

modules_dir = os.path.dirname(os.path.abspath(read_store.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


class TestReadStore(unittest.TestCase):
    def test_sort_file(self):
        '''test _sort_file'''
        infile = os.path.join(data_dir, 'read_store_test_sort_file.in')
        expected = os.path.join(data_dir, 'read_store_test_sort_file.out')
        tmpfile = 'tmp.read_store_test_sort_file.out'
        read_store.ReadStore._sort_file(infile, tmpfile)
        self.assertTrue(filecmp.cmp(expected, tmpfile, shallow=False))
        os.unlink(tmpfile)

