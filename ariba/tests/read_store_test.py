import unittest
import sys
import os
import shutil
import filecmp
import pyfastaq
from ariba import read_store

modules_dir = os.path.dirname(os.path.abspath(read_store.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


def file_to_list(infile):
    f = pyfastaq.utils.open_file_read(infile)
    lines = [x for x in f.readlines()]
    pyfastaq.utils.close(f)
    return lines


class TestReadStore(unittest.TestCase):
    def test_sort_file(self):
        '''test _sort_file'''
        infile = os.path.join(data_dir, 'read_store_test_sort_file.in')
        expected = os.path.join(data_dir, 'read_store_test_sort_file.out')
        tmpfile = 'tmp.read_store_test_sort_file.out'
        read_store.ReadStore._sort_file(infile, tmpfile)
        self.assertTrue(filecmp.cmp(expected, tmpfile, shallow=False))
        os.unlink(tmpfile)


    def test_compress_and_index_file(self):
        '''Test _compress_and_index_file'''
        infile = os.path.join(data_dir, 'read_store_test_compress_and_index_file.in')
        tmpfile = 'tmp.test_compress_and_index_file.in'
        tmpfile_gz = 'tmp.test_compress_and_index_file.in.gz'
        shutil.copyfile(infile, tmpfile)
        read_store.ReadStore._compress_and_index_file(tmpfile)
        self.assertTrue(os.path.exists(tmpfile_gz))
        expected_lines = file_to_list(infile)
        got_lines = file_to_list(tmpfile_gz)
        self.assertEqual(expected_lines, got_lines)
        self.assertTrue(os.path.exists(tmpfile_gz + '.tbi'))
        os.unlink(tmpfile)
        os.unlink(tmpfile_gz)
        os.unlink(tmpfile_gz + '.tbi')


