import unittest
import sys
import os
import shutil
import filecmp
import pyfastaq
from ariba import read_filter, read_store, external_progs

modules_dir = os.path.dirname(os.path.abspath(read_filter.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


class TestReadFilter(unittest.TestCase):
    def setUp(self):
        self.external_progs = external_progs.ExternalProgs()


    def test_run_cdhit_est_2d(self):
        '''test _run_cdhit_est_2d'''
        reads_in = os.path.join(data_dir, 'read_filter_test_run_cdhit_est_2d.reads.in.fa')
        ref_in = os.path.join(data_dir, 'read_filter_test_run_cdhit_est_2d.ref.in.fa')
        tmp_out = 'tmp.test_run_cdhit_est_2d.out'
        read_filter.ReadFilter._run_cdhit_est_2d(ref_in, reads_in, tmp_out, self.external_progs.exe('cdhit2d'))
        got_clustr = tmp_out + '.clstr'
        expected_clstr = os.path.join(data_dir, 'read_filter_test_run_cdhit_est_2d.expected.clstr')
        self.assertTrue(filecmp.cmp(got_clustr, expected_clstr, shallow=False))
        os.unlink(got_clustr)


    def test_cdhit_clstr_to_reads(self):
        '''test _cdhit_clstr_to_reads'''
        infile = os.path.join(data_dir, 'read_filter_cdhit_clstr_to_reads.in.clstr')
        expected = {1, 3, 5, 7, 9, 11}
        got = read_filter.ReadFilter._cdhit_clstr_to_reads(infile)
        self.assertEqual(expected, got)

