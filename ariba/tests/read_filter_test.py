import unittest
import sys
import os
import filecmp
from ariba import read_filter, read_store, external_progs

modules_dir = os.path.dirname(os.path.abspath(read_filter.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


class TestReadFilter(unittest.TestCase):
    def setUp(self):
        self.external_progs = external_progs.ExternalProgs()


    # skip this, as no longer using cdhit2d, but leave it here in case we want
    # to put it back in at a later date
    def _test_run_cdhit_est_2d(self):
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


    # skip this, as no longer using cdhit2d, but leave it here in case we want
    # to put it back in at a later date
    def _test_run(self):
        '''test run'''
        rstore_infile = os.path.join(data_dir, 'read_filter_test_run.in.read_store')
        ref_fasta = os.path.join(data_dir, 'read_filter_test_run.in.ref.fa')
        expected_reads1 = os.path.join(data_dir, 'read_filter_test_run.expected.reads_1.fq')
        expected_reads2 = os.path.join(data_dir, 'read_filter_test_run.expected.reads_2.fq')
        tmp_rstore_prefix = 'tmp.filter_test_run.read_store'
        tmp_reads1 = 'tmp.filter_test_run.reads_1.fq'
        tmp_reads2 = 'tmp.filter_test_run.reads_2.fq'
        rstore = read_store.ReadStore(rstore_infile, tmp_rstore_prefix)
        rfilter = read_filter.ReadFilter(rstore, ref_fasta, '1', sys.stdout)
        got_reads, got_bases = rfilter.run(tmp_reads1, tmp_reads2)
        self.assertEqual(12, got_reads)
        self.assertEqual(912, got_bases)
        self.assertTrue(filecmp.cmp(expected_reads1, tmp_reads1, shallow=False))
        self.assertTrue(filecmp.cmp(expected_reads2, tmp_reads2, shallow=False))
        os.unlink(tmp_reads1)
        os.unlink(tmp_reads2)
        rstore.clean()
