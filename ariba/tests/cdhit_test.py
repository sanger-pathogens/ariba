import unittest
import os
import filecmp
from ariba import cdhit

modules_dir = os.path.dirname(os.path.abspath(cdhit.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')

class TestCdhit(unittest.TestCase):
    def test_init_fail_infile_missing(self):
        '''test init_fail_infile_missing'''
        with self.assertRaises(cdhit.Error):
            r = cdhit.Runner('oopsnotafile', 'out')


    def test_get_ids(self):
        '''test _get_ids'''
        infile = os.path.join(data_dir, 'cdhit_test_get_ids.fa')
        expected = {'id1', 'id2', 'id3'}
        r = cdhit.Runner(infile, 'out')
        got = r._get_ids(infile)
        self.assertEqual(expected, got)


    def test_parse_cluster_info_file(self):
        '''test _parse_cluster_info_file'''
        cluster_representatives = {'seq1', 'seq4'}
        infile = os.path.join(data_dir, 'cdhit_test_parse_cluster_info_file.infile')
        got_clusters = cdhit.Runner._parse_cluster_info_file(infile, cluster_representatives)
        expected_clusters = {
            'seq1': {'seq1', 'seq2', 'seq3'},
            'seq4': {'seq4'}
        }
        self.assertEqual(expected_clusters, got_clusters)


    def test_run(self):
        '''test run'''
        infile = os.path.join(data_dir, 'cdhit_test_run.in.fa')
        expected_outfile = os.path.join(data_dir, 'cdhit_test_run.out.fa')
        tmpfile = 'tmp.cdhit_test_run.out.fa'
        r = cdhit.Runner(infile, tmpfile)
        clusters = r.run()
        expected_clusters = {
            'seq1': {'seq1', 'seq2', 'seq3'},
            'seq4': {'seq4'},
        }
        self.assertEqual(clusters, expected_clusters)
        self.assertTrue(filecmp.cmp(tmpfile, expected_outfile, shallow=False))
        os.unlink(tmpfile)


    def test_fake_run(self):
        '''test fake_run'''
        infile = os.path.join(data_dir, 'cdhit_test_fake_run.in.fa')
        expected_outfile = os.path.join(data_dir, 'cdhit_test_fake_run.out.fa')
        tmpfile = 'tmp.cdhit_test_fake_run.out.fa'
        r = cdhit.Runner(infile, tmpfile)
        clusters = r.fake_run()
        expected_clusters = {
            'seq1': {'seq1'},
            'seq2': {'seq2'},
            'seq3': {'seq3'},
            'seq4': {'seq4'},
        }
        self.assertEqual(clusters, expected_clusters)
        self.assertTrue(filecmp.cmp(tmpfile, expected_outfile, shallow=False))
        os.unlink(tmpfile)


    def test_fake_run_fail(self):
        '''test fake_run with non-unique names'''
        infile = os.path.join(data_dir, 'cdhit_test_fake_run.non-unique.in.fa')
        tmpfile = 'tmp.cdhit_test_fake_run.out.non-unique.fa'
        r = cdhit.Runner(infile, tmpfile)
        with self.assertRaises(cdhit.Error):
            clusters = r.fake_run()
        os.unlink(tmpfile)

