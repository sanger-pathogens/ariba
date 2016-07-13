import unittest
import os
import filecmp
from ariba import cdhit, external_progs

modules_dir = os.path.dirname(os.path.abspath(cdhit.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')
extern_progs = external_progs.ExternalProgs()

class TestCdhit(unittest.TestCase):
    def test_init_fail_infile_missing(self):
        '''test init_fail_infile_missing'''
        with self.assertRaises(cdhit.Error):
            r = cdhit.Runner('oopsnotafile', 'out')


    def test_get_clusters_from_bak_file(self):
        '''test _get_clusters_from_bak_file'''
        infile = os.path.join(data_dir, 'cdhit_test_get_clusters_from_bak_file.in')
        expected = {
            '0': {'seq1', 'seq2', 'seq3'},
            '1': {'seq4'},
            '2': {'seq5'}
        }
        got = cdhit.Runner._get_clusters_from_bak_file(infile)
        self.assertEqual(expected, got)

        expected = {
            '42': {'seq1', 'seq2', 'seq3'},
            '43': {'seq4'},
            '44': {'seq5'}
        }
        got = cdhit.Runner._get_clusters_from_bak_file(infile, min_cluster_number=42)
        self.assertEqual(expected, got)


    def test_run(self):
        '''test run'''
        infile = os.path.join(data_dir, 'cdhit_test_run.in.fa')
        expected_outfile = os.path.join(data_dir, 'cdhit_test_run.out.fa')
        r = cdhit.Runner(infile)
        clusters = r.run()
        expected_clusters = {
            '0': {'seq1', 'seq2', 'seq3'},
            '1': {'seq4'},
        }
        self.assertEqual(clusters, expected_clusters)


    def test_run_min_cluster_number_42(self):
        '''test run with min_cluster_number 42'''
        infile = os.path.join(data_dir, 'cdhit_test_run.in.fa')
        r = cdhit.Runner(infile, min_cluster_number=42)
        clusters = r.run()
        expected_clusters = {
            '42': {'seq1', 'seq2', 'seq3'},
            '43': {'seq4'},
        }
        self.assertEqual(clusters, expected_clusters)


    def test_fake_run(self):
        '''test fake_run'''
        infile = os.path.join(data_dir, 'cdhit_test_fake_run.in.fa')
        r = cdhit.Runner(infile)
        clusters = r.fake_run()
        expected_clusters = {
            '0': {'seq1'},
            '1': {'seq2'},
            '2': {'seq3'},
            '3': {'seq4'},
        }
        self.assertEqual(clusters, expected_clusters)


    def test_fake_run_cluster_min_42(self):
        '''test fake_run min_cluster 42'''
        infile = os.path.join(data_dir, 'cdhit_test_fake_run.in.fa')
        r = cdhit.Runner(infile, min_cluster_number=42)
        clusters = r.fake_run()
        expected_clusters = {
            '42': {'seq1'},
            '43': {'seq2'},
            '44': {'seq3'},
            '45': {'seq4'},
        }
        self.assertEqual(clusters, expected_clusters)


    def test_fake_run_fail(self):
        '''test fake_run with non-unique names'''
        infile = os.path.join(data_dir, 'cdhit_test_fake_run.non-unique.in.fa')
        r = cdhit.Runner(infile)
        with self.assertRaises(cdhit.Error):
            clusters = r.fake_run()


    def test_load_user_clusters_file_good_file(self):
        '''test _load_user_clusters_file with good input file'''
        infile = os.path.join(data_dir, 'cdhit_test_load_user_clusters_file.good')
        expected  = {
            '0': {'seq1', 'seq2', 'seq3'},
            '1': {'seq4'},
            '2': {'seq5', 'seq6'}
        }

        got = cdhit.Runner._load_user_clusters_file(infile)
        self.assertEqual(expected, got)


    def test_load_user_clusters_file_bad_file(self):
        '''test _load_user_clusters_file with bad input files'''
        infiles = [
            os.path.join(data_dir, 'cdhit_test_load_user_clusters_file.bad1'),
            os.path.join(data_dir, 'cdhit_test_load_user_clusters_file.bad2'),
            os.path.join(data_dir, 'cdhit_test_load_user_clusters_file.bad3')
        ]
        for filename in infiles:
            with self.assertRaises(cdhit.Error):
                cdhit.Runner._load_user_clusters_file(filename)


    def test_run_get_clusters_from_file(self):
        '''test run_get_clusters_from_file'''
        fa_infile = os.path.join(data_dir, 'cdhit_test_run_get_clusters_from_dict.in.fa')
        clusters_infile = os.path.join(data_dir, 'cdhit_test_run_get_clusters_from_dict.in.clusters')
        r = cdhit.Runner(fa_infile)
        clusters = r.run_get_clusters_from_file(clusters_infile)
        expected_clusters = {
            '0': {'seq1', 'seq2'},
            '1': {'seq3'},
        }
        self.assertEqual(clusters, expected_clusters)
