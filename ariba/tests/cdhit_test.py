import unittest
import os
import re
from ariba import cdhit, external_progs


modules_dir = os.path.dirname(os.path.abspath(cdhit.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')
extern_progs = external_progs.ExternalProgs()

class TestCdhit(unittest.TestCase):
    def test_init_fail_infile_missing(self):
        '''test init_fail_infile_missing'''
        with self.assertRaises(cdhit.Error):
            cdhit.Runner('oopsnotafile', 'out')


    def test_init_fail_invalid_memory(self):
        '''test_init_fail_invalid_memory'''
        infile = os.path.join(data_dir, 'cdhit_test_run.in.fa')
        with self.assertRaises(cdhit.Error):
            cdhit.Runner(infile, memory_limit=-10)


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
            r.fake_run()


    def test_load_user_clusters_file_good_file(self):
        '''test _load_user_clusters_file with good input file'''
        infile = os.path.join(data_dir, 'cdhit_test_load_user_clusters_file.good')
        expected  = {
            '0': {'seq1', 'seq2', 'seq3'},
            '1': {'seq4'},
            '2': {'seq5', 'seq6'}
        }

        got = cdhit.Runner._load_user_clusters_file(infile, {'seq' + str(i) for i in range(1,7,1)})
        self.assertEqual(expected, got)

        expected['2'] = {'seq5'}
        got = cdhit.Runner._load_user_clusters_file(infile, {'seq' + str(i) for i in range(1,6,1)})
        self.assertEqual(expected, got)


    def test_load_user_clusters_file_good_file_with_renaming(self):
        '''test _load_user_clusters_file with good input file with some renamed'''
        rename_dict = {'seq2': 'seq2_renamed', 'seq6': 'seq6_renamed'}
        infile = os.path.join(data_dir, 'cdhit_test_load_user_clusters_file.good')
        expected  = {
            '0': {'seq1', 'seq2_renamed', 'seq3'},
            '1': {'seq4'},
            '2': {'seq5', 'seq6_renamed'}
        }

        names = {'seq1', 'seq2_renamed', 'seq3', 'seq4', 'seq5', 'seq6_renamed'}
        got = cdhit.Runner._load_user_clusters_file(infile, names, rename_dict=rename_dict)
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
                cdhit.Runner._load_user_clusters_file(filename, {'seq1', 'seq2', 'seq3'})


    def test_run_get_clusters_from_file(self):
        '''test run_get_clusters_from_file'''
        fa_infile = os.path.join(data_dir, 'cdhit_test_run_get_clusters_from_dict.in.fa')
        clusters_infile = os.path.join(data_dir, 'cdhit_test_run_get_clusters_from_dict.in.clusters')
        r = cdhit.Runner(fa_infile)
        clusters = r.run_get_clusters_from_file(clusters_infile, {'seq1', 'seq2', 'seq3'})
        expected_clusters = {
            '0': {'seq1', 'seq2'},
            '1': {'seq3'},
        }
        self.assertEqual(clusters, expected_clusters)


    def test_run_get_clusters_from_file_with_renaming(self):
        '''test run_get_clusters_from_file with renaming'''
        rename_dict = {'seq2': 'seq2_renamed'}
        fa_infile = os.path.join(data_dir, 'cdhit_test_run_get_clusters_from_dict_rename.in.fa')
        clusters_infile = os.path.join(data_dir, 'cdhit_test_run_get_clusters_from_dict.in.clusters')
        r = cdhit.Runner(fa_infile)
        clusters = r.run_get_clusters_from_file(clusters_infile, {'seq1', 'seq2_renamed', 'seq3'}, rename_dict=rename_dict)
        expected_clusters = {
            '0': {'seq1', 'seq2_renamed'},
            '1': {'seq3'},
        }
        self.assertEqual(clusters, expected_clusters)


    def test_get_run_cmd_with_default_memory(self):
        '''test_get_run_cmd_with_default_memory'''
        fa_infile = os.path.join(data_dir, 'cdhit_test_run_get_clusters_from_dict_rename.in.fa')
        r = cdhit.Runner(fa_infile)
        run_cmd = r.get_run_cmd('foo/bar/file.out')
        match = re.search('^.+ -o foo/bar/file.out -c 0.9 -T 1 -s 0.0 -d 0 -bak 1$', run_cmd)
        self.assertIsNotNone(match, msg="Command output was " + run_cmd)


    def test_get_run_cmd_with_non_default_memory(self):
        '''test_get_run_cmd_with_non_default_memory'''
        fa_infile = os.path.join(data_dir, 'cdhit_test_run_get_clusters_from_dict_rename.in.fa')
        r = cdhit.Runner(fa_infile, memory_limit=900)
        run_cmd = r.get_run_cmd('foo/bar/file.out')
        match = re.search('^.+ -o foo/bar/file.out -c 0.9 -T 1 -s 0.0 -d 0 -bak 1 -M 900$', run_cmd)
        self.assertIsNotNone(match, msg="Command output was " + run_cmd)


    def test_get_run_cmd_with_unlimited_memory(self):
        '''test_get_run_cmd_with_unlimited_memory'''
        fa_infile = os.path.join(data_dir, 'cdhit_test_run_get_clusters_from_dict_rename.in.fa')
        r = cdhit.Runner(fa_infile, memory_limit=0)
        run_cmd = r.get_run_cmd('foo/bar/file.out')
        match = re.search('^.+ -o foo/bar/file.out -c 0.9 -T 1 -s 0.0 -d 0 -bak 1 -M 0$', run_cmd)
        self.assertIsNotNone(match, msg="Command output was " + run_cmd)
