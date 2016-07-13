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


    def test_get_ids(self):
        '''test _get_ids'''
        infile = os.path.join(data_dir, 'cdhit_test_get_ids.fa')
        expected = {'id1', 'id2', 'id3'}
        r = cdhit.Runner(infile, 'out')
        got = r._get_ids(infile)
        self.assertEqual(expected, got)


    def test_rename_clusters(self):
        '''test _rename_clusters'''
        infile = os.path.join(data_dir, 'cdhit_test_rename_clusters.in.fa')
        tmpfile = 'tmp.test_rename_clusters.out.fa'
        expected_file = os.path.join(data_dir, 'cdhit_test_rename_clusters.expected.fa')

        clusters_in = {
            'seq.foo': {'seq.foo', 'seq'},
            'seq.bar': {'seq.bar', 'seq3.spam'},
            'seq4.eggs': {'seq4.eggs'}
        }
        tmp_out = 'tmp.test_rename_clusters.out.fa'
        expected_clusters = {
            'seq.x': {'seq.foo', 'seq'},
            'seq.x.2': {'seq.bar', 'seq3.spam'},
            'seq4.x': {'seq4.eggs'}
        }
        got = cdhit.Runner._rename_clusters(clusters_in, infile, tmpfile)
        self.assertEqual(expected_clusters, got)
        self.assertTrue(filecmp.cmp(expected_file, tmpfile, shallow=False))
        os.unlink(tmpfile)


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
        tmpfile = 'tmp.cdhit_test_run.out.fa'
        r = cdhit.Runner(infile, tmpfile)
        clusters = r.run()
        expected_clusters = {
            '0': {'seq1', 'seq2', 'seq3'},
            '1': {'seq4'},
        }
        self.assertEqual(clusters, expected_clusters)


    def test_run_min_cluster_number_42(self):
        '''test run with min_cluster_number 42'''
        infile = os.path.join(data_dir, 'cdhit_test_run.in.fa')
        expected_outfile = os.path.join(data_dir, 'cdhit_test_run.out.fa')
        tmpfile = 'tmp.cdhit_test_run.out.fa'
        r = cdhit.Runner(infile, tmpfile, min_cluster_number=42)
        clusters = r.run()
        expected_clusters = {
            '42': {'seq1', 'seq2', 'seq3'},
            '43': {'seq4'},
        }
        self.assertEqual(clusters, expected_clusters)


    def test_fake_run(self):
        '''test fake_run'''
        infile = os.path.join(data_dir, 'cdhit_test_fake_run.in.fa')
        expected_outfile = os.path.join(data_dir, 'cdhit_test_fake_run.out.fa')
        tmpfile = 'tmp.cdhit_test_fake_run.out.fa'
        r = cdhit.Runner(infile, tmpfile)
        clusters = r.fake_run()
        expected_clusters = {
            'seq1.x': {'seq1'},
            'seq2.x': {'seq2'},
            'seq3.x': {'seq3'},
            'seq4.x': {'seq4'},
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


    def test_load_user_clusters_file_good_file(self):
        '''test _load_user_clusters_file with good input file'''
        infile = os.path.join(data_dir, 'cdhit_test_load_user_clusters_file.good')
        expected  = {
            'seq1': 'seq1',
            'seq2': 'seq1',
            'seq3': 'seq1',
            'seq4': 'seq4',
            'seq5': 'seq5',
            'seq6': 'seq5',
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
        expected_outfile = os.path.join(data_dir, 'cdhit_test_run_get_clusters_from_dict.out.fa')
        tmpfile = 'tmp.cdhit_test_run_get_clusters_from_dict.out.fa'
        r = cdhit.Runner(fa_infile, tmpfile)
        clusters = r.run_get_clusters_from_file(clusters_infile)
        expected_clusters = {
            'seq1.x': {'seq1', 'seq2'},
            'seq3.x': {'seq3'},
        }
        self.assertEqual(clusters, expected_clusters)
        self.assertTrue(filecmp.cmp(tmpfile, expected_outfile, shallow=False))
        os.unlink(tmpfile)
