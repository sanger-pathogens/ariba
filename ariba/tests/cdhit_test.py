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


    def test_enumerate_fasta(self):
        '''test _enumerate_fasta'''
        infile = os.path.join(data_dir, 'cdhit_test_enumerate_fasta.in.fa')
        expected_outfile = os.path.join(data_dir, 'cdhit_test_enumerate_fasta.out.fa')
        tmpfile = 'tmp.test_enumerate_fasta.out.fa'
        expected_dict = {'1': 'a', '2': 'b', '3': 'c'}
        r = cdhit.Runner(infile, 'out')
        got_dict = r._enumerate_fasta(infile, tmpfile)
        self.assertTrue(filecmp.cmp(expected_outfile, tmpfile, shallow=False))
        self.assertEqual(expected_dict, got_dict)
        os.unlink(tmpfile)


    def test_get_ids(self):
        '''test _get_ids'''
        infile = os.path.join(data_dir, 'cdhit_test_get_ids.fa')
        expected = {'id1', 'id2', 'id3'}
        r = cdhit.Runner(infile, 'out')
        got = r._get_ids(infile)
        self.assertEqual(expected, got)


    def test_parse_cluster_info_file(self):
        '''test _parse_cluster_info_file'''
        infile = os.path.join(data_dir, 'cdhit_test_parse_cluster_info_file.in.fa')
        r = cdhit.Runner(infile, 'out')
        names_dict = {str(i): 'seq' + str(i) for i in range(1,5)}
        cluster_representatives = {'1', '4'}
        cluster_file = os.path.join(data_dir, 'cdhit_test_parse_cluster_info_file.out.fa.bak.clstr')
        got_clusters, got_reps = r._parse_cluster_info_file(cluster_file, names_dict, cluster_representatives)
        expected_clusters = {
            '0': {'seq1', 'seq2', 'seq3'},
            '1': {'seq4'}
        }
        expected_reps = {'1': '0', '4': '1'}
        self.assertEqual(expected_clusters, got_clusters)
        self.assertEqual(expected_reps, got_reps)


    def test_rename_fasta(self):
        '''test _rename_fasta'''
        infile = os.path.join(data_dir, 'cdhit_test_rename_fasta.in.fa')
        tmpfile = 'tmp.rename_fasta.out.fa'
        expected = os.path.join(data_dir, 'cdhit_test_rename_fasta.out.fa')
        names_dict = {'a': 'seq1', 'b': 'seq2', 'c': 'seq3'}
        r = cdhit.Runner(infile, 'out')
        r._rename_fasta(infile, tmpfile, names_dict)
        self.assertTrue(filecmp.cmp(expected, tmpfile, shallow=False))
        os.unlink(tmpfile)


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
            '0': {'seq1'},
            '1': {'seq2'},
            '2': {'seq3'},
            '3': {'seq4'},
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

