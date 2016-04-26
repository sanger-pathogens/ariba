import unittest
import copy
import filecmp
import os
from ariba import flag, summary, summary_cluster, summary_sample

modules_dir = os.path.dirname(os.path.abspath(summary.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')

class TestSummary(unittest.TestCase):
    def test_init(self):
        '''Test init'''
        fofn = os.path.join(data_dir, 'summary_test_init.fofn')
        s = summary.Summary('out', fofn=fofn)
        self.assertEqual(s.filenames, ['file1', 'file2'])
        s = summary.Summary('out', filenames=['file42'])
        self.assertEqual(s.filenames, ['file42'])
        s = summary.Summary('out', fofn=fofn, filenames=['file42'])
        self.assertEqual(s.filenames, ['file42', 'file1', 'file2'])


    def test_load_input_files(self):
        '''Test _load_input_files'''
        file1 = os.path.join(data_dir, 'summary_test_load_input_files.1.tsv')
        file2 = os.path.join(data_dir, 'summary_test_load_input_files.2.tsv')
        sample1 = summary_sample.SummarySample(file1)
        sample2 = summary_sample.SummarySample(file2)
        sample1.run()
        sample2.run()
        got = summary.Summary._load_input_files([file1, file2], 90)
        expected = {file1: sample1, file2: sample2}
        self.assertEqual(expected, got)


    def test_get_all_cluster_names(self):
        '''Test _get_all_cluster_names'''
        file1 = os.path.join(data_dir, 'summary_test_get_all_cluster_names.1.tsv')
        file2 = os.path.join(data_dir, 'summary_test_get_all_cluster_names.2.tsv')
        samples = summary.Summary._load_input_files([file1, file2], 90)
        got = summary.Summary._get_all_cluster_names(samples)
        expected = {'cluster.n.1', 'cluster.v.1', 'cluster.p.1', 'cluster.p.2'}
        self.assertEqual(expected, got)


    def test_get_all_variant_columns(self):
        '''Test _get_all_variant_columns'''
        file1 = os.path.join(data_dir, 'summary_test_get_all_cluster_names.1.tsv')
        file2 = os.path.join(data_dir, 'summary_test_get_all_cluster_names.2.tsv')
        samples = summary.Summary._load_input_files([file1, file2], 90)
        got = summary.Summary._get_all_variant_columns(samples)
        expected = {
            'cluster.p.2': {('presence_absence1', 'A10V')},
            'cluster.n.1': {('noncoding1', 'A6G'), ('noncoding1', 'A14T')},
            'cluster.p.1': {('presence_absence1', 'A10V')},
        }
        self.assertEqual(expected, got)


    def test_gather_output_rows(self):
        '''Test _gather_output_rows'''
        infiles = [
            os.path.join(data_dir, 'summary_test_gather_output_rows.in.1.tsv'),
            os.path.join(data_dir, 'summary_test_gather_output_rows.in.2.tsv')
        ]
        s = summary.Summary('out', filenames=infiles)
        s.samples = summary.Summary._load_input_files(infiles, 90)
        expected = {
            infiles[0]: {
                'noncoding1': {
                    'assembled': 'yes',
                    'ref_seq': 'noncoding1',
                    'any_var': 'yes',
                    'pct_id': '98.33',
                },
                'presence_absence1': {
                    'assembled': 'yes',
                    'ref_seq': 'presence_absence1',
                    'any_var': 'yes',
                    'pct_id': '98.96',
                },
                'variants_only1': {
                    'assembled': 'no',
                    'ref_seq': 'NA',
                    'any_var': 'NA',
                    'pct_id': 'NA',
                }
            },
            infiles[1]: {
                'noncoding1': {
                    'assembled': 'yes',
                    'ref_seq': 'noncoding1',
                    'any_var': 'yes',
                    'pct_id': '98.33',
                },
                'presence_absence1': {
                    'assembled': 'yes',
                    'ref_seq': 'presence_absence1',
                    'pct_id': '98.96',
                    'any_var': 'yes',
                },
                'variants_only1': {
                    'assembled': 'no',
                    'ref_seq': 'NA',
                    'any_var': 'NA',
                    'pct_id': 'NA',
                }
            },
        }
        got = s._gather_output_rows()
        self.assertEqual(expected, got)

        s.include_all_variant_columns = True
        expected[infiles[0]]['noncoding1']['noncoding1.A14T'] = 'yes'
        expected[infiles[0]]['noncoding1']['noncoding1.A6G'] = 'no'
        expected[infiles[0]]['presence_absence1']['presence_absence1.A10V'] = 'yes'
        expected[infiles[1]]['noncoding1']['noncoding1.A14T'] = 'yes'
        expected[infiles[1]]['noncoding1']['noncoding1.A6G'] = 'yes'
        expected[infiles[1]]['presence_absence1']['presence_absence1.A10V'] = 'yes'
        got = s._gather_output_rows()
        self.assertEqual(expected, got)


    def test_write_csv(self):
        '''Test _write_csv'''
        tmp_out = 'tmp.out.tsv'
        rows = {
            'file1': {
                'cluster.n.1': {
                    'assembled': 'yes',
                    'ref_seq': 'noncoding1',
                    'any_var': 'yes',
                    'pct_id': '98.33',
                    'noncoding1.A14T': 'yes'
                },
                'cluster.p.1': {
                    'assembled': 'yes',
                    'ref_seq': 'presence_absence1',
                    'any_var': 'yes',
                    'pct_id': '98.96',
                    'presence_absence1.I42L': 'yes'
                },
                'cluster.v.1': {
                    'assembled': 'yes',
                    'ref_seq': 'varonly1',
                    'any_var': 'no',
                    'pct_id': '99.42',
                }
            },
            'file2': {
                'cluster.n.1': {
                    'assembled': 'yes',
                    'ref_seq': 'noncoding1',
                    'any_var': 'no',
                    'pct_id': '98.33',
                    'noncoding1.A14T': 'no'
                },
                'cluster.p.1': {
                    'assembled': 'yes',
                    'ref_seq': 'presence_absence1',
                    'pct_id': '98.96',
                    'any_var': 'no',
                    'presence_absence1.I42L': 'no'
                },
                'cluster.v.1': {
                    'assembled': 'no',
                    'ref_seq': 'NA',
                    'any_var': 'NA',
                    'pct_id': 'NA',
                }
            },
        }
        filenames = ['file1', 'file2']
        got_lines = summary.Summary._write_csv(filenames, rows, tmp_out, phandango=False)
        expected_file = os.path.join(data_dir, 'summary_test_write_tsv.out.not_phandango.csv')
        self.assertTrue(filecmp.cmp(tmp_out, expected_file, shallow=False))
        os.unlink(tmp_out)
        expected_lines = [
            'name,cluster.n.1,cluster.n.1.ref,cluster.n.1.idty,cluster.n.1.any_var,cluster.n.1.noncoding1.A14T,cluster.p.1,cluster.p.1.ref,cluster.p.1.idty,cluster.p.1.any_var,cluster.p.1.presence_absence1.I42L,cluster.v.1,cluster.v.1.ref,cluster.v.1.idty,cluster.v.1.any_var',
            'file1,yes,noncoding1,98.33,yes,yes,yes,presence_absence1,98.96,yes,yes,yes,varonly1,99.42,no',
            'file2,yes,noncoding1,98.33,no,no,yes,presence_absence1,98.96,no,no,no,NA,NA,NA'
]
        self.assertEqual(expected_lines, got_lines)

        got_lines = summary.Summary._write_csv(filenames, rows, tmp_out, phandango=True)
        expected_lines = [
            'name,cluster.n.1:o1,cluster.n.1.ref:o2,cluster.n.1.idty:c1,cluster.n.1.any_var:o1,cluster.n.1.noncoding1.A14T:o1,cluster.p.1:o1,cluster.p.1.ref:o2,cluster.p.1.idty:c1,cluster.p.1.any_var:o1,cluster.p.1.presence_absence1.I42L:o1,cluster.v.1:o1,cluster.v.1.ref:o2,cluster.v.1.idty:c1,cluster.v.1.any_var:o1',
            'file1,yes,noncoding1,98.33,yes,yes,yes,presence_absence1,98.96,yes,yes,yes,varonly1,99.42,no',
            'file2,yes,noncoding1,98.33,no,no,yes,presence_absence1,98.96,no,no,no,NA,NA,NA'
]
        self.assertEqual(expected_lines, got_lines)
        expected_file = os.path.join(data_dir, 'summary_test_write_tsv.out.phandango.csv')
        self.assertTrue(filecmp.cmp(tmp_out, expected_file, shallow=False))
        os.unlink(tmp_out)


    def test_write_phandango_csv(self):
        '''Test _write_phandango_csv'''
        tmp_out = 'tmp.test_write_phandango.csv'
        rows = [
            ['filename', 'seq1', 'seq1;var.p.I14L', 'seq1;var.p.P42Q', 'seq2', 'seq2;var.n.A14T'],
            ['file1', 3, 0, 1, 3, 1],
            ['file2', 3, 1, 0, 3, 0],
            ['file3', 1, 0, 0, 3, 0],
            ['file4', 2, 1, 0, 0, 0],
        ]
        summary.Summary._write_phandango_csv(rows, tmp_out)
        expected = os.path.join(data_dir, 'summary_test_write_phandango_csv.csv')
        self.assertTrue(filecmp.cmp(expected, tmp_out, shallow=False))
        os.unlink(tmp_out)


    def test_distance_score_bewteen_values(self):
        '''Test _distance_score_bewteen_values'''
        tests = [
            ((0, 0), 0),
            ((0, 1), 1),
            ((0, 2), 1),
            ((1, 0), 1),
            ((1, 1), 0),
            ((1, 2), 0),
            ((2, 0), 1),
            ((2, 1), 0),
            ((2, 2), 0),
        ]

        for (val1, val2), expected in tests:
            self.assertEqual(expected, summary.Summary._distance_score_between_values(val1, val2))

        # check distance calculation is commutative
        for val1 in range(5):
            for val2 in range(5):
                d1 = summary.Summary._distance_score_between_values(val1, val2)
                d2 = summary.Summary._distance_score_between_values(val2, val1)
                self.assertEqual(d1, d2)


    def test_distance_score_between_lists(self):
        '''Test _distance_score_between_lists'''
        list1 = ['na', 0, 0, 0, 1, 1, 1, 2, 2, 2]
        list2 = ['na', 0, 1, 2, 0, 1, 2, 1, 2, 2]
        self.assertEqual(3, summary.Summary._distance_score_between_lists(list1, list2))


    def test_write_distance_matrix(self):
        '''Test _write_distance_matrix'''
        rows = [
            ['filename', 'gene1', 'gene2', 'gene3'],
            ['file1', 0, 1, 0],
            ['file2', 1, 0, 3],
            ['file3', 0, 0, 4],
        ]

        tmp_distances = 'tmp.test.write_distance_matrix.distances'
        summary.Summary._write_distance_matrix(rows, tmp_distances)
        expected = os.path.join(data_dir, 'summary_test_write_distance_matrix.distances')
        self.assertTrue(filecmp.cmp(expected, tmp_distances, shallow=False))
        os.unlink(tmp_distances)


    def test_newick_from_dist_matrix(self):
        '''Test _newick_from_dist_matrix'''
        tmp_tree = 'tmp.test.newick_from_dist_matrix.tre'
        dist_file = os.path.join(data_dir, 'summary_test_newick_from_dist_matrix.distances')
        summary.Summary._newick_from_dist_matrix(dist_file, tmp_tree)
        expected = os.path.join(data_dir, 'summary_test_newick_from_dist_matrix.tre')
        self.assertTrue(filecmp.cmp(expected, tmp_tree, shallow=False))
        os.unlink(tmp_tree)


    def test_write_phandango_files(self):
        '''Test _write_phandango_files'''
        tmp_prefix = 'tmp.test.write_phandango_files'
        rows = [
            ['filename', 'seq1', 'seq1;var.p.I14L', 'seq1;var.p.P42Q', 'seq2', 'seq2;var.n.A14T'],
            ['file1', 3, 0, 1, 3, 1],
            ['file2', 3, 1, 0, 3, 0],
            ['file3', 1, 0, 0, 3, 0],
            ['file4', 2, 1, 0, 0, 0],
        ]
        summary.Summary._write_phandango_files(rows, tmp_prefix)
        expected_csv = os.path.join(data_dir, 'summary_test_write_phandango_files.csv')
        expected_tre = os.path.join(data_dir, 'summary_test_write_phandango_files.tre')
        self.assertTrue(filecmp.cmp(expected_csv, tmp_prefix + '.csv', shallow=False))
        self.assertTrue(filecmp.cmp(expected_tre, tmp_prefix + '.tre', shallow=False))
        os.unlink(tmp_prefix + '.csv')
        os.unlink(tmp_prefix + '.tre')

