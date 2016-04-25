import unittest
import copy
import filecmp
import os
from ariba import summary, flag

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


    def test_dict2key(self):
        '''Test _dict2key'''
        d = {
            'ref_name': 'ref',
            'var_type': '.',
            'known_var_change': '.',
            'ref_ctg_change': '.',
            'var_seq_type': '.'
        }

        self.assertEqual(('ref', '', ''), summary.Summary._dict2key(d))

        d['var_type'] = 'p'
        with self.assertRaises(summary.Error):
            summary.Summary._dict2key(d)

        d['known_var_change'] = 'I42L'
        d['var_seq_type'] = 'p'
        self.assertEqual(('ref', 'p', 'I42L'), summary.Summary._dict2key(d))

        d['ref_ctg_change'] = 'P43Q'
        with self.assertRaises(summary.Error):
            summary.Summary._dict2key(d)

        d['known_var_change'] = '.'
        self.assertEqual(('ref', 'p', 'P43Q'), summary.Summary._dict2key(d))


    def test_load_file(self):
        '''Test _load_file'''
        infile = os.path.join(data_dir, 'summary_test_load_file.in.tsv')
        with open(infile) as f:
             lines = [x.rstrip() for x in f]

        dicts = [summary.Summary._line2dict(x) for x in lines[1:]]
        expected = {
            'cluster.n': {
                'noncoding1': [dicts[0], dicts[1]],
                'noncoding2': [dicts[2]]
            },
            'cluster.p': {
                'presence_absence1': [dicts[3]],
                'presence_absence2': [dicts[4]]
            },
            'cluster.v': {
                'variants_only1': [dicts[5]]
            }
        }

        infile = os.path.join(data_dir, 'summary_test_load_file.in.tsv')
        got = summary.Summary._load_file(infile)
        self.assertEqual(expected, got)


    def test_pc_id_of_longest(self):
        '''Test _pc_id_of_longest'''
        d = {
            'seqname': {
                'key1': {'ref_base_assembled': 10, 'pc_ident': 90.0},
                'key2': {'ref_base_assembled': 20, 'pc_ident': 89.0},
                'key3': {'ref_base_assembled': 50, 'pc_ident': 95.1},
                'key4': {'ref_base_assembled': 42, 'pc_ident': 91.0},
            }
        }

        self.assertEqual(95.1, summary.Summary._pc_id_of_longest(d, 'seqname'))


    def test_to_summary_number_for_seq(self):
        '''Test _to_summary_number_for_seq'''
        tests = [
            (0, 0),
            (64, 0),
            (7, 1),
            (259, 1),
            (15, 1),
            (539, 2),
            (27, 3),
        ]

        for test_flag, expected in tests:
            data_dict = {'name': {
                'key1': {'flag': flag.Flag(test_flag), 'ref_base_assembled': 100, 'pc_ident': 99}
            }}

            self.assertEqual(expected, summary.Summary._to_summary_number_for_seq(data_dict, 'name', 90))


    def test_to_summary_number_for_variant(self):
        '''Test _to_summary_number_for_variant'''
        tests = [
            (1, {'known_var': '1', 'has_known_var': '1', 'ref_ctg_change': 'I42L'}),
            (1, {'known_var': '1', 'has_known_var': '1', 'ref_ctg_change': '.'}),
            (0, {'known_var': '1', 'has_known_var': '0', 'ref_ctg_change': 'I42L'}),
            (0, {'known_var': '1', 'has_known_var': '0', 'ref_ctg_change': '.'}),
            (1, {'known_var': '0', 'has_known_var': '0', 'ref_ctg_change': 'I42L'}),
            (0, {'known_var': '0', 'has_known_var': '0', 'ref_ctg_change': '.'}),
        ]

        for expected, data_dict in tests:
            self.assertEqual(expected, summary.Summary._to_summary_number_for_variant(data_dict))


    def test_gather_output_rows(self):
        '''Test _gather_output_rows'''
        infiles = [
            os.path.join(data_dir, 'summary_test_gather_output_rows.in.1.tsv'),
            os.path.join(data_dir, 'summary_test_gather_output_rows.in.2.tsv')
        ]
        got = summary.Summary._gather_output_rows(infiles, 90)
        expected = [
            ['filename', 'noncoding1', 'noncoding1;var.n.A14T', 'noncoding1;var.n.A6G', 'presence_absence1', 'presence_absence1;var.p.A10V', 'variants_only1'],
            [infiles[0], 1, 1, 0, 3, 1, 0],
            [infiles[1], 1, 1, 1, 3, 1, 0],
        ]

        self.assertEqual(expected, got)


    def test_filter_output_rows(self):
        '''Test _filter_output_rows'''
        rows = [
            ['filename', 'gene1', 'gene2', 'gene3'],
            ['file1', 0, 0, 0],
            ['file2', 1, 0, 3],
            ['file3', 2, 0, 4],
        ]

        expected = [
            ['filename', 'gene1', 'gene3'],
            ['file2', 1, 3],
            ['file3', 2, 4],
        ]

        got = summary.Summary._filter_output_rows(rows)
        self.assertEqual(expected, got)


    def test_write_tsv(self):
        '''Test _write_tsv'''
        tmp_out = 'tmp.out.tsv'
        rows = [
            ['filename', 'gene1', 'gene3'],
            ['file2', 1, 3],
            ['file3', 2, 4],
        ]
        summary.Summary._write_tsv(rows, tmp_out)
        expected = os.path.join(data_dir, 'summary_test_write_tsv.out.tsv')
        self.assertTrue(filecmp.cmp(tmp_out, expected, shallow=False))
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

