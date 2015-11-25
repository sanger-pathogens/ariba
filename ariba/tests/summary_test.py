import unittest
import copy
import filecmp
import os
from ariba import summary, flag

modules_dir = os.path.dirname(os.path.abspath(summary.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')

class TestSummry(unittest.TestCase):
    def test_init(self):
        '''Test init'''
        fofn = os.path.join(data_dir, 'summary_test_init.fofn')
        s = summary.Summary('out', fofn=fofn)
        self.assertEqual(s.filenames, ['file1', 'file2'])
        s = summary.Summary('out', filenames=['file42'])
        self.assertEqual(s.filenames, ['file42'])
        s = summary.Summary('out', fofn=fofn, filenames=['file42'])
        self.assertEqual(s.filenames, ['file42', 'file1', 'file2'])



    def test_line2dict(self):
        '''Test _line2dict'''
        line = '\t'.join(['gene1', '187', '42', '3', '750', '750', '98.93', 'SNP', 'SYN', '.', '66', '66', 'A', 'gene1.scaffold.1', '1047', '67', '67', 'C', '42', 'A', '22,20'])
        s = summary.Summary('out', filenames=['spam', 'eggs'])
        expected = {
            'gene': 'gene1',
            'flag':  flag.Flag(187),
            'reads': 42,
            'cluster': '3',
            'gene_len': 750,
            'assembled': 750,
            'pc_ident': 98.93,
            'var_type': 'SNP',
            'var_effect': 'SYN',
            'new_aa': '.',
            'gene_start': 66,
            'gene_end': 66,
            'gene_nt': 'A',
            'scaffold': 'gene1.scaffold.1',
            'scaff_len': 1047,
            'scaff_start': 67,
            'scaff_end': 67,
            'scaff_nt': 'C',
            'read_depth': 42,
            'alt_bases': 'A',
            'ref_alt_depth': '22,20'
        }
        self.assertEqual(s._line2dict(line), expected)


    def test_load_file(self):
        '''Test _load_file'''
        s = summary.Summary('out', filenames=['spam', 'eggs'])
        infile = os.path.join(data_dir, 'summary_test_load_file.in.tsv')

        lines = [
            ['gene1', '27', '42', '1', '822', '822', '100.0', '.', '.', '.', '.', '.', '.', 'gene1.scaffold.1', '1490', '.', '.', '.', '.', '.', '.'],
            ['gene2', '15', '44', '2', '780', '780', '100.0', '.', '.', '.', '.', '.', '.', 'gene2.scaffold.2', '1124', '.', '.', '.', '.', '.', '.'],
            ['gene2', '15', '46', '2', '780', '770', '99.0', '.', '.', '.', '.', '.', '.', 'gene2.scaffold.3', '1097', '.', '.', '.', '.', '.', '.'],
            ['gene3', '187', '48', '3', '750', '750', '98.93', 'SNP', 'SYN', '.', '318', '318', 'C', 'gene3.scaffold.1', '1047', '319', '319', 'G', '.', '.', '.']
]
        dicts = [s._line2dict('\t'.join(x)) for x in lines]
        expected = {'gene1': [dicts[0]], 'gene2': dicts[1:3], 'gene3': [dicts[3]]}
        got = s._load_file(infile)
        self.assertEqual(expected, got)


    def test_to_summary_number(self):
        '''Test _to_summary_number'''
        s = summary.Summary('out', filenames=['spam', 'eggs'])
        tests = [
            (0, 0),
            (64, 0),
            (7, 1),
            (259, 1),
            (15, 2),
            (539, 3),
            (27, 4),
        ]

        for t in tests:
            l = [{'flag': flag.Flag(t[0]), 'assembled': 42, 'pc_ident': 99}]
            self.assertEqual(s._to_summary_number(l), t[1])

        l = [{'flag': flag.Flag(27), 'assembled': 42, 'pc_ident': 89}]
        self.assertEqual(s._to_summary_number(l), 0)


    def test_gather_output_rows(self):
        '''Test _gather_output_rows'''
        infiles = [
            os.path.join(data_dir, 'summary_test_gather_output_rows.in.1.tsv'),
            os.path.join(data_dir, 'summary_test_gather_output_rows.in.2.tsv')
        ]
        s = summary.Summary('out', filenames=infiles)
        s._gather_output_rows()
        expected = [
            ['filename', 'gene1', 'gene2', 'gene3'],
            [infiles[0], 4, 2, 0],
            [infiles[1], 4, 0, 4],
        ]
        self.assertEqual(expected, s.rows_out)


    def test_filter_output_rows_filter_true(self):
        '''Test _filter_output_rows'''
        s = summary.Summary('out', filenames=['spam', 'eggs'])
        s.rows_out = [
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

        s._filter_output_rows()
        self.assertEqual(s.rows_out, expected)


    def test_filter_output_rows_filter_false(self):
        '''Test _filter_output_rows'''
        s = summary.Summary('out', filenames=['spam', 'eggs'], filter_output=False)
        rows_out = [
            ['filename', 'gene1', 'gene2', 'gene3'],
            ['file1', 0, 0, 0],
            ['file2', 1, 0, 3],
            ['file3', 2, 0, 4],
        ]

        s.rows_out = copy.copy(rows_out)

        s._filter_output_rows()
        self.assertEqual(s.rows_out, rows_out)


    def test_write_tsv(self):
        '''Test _write_tsv'''
        tmp_out = 'tmp.out.tsv'
        s = summary.Summary(tmp_out, filenames=['spam', 'eggs'])
        s.rows_out = [
            ['filename', 'gene1', 'gene3'],
            ['file2', 1, 3],
            ['file3', 2, 4],
        ]
        s._write_tsv()
        expected = os.path.join(data_dir, 'summary_test_write_tsv.out.tsv')
        self.assertTrue(filecmp.cmp(tmp_out, expected, shallow=False))
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
        s = summary.Summary('out', filenames=['spam', 'eggs'])
        s.rows_out = [
            ['filename', 'gene1', 'gene2', 'gene3'],
            ['file1', 0, 1, 0],
            ['file2', 1, 0, 3],
            ['file3', 0, 0, 4],
        ]

        tmp_distances = 'tmp.test.write_distance_matrix.distances'
        s._write_distance_matrix(tmp_distances)
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

