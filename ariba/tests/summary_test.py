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



    def test_line2dict(self):
        '''Test _line2dict'''
        line = 'refname\treftype\t19\t78\tcluster\t120\t120\t98.33\tctg_name\t279\t1\tSNP\tn\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1_n_A14T_N_ref has wild type, foo bar\tsome free text'

        expected = {
            'ref_name': 'refname',
            'ref_type': 'reftype',
            'flag': flag.Flag(19),
            'reads': 78,
            'cluster_rep': 'cluster',
            'ref_len': 120,
            'ref_base_assembled': 120,
            'pc_ident': 98.33,
            'ctg': 'ctg_name',
            'ctg_len': 279,
            'known_var': '1',
            'var_type': 'SNP',
            'var_seq_type': 'n',
            'known_var_change': 'A14T',
            'has_known_var': '1',
            'ref_ctg_change': 'A14T',
            'ref_ctg_effect': 'SNP',
            'ref_start': 13,
            'ref_end': 13,
            'ref_nt': 'A',
            'ctg_start': 84,
            'ctg_end': 84,
            'ctg_nt': 'T',
            'smtls_total_depth': 17,
            'smtls_alt_nt': '.',
            'smtls_alt_depth': '17',
            'var_description': 'noncoding1_n_A14T_N_ref has wild type, foo bar',
            'free_text': 'some free text'
        }

        self.assertEqual(summary.Summary._line2dict(line), expected)


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
        lines = [
            'noncoding1\tnon_coding\t19\t78\tnoncoding1\t120\t120\t98.33\tnoncoding1.scaffold.1\t279\t1\tSNP\tn\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1_n_A14T_N_ref has wild type, reads have variant so should report\tgeneric description of noncoding1',
            'noncoding1\tnon_coding\t19\t78\tnoncoding1\t120\t120\t98.33\tnoncoding1.scaffold.1\t279\t1\tSNP\tn\tA6G\t1\t.\t.\t6\t6\tG\t77\t77\tG\t18\t.\t18\tnoncoding1_n_A6G_N_variant in ref and reads so should report\tgeneric description of noncoding1',
            'presence_absence1\tpresence_absence\t27\t88\tpresence_absence1\t96\t96\t98.96\tpresence_absence1.scaffold.1\t267\t1\tSNP\tp\tA10V\t1\tA10V\tNONSYN\t28\t28\tC\t113\t113\tT\t29\t.\t29\tpresence_absence1_p_A10V_N_Ref has wild, reads have variant so report\tGeneric description of presence_absence1',
            'variants_only1\tvariants_only\t27\t64\tvariants_only1\t90\t90\t100.0\tvariants_only1.scaffold.1\t260\t1\tSNP\tp\tS5T\t1\t.\t.\t13\t15\tA;C;C\t96\t98\tA;C;C\t12;13;13\t.;.;.\t12;13;13\tvariants_only1_p_S5T_N_Ref and reads have variant so report\tGeneric description of variants_only1',
        ]

        dicts = [summary.Summary._line2dict(x) for x in lines]
        expected = {
            'noncoding1': {
                ('noncoding1', 'n', 'A14T'): dicts[0],
                ('noncoding1', 'n', 'A6G'): dicts[1],
            },
            'presence_absence1': {('presence_absence1', 'p', 'A10V'): dicts[2]},
            'variants_only1': {('variants_only1', 'p', 'S5T'): dicts[3]}
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
            ['filename', 'noncoding1', 'noncoding1.n.A14T', 'noncoding1.n.A6G', 'presence_absence1', 'presence_absence1.p.A10V', 'variants_only1'],
            [infiles[0], 1, 1, 0, 3, 1, 0],
            [infiles[1], 1, 1, 1, 3, 1, 0],
        ]
        self.assertEqual(expected, got)


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


    def test_write_js_candy_csv(self):
        '''Test _write_js_candy_csv'''
        tmp_out = 'tmp.test_write_js_candy.csv'
        s = summary.Summary(tmp_out, filenames=['spam', 'eggs'])
        s.rows_out = [
            ['filename', 'gene1', 'gene3'],
            ['file1', 1, 3],
            ['file2', 2, 4],
        ]
        s._write_js_candy_csv(tmp_out)
        expected = os.path.join(data_dir, 'summary_test_write_js_candy_csv.csv')
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


    def test_write_js_candy_files(self):
        '''Test _write_js_candy_files'''
        tmp_prefix = 'tmp.test.write_js_candy_files'
        s = summary.Summary('out', filenames=['spam', 'eggs'])
        s.rows_out = [
            ['filename', 'gene1', 'gene2', 'gene3'],
            ['file1', 0, 1, 0],
            ['file2', 1, 0, 3],
            ['file3', 0, 0, 4],
        ]
        s._write_js_candy_files(tmp_prefix)
        expected_csv = os.path.join(data_dir, 'summary_test_write_js_candy_files.csv')
        expected_tre = os.path.join(data_dir, 'summary_test_write_js_candy_files.tre')
        self.assertTrue(filecmp.cmp(expected_csv, tmp_prefix + '.csv', shallow=False))
        self.assertTrue(filecmp.cmp(expected_tre, tmp_prefix + '.tre', shallow=False))
        os.unlink(tmp_prefix + '.csv')
        os.unlink(tmp_prefix + '.tre')

