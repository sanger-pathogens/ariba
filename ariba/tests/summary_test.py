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


    def test_determine_cluster_cols(self):
        col_strings = [
            'assembled,has_res,ref_seq,idty,known_var,novel_var',
            'ref_seq,idty,known_var,novel_var',
            'assembled,idty,known_var,novel_var',
            'assembled',
            '',
            None,
        ]

        expected = [
            {'assembled': True, 'has_res': True, 'ref_seq': True, 'idty': True, 'known_var': True, 'novel_var': True},
            {'assembled': False, 'has_res': False, 'ref_seq': True, 'idty': True, 'known_var': True, 'novel_var': True},
            {'assembled': True, 'has_res': False, 'ref_seq': False, 'idty': True, 'known_var': True, 'novel_var': True},
            {'assembled': True, 'has_res': False, 'ref_seq': False, 'idty': False, 'known_var': False, 'novel_var': False},
            {'assembled': False, 'has_res': False, 'ref_seq': False, 'idty': False, 'known_var': False, 'novel_var': False},
            {'assembled': False, 'has_res': False, 'ref_seq': False, 'idty': False, 'known_var': False, 'novel_var': False},
        ]

        assert len(col_strings) == len(expected)

        for i in range(len(col_strings)):
            self.assertEqual(expected[i], summary.Summary._determine_cluster_cols(col_strings[i]))


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
            'cluster.p.2': {('presence_absence1', 'A10V', 'known')},
            'cluster.n.1': {('noncoding1', 'A6G', 'known'), ('noncoding1', 'A14T', 'known')},
            'cluster.p.1': {('presence_absence1', 'A10V', 'known')},
        }
        self.assertEqual(expected, got)


    def test_gather_output_rows(self):
        '''Test _gather_output_rows'''
        infiles = [
            os.path.join(data_dir, 'summary_test_gather_output_rows.in.1.tsv'),
            os.path.join(data_dir, 'summary_test_gather_output_rows.in.2.tsv')
        ]
        s = summary.Summary('out', filenames=infiles, include_all_known_variant_columns=False)
        s.samples = summary.Summary._load_input_files(infiles, 90)
        expected = {
            infiles[0]: {
                'noncoding1': {
                    'assembled': 'yes',
                    'has_res': 'yes',
                    'ref_seq': 'noncoding1',
                    'known_var': 'yes',
                    'novel_var': 'no',
                    'pct_id': '98.33',
                },
                'presence_absence1': {
                    'assembled': 'yes',
                    'has_res': 'yes',
                    'ref_seq': 'presence_absence1',
                    'known_var': 'no',
                    'novel_var': 'yes',
                    'pct_id': '98.96',
                },
                'variants_only1': {
                    'assembled': 'no',
                    'has_res': 'NA',
                    'ref_seq': 'NA',
                    'known_var': 'NA',
                    'novel_var': 'NA',
                    'pct_id': 'NA',
                }
            },
            infiles[1]: {
                'noncoding1': {
                    'assembled': 'yes',
                    'has_res': 'yes',
                    'ref_seq': 'noncoding1',
                    'known_var': 'yes',
                    'novel_var': 'no',
                    'pct_id': '98.33',
                },
                'presence_absence1': {
                    'assembled': 'yes',
                    'has_res': 'yes',
                    'ref_seq': 'presence_absence1',
                    'pct_id': '98.96',
                    'known_var': 'no',
                    'novel_var': 'yes',
                },
                'variants_only1': {
                    'assembled': 'no',
                    'has_res': 'NA',
                    'ref_seq': 'NA',
                    'known_var': 'NA',
                    'novel_var': 'NA',
                    'pct_id': 'NA',
                }
            },
        }
        got = s._gather_output_rows()
        self.assertEqual(expected, got)

        s.include_all_known_variant_columns = True
        expected[infiles[0]]['noncoding1']['noncoding1.A14T'] = 'yes'
        expected[infiles[0]]['noncoding1']['noncoding1.A6G'] = 'no'
        expected[infiles[1]]['noncoding1']['noncoding1.A14T'] = 'yes'
        expected[infiles[1]]['noncoding1']['noncoding1.A6G'] = 'yes'
        got = s._gather_output_rows()
        self.assertEqual(expected, got)

        s.include_all_novel_variant_columns = True
        expected[infiles[0]]['presence_absence1']['presence_absence1.A10V'] = 'yes'
        expected[infiles[1]]['presence_absence1']['presence_absence1.A10V'] = 'yes'
        got = s._gather_output_rows()
        self.assertEqual(expected, got)

        for filename in expected:
            for gene_type in expected[filename]:
                del expected[filename][gene_type]['ref_seq']

        s = summary.Summary('out', filenames=infiles, cluster_cols='assembled,has_res,idty,known_var,novel_var', include_all_novel_variant_columns=True)
        s.samples = summary.Summary._load_input_files(infiles, 90)
        s.include_all_variant_columns = True
        got = s._gather_output_rows()
        self.assertEqual(expected, got)


    def test_filter_clusters(self):
        '''Test _filter_clusters'''
        rows = {
            'file1': {
                'cluster1': {'assembled': 'yes'},
                'cluster2': {'assembled': 'yes'},
                'cluster3': {'assembled': 'no'},
                'cluster4': {'assembled': 'no'},
            },
            'file2': {
                'cluster1': {'assembled': 'yes'},
                'cluster2': {'assembled': 'no'},
                'cluster3': {'assembled': 'yes'},
                'cluster4': {'assembled': 'no'},
            }
        }

        expected = {
            'file1': {
                'cluster1': {'assembled': 'yes'},
                'cluster2': {'assembled': 'yes'},
                'cluster3': {'assembled': 'no'},
            },
            'file2': {
                'cluster1': {'assembled': 'yes'},
                'cluster2': {'assembled': 'no'},
                'cluster3': {'assembled': 'yes'},
            }
        }

        got = summary.Summary._filter_clusters(rows)
        self.assertEqual((expected, 3), got)


    def test_write_csv(self):
        '''Test _write_csv'''
        tmp_out = 'tmp.out.tsv'
        rows = {
            'file1': {
                'cluster.n.1': {
                    'assembled': 'yes',
                    'ref_seq': 'noncoding1',
                    'known_var': 'yes',
                    'novel_var': 'no',
                    'pct_id': '98.33',
                    'noncoding1.A14T': 'yes'
                },
                'cluster.p.1': {
                    'assembled': 'yes',
                    'ref_seq': 'presence_absence1',
                    'known_var': 'yes',
                    'novel_var': 'no',
                    'pct_id': '98.96',
                    'presence_absence1.I42L': 'yes'
                },
                'cluster.v.1': {
                    'assembled': 'yes',
                    'ref_seq': 'varonly1',
                    'known_var': 'no',
                    'novel_var': 'no',
                    'pct_id': '99.42',
                }
            },
            'file2': {
                'cluster.n.1': {
                    'assembled': 'yes',
                    'ref_seq': 'noncoding1',
                    'known_var': 'no',
                    'novel_var': 'no',
                    'pct_id': '98.33',
                    'noncoding1.A14T': 'no'
                },
                'cluster.p.1': {
                    'assembled': 'yes',
                    'ref_seq': 'presence_absence1',
                    'pct_id': '98.96',
                    'known_var': 'no',
                    'novel_var': 'no',
                    'presence_absence1.I42L': 'no'
                },
                'cluster.v.1': {
                    'assembled': 'no',
                    'ref_seq': 'NA',
                    'known_var': 'NA',
                    'novel_var': 'NA',
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
            'name,cluster.n.1.assembled,cluster.n.1.ref,cluster.n.1.idty,cluster.n.1.known_var,cluster.n.1.novel_var,cluster.n.1.noncoding1.A14T,cluster.p.1.assembled,cluster.p.1.ref,cluster.p.1.idty,cluster.p.1.known_var,cluster.p.1.novel_var,cluster.p.1.presence_absence1.I42L,cluster.v.1.assembled,cluster.v.1.ref,cluster.v.1.idty,cluster.v.1.known_var,cluster.v.1.novel_var',
            'file1,yes,noncoding1,98.33,yes,no,yes,yes,presence_absence1,98.96,yes,no,yes,yes,varonly1,99.42,no,no',
            'file2,yes,noncoding1,98.33,no,no,no,yes,presence_absence1,98.96,no,no,no,no,NA,NA,NA,NA'
]
        expected_lines = [x.split(',') for x in expected_lines]
        self.assertEqual(expected_lines, got_lines)

        got_lines = summary.Summary._write_csv(filenames, rows, tmp_out, phandango=True)
        expected_lines = [
            'name,cluster.n.1.assembled:o1,cluster.n.1.ref:o2,cluster.n.1.idty:c1,cluster.n.1.known_var:o1,cluster.n.1.novel_var:o1,cluster.n.1.noncoding1.A14T:o1,cluster.p.1.assembled:o1,cluster.p.1.ref:o2,cluster.p.1.idty:c1,cluster.p.1.known_var:o1,cluster.p.1.novel_var:o1,cluster.p.1.presence_absence1.I42L:o1,cluster.v.1.assembled:o1,cluster.v.1.ref:o2,cluster.v.1.idty:c1,cluster.v.1.known_var:o1,cluster.v.1.novel_var:o1',
            'file1,yes,noncoding1,98.33,yes,no,yes,yes,presence_absence1,98.96,yes,no,yes,yes,varonly1,99.42,no,no',
            'file2,yes,noncoding1,98.33,no,no,no,yes,presence_absence1,98.96,no,no,no,no,NA,NA,NA,NA'
]
        expected_lines = [x.split(',') for x in expected_lines]
        self.assertEqual(expected_lines, got_lines)
        expected_file = os.path.join(data_dir, 'summary_test_write_tsv.out.phandango.csv')
        self.assertTrue(filecmp.cmp(tmp_out, expected_file, shallow=False))
        os.unlink(tmp_out)


    def test_distance_score_bewteen_values(self):
        '''Test _distance_score_bewteen_values'''
        tests = [
            (('no', 'no'), 0),
            (('no', 'yes'), 1),
            (('no', 'yes_nonunique'), 1),
            (('no', 'fragmented'), 1),
            (('yes', 'no'), 1),
            (('yes', 'yes'), 0),
            (('yes', 'yes_nonunique'), 1),
            (('yes', 'fragmented'), 1),
            (('yes_nonunique', 'no'), 1),
            (('yes_nonunique', 'yes'), 1),
            (('yes_nonunique', 'yes_nonunique'), 0),
            (('yes_nonunique', 'fragmented'), 1),
            (('fragmented', 'no'), 1),
            (('fragmented', 'yes'), 1),
            (('fragmented', 'yes_nonunique'), 1),
            (('fragmented', 'fragmented'), 0),
            (('NA', 'no'), 0),
            (('NA', 'yes'), 1),
            (('NA', 'yes_nonunique'), 1),
            (('NA', 'fragmented'), 1),
        ]

        for (val1, val2), expected in tests:
            self.assertEqual(expected, summary.Summary._distance_score_between_values(val1, val2))
            self.assertEqual(expected, summary.Summary._distance_score_between_values(val2, val1))



    def test_distance_score_between_lists(self):
        '''Test _distance_score_between_lists'''
        list1 = ['NA', 'no', 'yes']
        list2 = ['NA', 'no', 'no']
        self.assertEqual(1, summary.Summary._distance_score_between_lists(list1, list2))


    def test_write_distance_matrix(self):
        '''Test _write_distance_matrix'''
        rows = [
            ['filename', 'gene1', 'gene2', 'gene3'],
            ['file1', 'no', 'yes', 'no'],
            ['file2', 'yes', 'no', 'yes'],
            ['file3', 'no', 'no', 'yes'],
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

