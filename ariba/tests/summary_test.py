import unittest
import filecmp
import os
from ariba import summary, summary_sample

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
            'assembled,match,ref_seq,pct_id,known_var,novel_var',
            'ref_seq,pct_id,known_var,novel_var',
            'assembled,pct_id,known_var,novel_var',
            'assembled',
            '',
            None,
        ]

        expected = [
            {'assembled': True, 'match': True, 'ref_seq': True, 'pct_id': True, 'known_var': True, 'novel_var': True},
            {'assembled': False, 'match': False, 'ref_seq': True, 'pct_id': True, 'known_var': True, 'novel_var': True},
            {'assembled': True, 'match': False, 'ref_seq': False, 'pct_id': True, 'known_var': True, 'novel_var': True},
            {'assembled': True, 'match': False, 'ref_seq': False, 'pct_id': False, 'known_var': False, 'novel_var': False},
            {'assembled': False, 'match': False, 'ref_seq': False, 'pct_id': False, 'known_var': False, 'novel_var': False},
            {'assembled': False, 'match': False, 'ref_seq': False, 'pct_id': False, 'known_var': False, 'novel_var': False},
        ]

        assert len(col_strings) == len(expected)

        for i in range(len(col_strings)):
            self.assertEqual(expected[i], summary.Summary._determine_cluster_cols(col_strings[i]))


    def test_determine_var_cols(self):
        col_strings = [
            'groups,grouped,ungrouped,novel',
            'groups,grouped,ungrouped',
            'grouped,novel',
            'ungrouped,novel',
            'grouped',
            'ungrouped',
            'novel',
            ''
        ]

        expected = [
            {'groups': True, 'grouped': True, 'ungrouped': True, 'novel': True},
            {'groups': True, 'grouped': True, 'ungrouped': True, 'novel': False},
            {'groups': False, 'grouped': True, 'ungrouped': False, 'novel': True},
            {'groups': False, 'grouped': False, 'ungrouped': True, 'novel': True},
            {'groups': False, 'grouped': True, 'ungrouped': False, 'novel': False},
            {'groups': False, 'grouped': False, 'ungrouped': True, 'novel': False},
            {'groups': False, 'grouped': False, 'ungrouped': False, 'novel': True},
            {'groups': False, 'grouped': False, 'ungrouped': False, 'novel': False},
        ]

        assert len(col_strings) == len(expected)

        for i in range(len(col_strings)):
            self.assertEqual(expected[i], summary.Summary._determine_var_cols(col_strings[i]))


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

        sample1 = summary_sample.SummarySample(file1, only_clusters={'noncoding1'})
        sample2 = summary_sample.SummarySample(file2, only_clusters={'noncoding1'})
        sample1.run()
        sample2.run()
        expected = {file1: sample1, file2: sample2}
        got = summary.Summary._load_input_files([file1, file2], 90, only_clusters={'noncoding1'})
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
            'cluster.p.2': {('presence_absence1', 'A10V', 'grouped', 'id3')},
            'cluster.n.1': {('noncoding1', 'A6G', 'grouped', 'id2'), ('noncoding1', 'A14T', 'grouped', 'id1')},
            'cluster.p.1': {('presence_absence1', 'A10V', 'grouped', 'id2')},
        }
        self.assertEqual(expected, got)


    def test_get_all_het_snps(self):
        '''test _get_all_het_snps'''
        file1 = os.path.join(data_dir, 'summary_test_get_all_het_snps.1.tsv')
        file2 = os.path.join(data_dir, 'summary_test_get_all_het_snps.2.tsv')
        samples = summary.Summary._load_input_files([file1, file2], 90)
        got = summary.Summary._get_all_het_snps(samples)
        expected = {('noncoding1', 'A14T')}
        self.assertEqual(expected, got)


    def test_get_all_var_groups(self):
        '''test _get_all_var_groups'''
        file1 = os.path.join(data_dir, 'summary_test_get_all_var_groups.1.tsv')
        file2 = os.path.join(data_dir, 'summary_test_get_all_var_groups.2.tsv')
        samples = summary.Summary._load_input_files([file1, file2], 90)
        got = summary.Summary._get_all_var_groups(samples)
        expected = {
            'cluster.p.1': {'id4'},
            'cluster.p.2': {'id3'},
            'cluster.v.1': set(),
            'cluster.n.1': {'id1', 'id2'}
        }
        self.assertEqual(expected, got)


    def test_gather_output_rows(self):
        '''Test _gather_output_rows'''
        infiles = [
            os.path.join(data_dir, 'summary_test_gather_output_rows.in.1.tsv'),
            os.path.join(data_dir, 'summary_test_gather_output_rows.in.2.tsv')
        ]
        s = summary.Summary('out', filenames=infiles, variant_cols=None)
        s.samples = summary.Summary._load_input_files(infiles, 90)
        expected = {
            infiles[0]: {
                'noncoding1': {
                    'assembled': 'yes',
                    'match': 'yes',
                    'ref_seq': 'noncoding1',
                    'known_var': 'yes',
                    'novel_var': 'no',
                    'pct_id': '98.33',
                },
                'presence_absence1': {
                    'assembled': 'yes',
                    'match': 'yes',
                    'ref_seq': 'presence_absence1',
                    'known_var': 'no',
                    'novel_var': 'yes',
                    'pct_id': '98.96',
                },
                'variants_only1': {
                    'assembled': 'no',
                    'match': 'no',
                    'ref_seq': 'NA',
                    'known_var': 'NA',
                    'novel_var': 'NA',
                    'pct_id': 'NA',
                }
            },
            infiles[1]: {
                'noncoding1': {
                    'assembled': 'yes',
                    'match': 'yes',
                    'ref_seq': 'noncoding1',
                    'known_var': 'yes',
                    'novel_var': 'no',
                    'pct_id': '98.33',
                },
                'presence_absence1': {
                    'assembled': 'yes',
                    'match': 'yes',
                    'ref_seq': 'presence_absence1',
                    'pct_id': '98.96',
                    'known_var': 'no',
                    'novel_var': 'yes',
                },
                'variants_only1': {
                    'assembled': 'no',
                    'match': 'no',
                    'ref_seq': 'NA',
                    'known_var': 'NA',
                    'novel_var': 'NA',
                    'pct_id': 'NA',
                }
            },
        }
        got = s._gather_output_rows()
        self.assertEqual(expected, got)

        s.var_columns['groups'] = True
        expected[infiles[0]]['noncoding1']['vgroup.id1'] = 'yes'
        expected[infiles[0]]['noncoding1']['vgroup.id3'] = 'no'
        expected[infiles[1]]['noncoding1']['vgroup.id1'] = 'yes'
        expected[infiles[1]]['noncoding1']['vgroup.id3'] = 'yes'
        got = s._gather_output_rows()
        self.assertEqual(expected, got)


        s.var_columns['grouped'] = True
        s.var_columns['ungrouped'] = True
        expected[infiles[0]]['noncoding1']['noncoding1.A14T'] = 'yes'
        expected[infiles[0]]['noncoding1']['noncoding1.A6G'] = 'no'
        expected[infiles[1]]['noncoding1']['noncoding1.A14T'] = 'yes'
        expected[infiles[1]]['noncoding1']['noncoding1.A6G'] = 'yes'
        self.maxDiff = None
        got = s._gather_output_rows()
        self.assertEqual(expected, got)

        s.var_columns['novel'] = True
        expected[infiles[0]]['presence_absence1']['presence_absence1.A10V'] = 'yes'
        expected[infiles[1]]['presence_absence1']['presence_absence1.A10V'] = 'yes'
        got = s._gather_output_rows()
        self.assertEqual(expected, got)

        s.show_known_het = True
        expected[infiles[0]]['noncoding1']['vgroup.id1.%'] = 'NA'
        expected[infiles[0]]['noncoding1']['vgroup.id3.%'] = 'NA'
        expected[infiles[1]]['noncoding1']['vgroup.id1'] = 'het'
        expected[infiles[1]]['noncoding1']['vgroup.id1.%'] = 80.0
        expected[infiles[1]]['noncoding1']['vgroup.id3'] = 'yes'
        expected[infiles[1]]['noncoding1']['vgroup.id3.%'] = 'NA'
        expected[infiles[0]]['noncoding1']['noncoding1.A14T.%'] = 'NA'
        expected[infiles[1]]['noncoding1']['noncoding1.A14T'] = 'het'
        expected[infiles[1]]['noncoding1']['noncoding1.A14T.%'] = 80.0
        got = s._gather_output_rows()
        self.assertEqual(expected, got)

        for filename in expected:
            del expected[filename]['noncoding1']['vgroup.id1']
            del expected[filename]['noncoding1']['vgroup.id3']
            del expected[filename]['noncoding1']['vgroup.id1.%']
            del expected[filename]['noncoding1']['vgroup.id3.%']
            for gene_type in expected[filename]:
                del expected[filename][gene_type]['ref_seq']

        s = summary.Summary('out', filenames=infiles, cluster_cols='assembled,match,pct_id,known_var,novel_var', variant_cols='ungrouped,grouped,novel')
        s.samples = summary.Summary._load_input_files(infiles, 90)
        s.include_all_variant_columns = True
        s.show_known_het = True
        got = s._gather_output_rows()
        self.assertEqual(expected, got)


    def test_to_matrix(self):
        '''Test _to_matrix'''
        rows = {
            'file1': {
                'cluster.n.1': {
                    'assembled': 'yes',
                    'match': 'yes',
                    'ref_seq': 'noncoding1',
                    'known_var': 'yes',
                    'novel_var': 'no',
                    'pct_id': '98.33',
                    'noncoding1.A14T': 'yes'
                },
                'cluster.p.1': {
                    'assembled': 'yes',
                    'match': 'yes',
                    'ref_seq': 'presence_absence1',
                    'known_var': 'yes',
                    'novel_var': 'no',
                    'pct_id': '98.96',
                    'presence_absence1.I42L': 'yes'
                },
                'cluster.v.1': {
                    'assembled': 'yes',
                    'match': 'yes',
                    'ref_seq': 'varonly1',
                    'known_var': 'no',
                    'novel_var': 'no',
                    'pct_id': '99.42',
                }
            },
            'file2': {
                'cluster.n.1': {
                    'assembled': 'yes',
                    'match': 'yes',
                    'ref_seq': 'noncoding1',
                    'known_var': 'no',
                    'novel_var': 'no',
                    'pct_id': '98.33',
                    'noncoding1.A14T': 'no'
                },
                'cluster.p.1': {
                    'assembled': 'yes',
                    'match': 'yes',
                    'ref_seq': 'presence_absence1',
                    'pct_id': '98.96',
                    'known_var': 'no',
                    'novel_var': 'no',
                    'presence_absence1.I42L': 'no'
                },
                'cluster.v.1': {
                    'assembled': 'no',
                    'match': 'NA',
                    'ref_seq': 'NA',
                    'known_var': 'NA',
                    'novel_var': 'NA',
                    'pct_id': 'NA',
                }
            },
        }
        filenames = ['file1', 'file2']
        cluster_cols = {'assembled': True, 'match': True, 'ref_seq': False, 'pct_id': False, 'known_var': False, 'novel_var': False}
        got_phandago_header, got_csv_header, got_lines  = summary.Summary._to_matrix(filenames, rows, cluster_cols)
        expected_phandango_header = ['name', 'cluster.n.1.assembled:o1', 'cluster.n.1.match:o1', 'cluster.n.1.noncoding1.A14T:o1', 'cluster.p.1.assembled:o1', 'cluster.p.1.match:o1', 'cluster.p.1.presence_absence1.I42L:o1', 'cluster.v.1.assembled:o1', 'cluster.v.1.match:o1']
        expected_csv_header = ['name', 'cluster.n.1.assembled', 'cluster.n.1.match', 'cluster.n.1.noncoding1.A14T', 'cluster.p.1.assembled', 'cluster.p.1.match', 'cluster.p.1.presence_absence1.I42L', 'cluster.v.1.assembled', 'cluster.v.1.match']
        expected_lines = [
            ['file1', 'yes', 'yes', 'yes', 'yes', 'yes', 'yes', 'yes', 'yes'],
            ['file2', 'yes', 'yes', 'no', 'yes', 'yes', 'no', 'no', 'NA']
        ]
        self.assertEqual(expected_phandango_header, got_phandago_header)
        self.assertEqual(expected_csv_header, got_csv_header)
        self.assertEqual(expected_lines, got_lines)


    def test_filter_matrix_rows(self):
        '''Test _filter_matrix_rows'''
        matrix = [
            ['yes', 'yes'],
            ['yes', 'no'],
            ['no', 'no'],
            ['yes_nonunique', 'no'],
            ['NA', 'no'],
            ['no', 'NA'],
            ['NA', 'NA']
        ]

        expected = [
            ['yes', 'yes'],
            ['yes', 'no'],
            ['yes_nonunique', 'no'],
        ]
        got = summary.Summary._filter_matrix_rows(matrix)
        self.assertEqual(expected, got)


    def test_filter_matrix_columns(self):
        '''Test _filter_matrix_columns'''
        matrix = [
            ['yes', 'yes', 'no', 'yes_nonunique', 'NA', 'no', 'NA'],
            ['yes', 'no', 'no', 'no', 'no', 'NA', 'NA']
        ]
        phandango_header = ['p1', 'p2', 'p3', 'p4', 'p5', 'p6', 'p7']
        csv_header = ['h1', 'h2', 'h3', 'h4', 'h5', 'h6', 'h7']

        got_phandago_header, got_csv_header, got_matrix  = summary.Summary._filter_matrix_columns(matrix, phandango_header, csv_header)
        expected_phandango_header = ['p1', 'p2', 'p4']
        expected_csv_header = ['h1', 'h2', 'h4']
        expected_matrix = [
            ['yes', 'yes', 'yes_nonunique'],
            ['yes', 'no', 'no'],
        ]
        self.assertEqual(expected_phandango_header, got_phandago_header)
        self.assertEqual(expected_csv_header, got_csv_header)
        self.assertEqual(expected_matrix, got_matrix)


    def test_add_phandango_colour_columns(self):
        '''Test _add_phandango_colour_columns'''
        header = ['head1', 'head2:o1', 'head3:o1', 'head4', 'head5:o1']
        matrix = [
            ['yes', 'yes', 'yes_nonunique', 'yes', 'no'],
            ['yes', 'yes_nonunique', 'no', 'yes', 'NA'],
            ['yes', 'no', 'NA', 'yes', 'yes'],
            ['yes', 'NA', 'yes', 'yes', 'yes_nonunique'],
        ]

        expected_header = ['head1', 'head2', 'head2:colour', 'head3', 'head3:colour', 'head4', 'head5', 'head5:colour']
        expected_matrix = [
            ['yes', 'yes', '#33a02c', 'yes_nonunique', '#b2df8a', 'yes', 'no', '#fb9a99'],
            ['yes', 'yes_nonunique', '#b2df8a', 'no', '#fb9a99', 'yes', 'NA', '#d3d3d3'],
            ['yes', 'no', '#fb9a99', 'NA', '#d3d3d3', 'yes', 'yes', '#33a02c'],
            ['yes', 'NA', '#d3d3d3', 'yes', '#33a02c', 'yes', 'yes_nonunique', '#b2df8a']
        ]
        got_header, got_matrix = summary.Summary._add_phandango_colour_columns(header, matrix)
        self.assertEqual(expected_header, got_header)
        self.assertEqual(expected_matrix, got_matrix)


    def test_matrix_to_csv(self):
        '''Test _matrix_to_csv'''
        matrix = [
            ['line1_1', 'line1_2'],
            ['line2_1', 'line2_2'],
        ]
        header = ['head1', 'head2']
        tmpfile = 'tmp.test.matrix_to_csv.csv'
        summary.Summary._matrix_to_csv(matrix, header, tmpfile)
        with open(tmpfile) as f:
            got = f.read()

        expected = 'head1,head2\nline1_1,line1_2\nline2_1,line2_2\n'
        self.assertEqual(expected, got)
        os.unlink(tmpfile)


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
        # the exact ordering of the nodes is not predictable, so we'll trust dendropy
        # and just check that an output file got written
        self.assertTrue(os.path.exists(tmp_tree))
        os.unlink(tmp_tree)

