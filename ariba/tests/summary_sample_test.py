import unittest
import os
from ariba import summary_cluster, summary_sample

modules_dir = os.path.dirname(os.path.abspath(summary_sample.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


class TestSummarySample(unittest.TestCase):
    def test_load_file(self):
        '''Test _load_file'''
        infile = os.path.join(data_dir, 'summary_sample_test_load_file.in.tsv')
        with open(infile) as f:
             lines = [x.rstrip() for x in f]

        dicts = [summary_cluster.SummaryCluster.line2dict(x) for x in lines[1:]]
        cluster1 = summary_cluster.SummaryCluster()
        cluster1.add_data_dict(dicts[0])
        cluster1.add_data_dict(dicts[1])
        cluster1.add_data_dict(dicts[2])
        cluster1.gather_data()
        cluster2 = summary_cluster.SummaryCluster()
        cluster2.add_data_dict(dicts[3])
        cluster2.add_data_dict(dicts[4])
        cluster2.gather_data()
        cluster3 = summary_cluster.SummaryCluster()
        cluster3.add_data_dict(dicts[5])
        cluster3.gather_data()

        expected = {
            'cluster.n': cluster1,
            'cluster.p': cluster2,
            'cluster.v': cluster3
        }

        got = summary_sample.SummarySample._load_file(infile, 90)
        self.assertEqual(expected, got)

        got = summary_sample.SummarySample._load_file(infile, 90, only_clusters={'cluster.n'})
        expected = {'cluster.n': cluster1}
        self.assertEqual(expected, got)

    def test_column_summary_data(self):
        '''Test _column_summary_data'''
        infile = os.path.join(data_dir, 'summary_sample_test_column_summary_data.tsv')
        sample_summary = summary_sample.SummarySample(infile)
        sample_summary.clusters = sample_summary._load_file(infile, 90)
        expected = {
            'cluster.n': {
                'assembled': 'yes',
                'match': 'yes',
                'ref_seq': 'noncoding1',
                'known_var': 'yes',
                'novel_var': 'yes',
                'pct_id': '98.33',
                'ctg_cov': '35.4',
            },
            'cluster.p': {
                'assembled': 'yes',
                'match': 'yes',
                'ref_seq': 'presence_absence1',
                'known_var': 'yes',
                'novel_var': 'no',
                'pct_id': '98.96',
                'ctg_cov': '35.1',
            },
            'cluster.v': {
                'assembled': 'yes',
                'match': 'yes',
                'ref_seq': 'variants_only1',
                'known_var': 'yes',
                'novel_var': 'no',
                'pct_id': '100.0',
                'ctg_cov': '42.4',
            }
        }
        self.maxDiff = None
        got = sample_summary._column_summary_data()
        self.assertEqual(expected, got)


    def test_var_groups(self):
        '''test _var_groups'''
        infile = os.path.join(data_dir, 'summary_sample_test_var_groups.tsv')
        sample_summary = summary_sample.SummarySample(infile)
        sample_summary.clusters = sample_summary._load_file(infile, 90)
        got = sample_summary._var_groups()
        expected = {
            'cluster.n': {'id1', 'id2'},
            'cluster.p': {'id3'},
            'cluster.v': {'id4'}
        }
        self.assertEqual(expected, got)


    def test_variant_column_names_tuples_and_het_snps(self):
        '''Test _variant_column_names_tuples_and_het_snps'''
        infile = os.path.join(data_dir, 'summary_sample_test_column_names_tuples_and_het_snps.tsv')
        sample_summary = summary_sample.SummarySample(infile)
        sample_summary.clusters = sample_summary._load_file(infile, 90)
        sample_summary.column_summary_data = sample_summary._column_summary_data()
        expected_column_names = {
            'cluster.v': {('variants_only1', 'S5T', 'ungrouped', None)},
            'cluster.n': {
                ('noncoding1', 'A6G', 'grouped', 'id2'),
                ('noncoding1', 'A14T', 'ungrouped', None),
                ('noncoding1', 'G15T', 'novel', None)
             },
            'cluster.p': {('presence_absence1', 'A10V', 'grouped', 'id3')}
        }
        got_column_names, got_het_snps = sample_summary._variant_column_names_tuples_and_het_snps()
        self.assertEqual(expected_column_names, got_column_names)

        expected_het_snps = {
            'cluster.v': {},
            'cluster.n': {'.': {'A14T': 80.0}},
            'cluster.p': {},
        }
        self.assertEqual(expected_het_snps, got_het_snps)

