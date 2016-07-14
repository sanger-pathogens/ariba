import unittest
import copy
import filecmp
import os
from ariba import flag, summary_cluster

modules_dir = os.path.dirname(os.path.abspath(summary_cluster.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')

class TestSummaryCluster(unittest.TestCase):
    def test_line2dict(self):
        '''Test _line2dict'''
        line = 'refname\t1\t0\t19\t78\tcluster\t120\t120\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1:n:A14T:var_group1:ref has wild type, foo bar\tsome free text'

        expected = {
            'ref_name': 'refname',
            'gene': '1',
            'var_only' : '0',
            'flag': flag.Flag(19),
            'reads': 78,
            'cluster': 'cluster',
            'ref_len': 120,
            'ref_base_assembled': 120,
            'pc_ident': 98.33,
            'ctg': 'ctg_name',
            'ctg_len': 279,
            'ctg_cov': '24.4',
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
            'smtls_total_depth': '17',
            'smtls_alt_nt': '.',
            'smtls_alt_depth': '17',
            'var_description': 'noncoding1:n:A14T:var_group1:ref has wild type, foo bar',
            'var_group': 'var_group1',
            'free_text': 'some free text'
        }

        self.assertEqual(summary_cluster.SummaryCluster.line2dict(line), expected)


    def test_add_data_dict(self):
        '''Test add_data_dict'''
        cluster = summary_cluster.SummaryCluster()
        self.assertTrue(cluster.name is None)
        line1 = 'refname\t1\t0\t19\t78\tcluster\t120\t120\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1:n:A14T:id1:ref has wild type, foo bar\tsome free text'
        line2 = 'refname\t1\t0\t19\t78\tcluster2\t120\t120\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1:n:A14T:id2:ref has wild type, foo bar\tsome free text'
        line3 = 'refname2\t1\t0\t19\t78\tcluster\t120\t120\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1:n:A14T:id3:ref has wild type, foo bar\tsome free text'
        data_dict1 = summary_cluster.SummaryCluster.line2dict(line1)
        data_dict2 = summary_cluster.SummaryCluster.line2dict(line2)
        data_dict3 = summary_cluster.SummaryCluster.line2dict(line3)
        cluster.add_data_dict(data_dict1)
        self.assertEqual(cluster.name, data_dict1['cluster'])
        self.assertEqual(cluster.data,[data_dict1])
        with self.assertRaises(summary_cluster.Error):
            cluster.add_data_dict(data_dict2)

        with self.assertRaises(summary_cluster.Error):
            cluster.add_data_dict(data_dict3)


    def test_pc_id_of_longest(self):
        '''Test pc_id_of_longest'''
        cluster = summary_cluster.SummaryCluster()
        self.assertTrue(cluster.name is None)
        line1 = 'refname\t1\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1:n:A14T:id1:ref has wild type, foo bar\tsome free text'
        line2 = 'refname\t1\t0\t19\t78\tcluster\t120\t119\t98.20\tctg_name2\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1:n:A14T:id1:ref has wild type, foo bar\tsome free text'
        line3 = 'refname\t1\t0\t19\t78\tcluster\t120\t114\t98.32\tctg_name3\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1:n:A14T:id1:ref has wild type, foo bar\tsome free text'
        data_dict1 = summary_cluster.SummaryCluster.line2dict(line1)
        data_dict2 = summary_cluster.SummaryCluster.line2dict(line2)
        data_dict3 = summary_cluster.SummaryCluster.line2dict(line3)
        cluster.add_data_dict(data_dict1)
        cluster.add_data_dict(data_dict2)
        cluster.add_data_dict(data_dict3)
        self.assertEqual(98.2, cluster.pc_id_of_longest())


    def test_to_cluster_summary_number(self):
        '''Test _to_cluster_summary_assembled'''
        line = 'refname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1:n:A14T:id1:ref has wild type, foo bar\tsome free text'
        data_dict = summary_cluster.SummaryCluster.line2dict(line)

        tests = [
            ('0', 0, 'no'),
            ('0', 64, 'no'),
            ('0', 1024, 'no'),
            ('0', 1, 'fragmented'),
            ('0', 3, 'yes_nonunique'),
            ('0', 19, 'yes'),
            ('0', 23, 'yes_nonunique'),
            ('0', 51, 'yes_nonunique'),
            ('0', 147, 'yes_nonunique'),
            ('0', 275, 'yes_nonunique'),
            ('1', 0, 'no'),
            ('1', 64, 'no'),
            ('1', 1024, 'no'),
            ('1', 1, 'fragmented'),
            ('1', 11, 'yes_nonunique'),
            ('1', 27, 'yes'),
            ('1', 29, 'fragmented'),
            ('1', 59, 'yes_nonunique'),
            ('1', 155, 'yes_nonunique'),
            ('1', 283, 'yes_nonunique'),
        ]

        for gene, f, expected in tests:
            cluster = summary_cluster.SummaryCluster()
            data_dict['gene'] = gene
            data_dict['flag'] = flag.Flag(f)
            cluster.add_data_dict(data_dict)
            self.assertEqual(expected, cluster._to_cluster_summary_assembled())


    def test_has_known_variant(self):
        '''Test _has_known_variant'''
        lines = [
            'refname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1:n:A14T:id1:ref has wild type, foo bar\tsome free text',
            'refname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t0\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1:n:A14T:id1:ref has wild type, foo bar\tsome free text',
            'refname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t0\tSNP\tn\t.\t.\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1:n:A14T:id1:ref has wild type, foo bar\tsome free text',
            'refname\t1\t0\t528\t2814\tcluster\t1188\t1009\t90.49\tctg_name\t2470\t141.8\t0\t.\tp\t.\t0\t.\tMULTIPLE\t594\t594\tC;T\t1195\t1195\t.;C\t207;204\t.;.\t207;204\t.\t.',
            'refname\t1\t0\t528\t2814\tcluster\t1188\t1009\t90.49\tctg_name\t2470\t141.8\t0\t.\tp\t.\t0\t.\tINDELS\t594\t594\tC;T\t1195\t1195\t.;C\t207;204\t.;.\t207;204\t.\t.'
        ]

        dicts = [summary_cluster.SummaryCluster.line2dict(x) for x in lines]
        expected = [True, False, False, False, False]
        assert len(dicts) == len(expected)

        for i in range(len(dicts)):
            self.assertEqual(expected[i], summary_cluster.SummaryCluster._has_known_variant(dicts[i]))


    def test_has_any_known_variant(self):
        lines = [
            'refname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1:n:A14T:id1:ref has wild type, foo bar\tsome free text',
            'refname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t0\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1:n:A14T:id1:ref has wild type, foo bar\tsome free text',
            'refname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t0\tSNP\tn\t.\t.\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1:n:A14T:id1:ref has wild type, foo bar\tsome free text',
            'refname\t1\t0\t528\t2814\tcluster\t1188\t1009\t90.49\tctg_name\t2470\t141.8\t0\t.\tp\t.\t0\t.\tMULTIPLE\t594\t594\tC;T\t1195\t1195\t.;C\t207;204\t.;.\t207;204\t.\t.',
            'refname\t1\t0\t528\t2814\tcluster\t1188\t1009\t90.49\tctg_name\t2470\t141.8\t0\t.\tp\t.\t0\t.\tINDELS\t594\t594\tC;T\t1195\t1195\t.;C\t207;204\t.;.\t207;204\t.\t.'
        ]

        expected = ['yes', 'no', 'no', 'no', 'no']
        assert len(lines) == len(expected)

        for i in range(len(lines)):
            data_dict = summary_cluster.SummaryCluster.line2dict(lines[i])
            cluster = summary_cluster.SummaryCluster()
            cluster.add_data_dict(data_dict)
            self.assertEqual(expected[i], cluster._has_any_known_variant())


    def test_has_nonsynonymous(self):
        '''Test _has_nonsynonymous'''
        lines = [
            'refname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSYN\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1:n:A14T:id1:ref has wild type, foo bar\tsome free text',
            'refname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1:n:A14T:id1:ref has wild type, foo bar\tsome free text',
            'refname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t0\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1:n:A14T:id1:ref has wild type, foo bar\tsome free text',
            'refname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t0\tSNP\tn\t.\t.\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1:n:A14T:id1:ref has wild type, foo bar\tsome free text',
            'refname\t1\t0\t528\t2814\tcluster\t1188\t1009\t90.49\tctg_name\t2470\t141.8\t0\t.\tp\t.\t0\t.\tMULTIPLE\t594\t594\tC;T\t1195\t1195\t.;C\t207;204\t.;.\t207;204\t.\t.',
            'refname\t1\t0\t528\t2814\tcluster\t1188\t1009\t90.49\tctg_name\t2470\t141.8\t0\t.\tp\t.\t0\t.\tINDELS\t594\t594\tC;T\t1195\t1195\t.;C\t207;204\t.;.\t207;204\t.\t.'
        ]

        dicts = [summary_cluster.SummaryCluster.line2dict(x) for x in lines]
        expected = [False, True, False, True, True, True]
        assert len(dicts) == len(expected)

        for i in range(len(dicts)):
            self.assertEqual(expected[i], summary_cluster.SummaryCluster._has_nonsynonymous(dicts[i]))


    def test_has_any_nonsynonymous(self):
        '''Test _has_any_nonsynonymous'''
        lines = [
            'refname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSYN\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1:n:A14T:N_ref has wild type, foo bar\tsome free text',
            'refname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1:n:A14T:id1:ref has wild type, foo bar\tsome free text',
            'refname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t0\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1:n:A14T:id1:ref has wild type, foo bar\tsome free text',
            'refname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t0\tSNP\tn\t.\t.\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1:n:A14T:id1:ref has wild type, foo bar\tsome free text',
            'refname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t0\tSNP\tn\t.\t.\t.\tMULTIPLE\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1:n:A14T:id1:ref has wild type, foo bar\tsome free text',
        ]

        expected = ['no', 'yes', 'no', 'yes', 'yes']
        assert len(lines) == len(expected)

        for i in range(len(lines)):
            data_dict = summary_cluster.SummaryCluster.line2dict(lines[i])
            cluster = summary_cluster.SummaryCluster()
            cluster.add_data_dict(data_dict)
            self.assertEqual(expected[i], cluster._has_any_nonsynonymous())


    def test_has_novel_nonsynonymous(self):
        '''Test _has_novel_nonsynonymous'''
        lines = [
            'refname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1:n:A14T:id1:ref has wild type, foo bar\tsome free text',
            'refname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t0\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1:n:A14T:id1:ref has wild type, foo bar\tsome free text',
            'refname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t0\tSNP\tn\t.\t.\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1:n:A14T:id1:ref has wild type, foo bar\tsome free text',
            'refname\t1\t0\t528\t2814\tcluster\t1188\t1009\t90.49\tctg_name\t2470\t141.8\t0\t.\tp\t.\t0\t.\tMULTIPLE\t594\t594\tC;T\t1195\t1195\t.;C\t207;204\t.;.\t207;204\t.\t.',
            'refname\t1\t0\t528\t2814\tcluster\t1188\t1009\t90.49\tctg_name\t2470\t141.8\t0\t.\tp\t.\t0\t.\tINDELS\t594\t594\tC;T\t1195\t1195\t.;C\t207;204\t.;.\t207;204\t.\t.'
        ]

        dicts = [summary_cluster.SummaryCluster.line2dict(x) for x in lines]
        expected = [False, False, True, True, True]
        assert len(dicts) == len(expected)

        for i in range(len(dicts)-1):
            self.assertEqual(expected[i], summary_cluster.SummaryCluster._has_novel_nonsynonymous(dicts[i]))


    def test_has_any_novel_nonsynonymous(self):
        '''Test _has_any_novel_nonsynonymous'''
        lines = [
            'refname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1:n:A14T:id1:ref has wild type, foo bar\tsome free text',
            'refname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t0\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1:n:A14T:id1:ref has wild type, foo bar\tsome free text',
            'refname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t0\tSNP\tn\t.\t.\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1:n:A14T:id1:ref has wild type, foo bar\tsome free text',
            'refname\t1\t0\t528\t2814\tcluster\t1188\t1009\t90.49\tctg_name\t2470\t141.8\t0\t.\tp\t.\t0\t.\tMULTIPLE\t594\t594\tC;T\t1195\t1195\t.;C\t207;204\t.;.\t207;204\t.\t.',
            'refname\t1\t0\t528\t2814\tcluster\t1188\t1009\t90.49\tctg_name\t2470\t141.8\t0\t.\tp\t.\t0\t.\tINDELS\t594\t594\tC;T\t1195\t1195\t.;C\t207;204\t.;.\t207;204\t.\t.'
        ]

        expected = ['no', 'no', 'yes', 'yes', 'yes']
        assert len(lines) == len(expected)

        for i in range(len(lines)):
            data_dict = summary_cluster.SummaryCluster.line2dict(lines[i])
            cluster = summary_cluster.SummaryCluster()
            cluster.add_data_dict(data_dict)
            self.assertEqual(expected[i], cluster._has_any_novel_nonsynonymous())


    def test_to_cluster_summary_has_known_nonsynonymous(self):
        '''Test _to_cluster_summary_has_known_nonsynonymous'''
        lines = [
            'refname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSYN\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1:n:A14T:id1:ref has wild type, foo bar\tsome free text',
            'refname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1:n:A14T:id1:ref has wild type, foo bar\tsome free text',
            'refname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t0\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1:n:A14T:id1:ref has wild type, foo bar\tsome free text',
            'refname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t0\tSNP\tn\tA14T\t.\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1:n:A14T:id1:ref has wild type, foo bar\tsome free text',
            'refname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t0\t.\tn\t.\t.\t.\tMULTIPLE\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1:n:A14T:id1:ref has wild type, foo bar\tsome free text',
        ]

        expected = ['yes', 'yes', 'no', 'no', 'no']

        for i in range(len(lines)):
            data_dict = summary_cluster.SummaryCluster.line2dict(lines[i])
            cluster = summary_cluster.SummaryCluster()
            cluster.add_data_dict(data_dict)
            for assembled_summary in ['yes', 'fragmented', 'yes_nonunique']:
                self.assertEqual(expected[i], cluster._to_cluster_summary_has_known_nonsynonymous(assembled_summary))
            self.assertEqual('NA', cluster._to_cluster_summary_has_known_nonsynonymous('no'))


    def test_to_cluster_summary_has_novel_nonsynonymous(self):
        '''Test _to_cluster_summary_has_novel_nonsynonymous'''
        lines = [
            'refname\tnon_coding\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSYN\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1:n:A14T:id1:ref has wild type, foo bar\tsome free text',
            'refname\tnon_coding\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1:n:A14T:id1:ref has wild type, foo bar\tsome free text',
            'refname\tnon_coding\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t0\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1:n:A14T:id1:ref has wild type, foo bar\tsome free text',
            'refname\tnon_coding\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t0\tSNP\tn\tA14T\t.\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1:n:A14T:id1:ref has wild type, foo bar\tsome free text',
            'refname\tnon_coding\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t0\t.\tn\t.\t.\t.\tMULTIPLE\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1:n:A14T:id1:ref has wild type, foo bar\tsome free text',
        ]

        expected = ['no', 'no', 'no', 'yes', 'yes']

        for i in range(len(lines)):
            data_dict = summary_cluster.SummaryCluster.line2dict(lines[i])
            cluster = summary_cluster.SummaryCluster()
            cluster.add_data_dict(data_dict)
            for assembled_summary in ['yes', 'fragmented', 'yes_nonunique']:
                self.assertEqual(expected[i], cluster._to_cluster_summary_has_novel_nonsynonymous(assembled_summary))
            self.assertEqual('NA', cluster._to_cluster_summary_has_novel_nonsynonymous('no'))


    def test_to_cluster_summary_has_nonsynonymous(self):
        '''Test _to_cluster_summary_has_nonsynonymous'''
        lines = [
            'refname\tnon_coding\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSYN\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1:n:A14T:id1:ref has wild type, foo bar\tsome free text',
            'refname\tnon_coding\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1:n:A14T:id1:ref has wild type, foo bar\tsome free text',
            'refname\tnon_coding\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t0\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1:n:A14T:id1:ref has wild type, foo bar\tsome free text',
            'refname\tnon_coding\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t0\tSNP\tn\tA14T\t.\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1:n:A14T:id1:ref has wild type, foo bar\tsome free text',
            'refname\tnon_coding\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t0\t.\tn\t.\t.\t.\tMULTIPLE\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1:n:A14T:id1:ref has wild type, foo bar\tsome free text',
        ]

        expected = ['no', 'yes', 'no', 'yes', 'yes']

        for i in range(len(lines)):
            data_dict = summary_cluster.SummaryCluster.line2dict(lines[i])
            cluster = summary_cluster.SummaryCluster()
            cluster.add_data_dict(data_dict)
            for assembled_summary in ['yes', 'fragmented', 'yes_nonunique']:
                self.assertEqual(expected[i], cluster._to_cluster_summary_has_nonsynonymous(assembled_summary))
            self.assertEqual('NA', cluster._to_cluster_summary_has_nonsynonymous('no'))


    def test_get_nonsynonymous_var(self):
        '''Test _get_nonsynonymous_var'''
        d = {
            'ref_name': 'ref',
            'var_type': '.',
            'known_var_change': '.',
            'has_known_var': '.',
            'known_var': '0',
            'ref_ctg_change': '.',
            'ref_ctg_effect': '.',
            'var_seq_type': '.',
            'var_group': '.',
        }

        self.assertEqual(None, summary_cluster.SummaryCluster._get_nonsynonymous_var(d))

        d['var_type'] = 'p'
        d['known_var'] = '1'
        d['has_known_var'] = '1'
        with self.assertRaises(summary_cluster.Error):
            summary_cluster.SummaryCluster._get_nonsynonymous_var(d)

        d['known_var_change'] = 'I42L'
        self.assertEqual(('ref', 'I42L', 'ungrouped', None), summary_cluster.SummaryCluster._get_nonsynonymous_var(d))

        d['var_group'] = 'vgroup'
        self.assertEqual(('ref', 'I42L', 'grouped', 'vgroup'), summary_cluster.SummaryCluster._get_nonsynonymous_var(d))
        d['var_group'] = '.'

        d['ref_ctg_change'] = 'P43Q'
        with self.assertRaises(summary_cluster.Error):
            summary_cluster.SummaryCluster._get_nonsynonymous_var(d)

        d['known_var_change'] = '.'
        self.assertEqual(('ref', 'P43Q', 'novel', None), summary_cluster.SummaryCluster._get_nonsynonymous_var(d))

        d['ref_ctg_change'] = '.'
        with self.assertRaises(summary_cluster.Error):
            summary_cluster.SummaryCluster._get_nonsynonymous_var(d)

        d['ref_ctg_effect'] = 'MULTIPLE'
        self.assertEqual(('ref', 'MULTIPLE', 'novel', None), summary_cluster.SummaryCluster._get_nonsynonymous_var(d))


    def test_has_resistance(self):
        '''Test _has_resistance'''
        lines = [
            'refname\tnon_coding\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSYN\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1:n:A14T:id1:ref has wild type, foo bar\tsome free text',
            'refname\tnon_coding\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t0\tSNP\tn\tA14T\t.\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1:n:A14T:id1:ref has wild type, foo bar\tsome free text',
            'refname\tpresence_absence\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tp\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1:n:A14T:id1:ref has wild type, foo bar\tsome free text',
            'refname\tpresence_absence\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tp\tA14T\t0\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1:n:A14T:id1:ref has wild type, foo bar\tsome free text',
            'refname\tpresence_absence\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t0\tSNP\tp\tA14T\t.\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1:n:A14T:id1:ref has wild type, foo bar\tsome free text',
            'refname\tvariants_only\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tp\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1:n:A14T:id1:ref has wild type, foo bar\tsome free text',
            'refname\tvariants_only\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tp\tA14T\t0\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1:n:A14T:id1:ref has wild type, foo bar\tsome free text',
            'refname\tvariants_only\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t0\tSNP\tp\tA14T\t.\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1:n:A14T:id1:ref has wild type, foo bar\tsome free text',
        ]

        expected = ['yes', 'yes', 'yes', 'yes', 'yes', 'yes', 'no', 'no']

        for i in range(len(lines)):
            data_dict = summary_cluster.SummaryCluster.line2dict(lines[i])
            cluster = summary_cluster.SummaryCluster()
            cluster.add_data_dict(data_dict)
            for assembled_summary in ['yes', 'yes_nonunique']:
                self.assertEqual(expected[i], cluster._has_resistance(assembled_summary))
            for assembled_summary in ['no', 'fragmented']:
                self.assertEqual('no', cluster._has_resistance(assembled_summary))


    def test_has_var_groups(self):
        '''Test has_var_groups'''
        lines = [
            'refname\tnon_coding\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSYN\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1:n:A14T:id1:ref has wild type, foo bar\tsome free text',
            'refname\tnon_coding\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t0\tSNP\tn\tA14T\t.\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1:n:A14T:id2:ref has wild type, foo bar\tsome free text',
            'refname\tpresence_absence\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tp\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1:n:A14T:id3:ref has wild type, foo bar\tsome free text',
            'refname\tpresence_absence\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tp\tA14T\t0\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1:n:A14T:id4:ref has wild type, foo bar\tsome free text',
            'refname\tpresence_absence\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t0\tSNP\tp\tA14T\t.\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1:n:A14T:id5:ref has wild type, foo bar\tsome free text',
            'refname\tvariants_only\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tp\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1:n:A14T:id6:ref has wild type, foo bar\tsome free text',
            'refname\tvariants_only\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tp\tA14T\t0\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1:n:A14T:id7:ref has wild type, foo bar\tsome free text',
            'refname\tvariants_only\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t0\tSNP\tp\tA14T\t.\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1:n:A14T:id7:ref has wild type, foo bar\tsome free text',
        ]
        dicts = [summary_cluster.SummaryCluster.line2dict(line) for line in lines]
        cluster = summary_cluster.SummaryCluster()
        for d in dicts:
            cluster.add_data_dict(d)
        got = cluster.has_var_groups()
        expected = {'id1', 'id3', 'id6'}
        self.assertEqual(expected, got)


    def test_column_summary_data(self):
        '''Test column_summary_data'''
        line1 = 'ref1\tnon_coding\t531\t78\tcluster1\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnon_coding1:n:A14T:id1:foo_bar\tspam eggs'
        line2 = 'ref1\tnon_coding\t531\t78\tcluster1\t120\t95\t98.42\tctg_name\t279\t24.4\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\tsome free text'

        data_dict1 = summary_cluster.SummaryCluster.line2dict(line1)
        data_dict2 = summary_cluster.SummaryCluster.line2dict(line2)
        cluster = summary_cluster.SummaryCluster()
        cluster.add_data_dict(data_dict1)
        cluster.add_data_dict(data_dict2)
        expected = {
            'assembled': 'yes',
            'has_res': 'yes',
            'ref_seq': 'ref1',
            'novel_var': 'no',
            'known_var': 'yes',
            'pct_id': '98.33',
        }
        got = cluster.column_summary_data()
        self.assertEqual(expected, got)


    def test_non_synon_variants(self):
        '''Test non_synon_variants'''
        line1 = 'ref1\tnon_coding\t531\t78\tcluster1\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnon_coding1:n:A14T:id1:foo_bar\tspam eggs'
        line2 = 'ref1\tnon_coding\t531\t78\tcluster1\t120\t95\t98.42\tctg_name\t279\t24.4\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\tsome free text'

        data_dict1 = summary_cluster.SummaryCluster.line2dict(line1)
        data_dict2 = summary_cluster.SummaryCluster.line2dict(line2)
        cluster = summary_cluster.SummaryCluster()
        cluster.add_data_dict(data_dict1)
        cluster.add_data_dict(data_dict2)
        got = cluster.non_synon_variants()
        expected = {('ref1', 'A14T', 'grouped', 'id1')}
        self.assertEqual(expected, got)
