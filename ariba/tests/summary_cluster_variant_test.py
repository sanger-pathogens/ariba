import unittest
import os
from ariba import summary_cluster, summary_cluster_variant


class TestSummaryClusterVariant(unittest.TestCase):
    def test_has_nonsynonymous(self):
        '''Test _has_nonsynonymous'''
        lines = [
            'refname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSYN\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1:0:0:A14T:id1:ref has wild type, foo bar\tsome free text',
            'refname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1:0:0:A14T:id1:ref has wild type, foo bar\tsome free text',
            'refname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t0\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1:0:0:A14T:id1:ref has wild type, foo bar\tsome free text',
            'refname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t0\tSNP\tn\t.\t.\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1:0:0:A14T:id1:ref has wild type, foo bar\tsome free text',
            'refname\t1\t0\t528\t2814\tcluster\t1188\t1009\t90.49\tctg_name\t2470\t141.8\t0\t.\tp\t.\t0\t.\tMULTIPLE\t594\t594\tC;T\t1195\t1195\t.;C\t207;204\t.;.\t207;204\t.\t.',
            'refname\t1\t0\t528\t2814\tcluster\t1188\t1009\t90.49\tctg_name\t2470\t141.8\t0\t.\tp\t.\t0\t.\tINDELS\t594\t594\tC;T\t1195\t1195\t.;C\t207;204\t.;.\t207;204\t.\t.'
        ]

        dicts = [summary_cluster.SummaryCluster.line2dict(x) for x in lines]
        expected = [False, True, False, True, True, True]
        assert len(dicts) == len(expected)

        for i in range(len(dicts)):
            self.assertEqual(expected[i], summary_cluster_variant.SummaryClusterVariant._has_nonsynonymous(dicts[i]))


    def test_depths_look_het(self):
        '''test _depths_look_het'''
        tests = [
            (['1'], False),
            (['2'], False),
            (['3'], False),
            (['4'], False),
            (['5'], False),
            (['90', '1'], False),
            (['90', '9'], False),
            (['90', '10'], True),
            (['9', '1'], False),
            (['9', '2'], True),
            (['1', '2'], True),
            (['90', '5', '5'], True),
            (['90', '2', '1', '1'], False),
            (['97', '2', '1'], False),
        ]

        for depths, expected in tests:
            self.assertEqual(expected, summary_cluster_variant.SummaryClusterVariant._depths_look_het(depths))


    def  test_get_het_percent(self):
        '''test _get_het_percent'''
        lines = [
            'ref1\t0\t0\t531\t78\tcluster1\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnon_coding1:0:0:A14T:id1:foo_bar\tspam eggs',
            'ref1\t0\t0\t531\t78\tcluster1\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA42T\t1\tA42T\tSNP\t42\t42\tA\t84\t84\tT\t40\tA\t10,30\tnon_coding1:0:0:A42T:id1:foo_bar\tspam eggs',
            'ref1\t0\t0\t531\t78\tcluster1\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA62T\t1\tA62T\tSNP\t62\t62\tA\t84\t84\tA\t40\tT\t10,30\tnon_coding1:0:0:A62T:id2:foo_bar\tspam eggs',
            'ref1\t0\t0\t531\t78\tcluster1\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA82T\t1\tA82T\tSNP\t82\t82\tA\t84\t84\tA\t100\tT,G\t10,40,50\tnon_coding1:0:0:A82T:.:foo_bar\tspam eggs',
            'ref1\t0\t0\t531\t78\tcluster1\t120\t100\t98.33\tctg_name\t279\t24.4\t0\tHET\t.\t.\t.\t.\t.\t.\t.\t.\t84\t84\tA\t50\tT\t40,10\t.\t.'
        ]

        expected = [None, 25.0, 75.0, 40.0, 20.0]
        assert len(lines) == len(expected)

        for i in range(len(lines)):
            data_dict = summary_cluster.SummaryCluster.line2dict(lines[i])
            got = summary_cluster_variant.SummaryClusterVariant._get_het_percent(data_dict)
            self.assertEqual(expected[i], got)


    def test_init(self):
        '''test __init__'''
        lines = [
            'ref1\t1\t0\t531\t78\tcluster1\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tI14L\t1\tI14L\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnon_coding1:0:0:I14L:.:foo_bar\tspam eggs',
            'ref1\t0\t0\t531\t78\tcluster1\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnon_coding1:0:0:A14T:id1:foo_bar\tspam eggs',
            'ref1\t0\t0\t531\t78\tcluster1\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t40\tA\t10,30\tnon_coding1:0:0:A14T:id1:foo_bar\tspam eggs',
            'ref1\t0\t0\t531\t78\tcluster1\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t40\tA,G\t20,10,10\tnon_coding1:0:0:A14T:id1:foo_bar\tspam eggs',
        ]

        expected = [
            {'coding': True, 'known': True, 'var_string': 'I14L', 'var_group': '.', 'het_percent': None},
            {'coding': False, 'known': True, 'var_string': 'A14T', 'var_group': 'id1', 'het_percent': None},
            {'coding': False, 'known': True, 'var_string': 'A14T', 'var_group': 'id1', 'het_percent': 25.0},
            {'coding': False, 'known': True, 'var_string': 'A14T', 'var_group': 'id1', 'het_percent': 50.0},
        ]
        assert len(lines) == len(expected)

        for i in range(len(lines)):
            data_dict = summary_cluster.SummaryCluster.line2dict(lines[i])
            cluster_var = summary_cluster_variant.SummaryClusterVariant(data_dict)
            for key in expected[i]:
                got_value = eval('cluster_var.' + key)
                self.assertEqual(expected[i][key], got_value)

