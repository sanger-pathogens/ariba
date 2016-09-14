import unittest
import os
from ariba import summary_cluster, summary_cluster_variant


class TestSummaryClusterVariant(unittest.TestCase):
    def test_has_nonsynonymous(self):
        '''Test _has_nonsynonymous'''
        lines = [
            'refname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSYN\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnoncoding1:0:0:A14T:id1:ref has wild type, foo bar\tsome free text',
            'refname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnoncoding1:0:0:A14T:id1:ref has wild type, foo bar\tsome free text',
            'refname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t0\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnoncoding1:0:0:A14T:id1:ref has wild type, foo bar\tsome free text',
            'refname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t0\tSNP\tn\t.\t.\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnoncoding1:0:0:A14T:id1:ref has wild type, foo bar\tsome free text',
            'refname\t1\t0\t528\t2814\tcluster\t1188\t1009\t90.49\tctg_name\t2470\t141.8\t0\t.\tp\t.\t0\t.\tMULTIPLE\t594\t594\tC;T\t1195\t1195\t.;C\t207;204\tC;C\t207;204\t.\t.',
            'refname\t1\t0\t528\t2814\tcluster\t1188\t1009\t90.49\tctg_name\t2470\t141.8\t0\t.\tp\t.\t0\t.\tINDELS\t594\t594\tC;T\t1195\t1195\t.;C\t207;204\tC;C\t207;204\t.\t.'
        ]

        dicts = [summary_cluster.SummaryCluster.line2dict(x) for x in lines]
        expected = [False, True, False, True, True, True]
        assert len(dicts) == len(expected)

        for i in range(len(dicts)):
            self.assertEqual(expected[i], summary_cluster_variant.SummaryClusterVariant._has_nonsynonymous(dicts[i]))


    def test_depths_look_het(self):
        '''test _depths_look_het'''
        tests = [
            ([1], False),
            ([2], False),
            ([3], False),
            ([4], False),
            ([5], False),
            ([90, 1], False),
            ([90, 9], False),
            ([90, 10], True),
            ([9, 1], False),
            ([9, 2], True),
            ([1, 2], True),
            ([90, 5, 5], True),
            ([90, 2, 1, 1], False),
            ([97, 2, 1], False),
        ]

        for depths, expected in tests:
            self.assertEqual(expected, summary_cluster_variant.SummaryClusterVariant._depths_look_het(depths))


    def  test_get_is_het_and_percent(self):
        '''test _get_is_het_and_percent'''
        tests = [
            ('ref1\t0\t0\t531\t78\tcluster1\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnon_coding1:0:0:A14T:id1:foo_bar\tspam eggs', (False, 100.0)),
            ('ref1\t0\t0\t531\t78\tcluster1\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA42T\t1\tA42T\tSNP\t42\t42\tA\t84\t84\tT\t40\tT,A\t10,30\tnon_coding1:0:0:A42T:id1:foo_bar\tspam eggs', (True, 25.0)),
            ('ref1\t0\t0\t531\t78\tcluster1\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA62T\t1\tA62T\tSNP\t62\t62\tA\t84\t84\tA\t40\tA,T\t10,30\tnon_coding1:0:0:A62T:id2:foo_bar\tspam eggs', (True, 75.0)),
            ('ref1\t0\t0\t531\t78\tcluster1\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA82T\t1\tA82T\tSNP\t82\t82\tA\t84\t84\tA\t100\tA,T,G\t10,40,50\tnon_coding1:0:0:A82T:.:foo_bar\tspam eggs', (True, 40.0)),
            ('ref1\t0\t0\t531\t78\tcluster1\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA82T\t1\tA82T\tSNP\t82\t82\tA\t84\t84\tA\t100\tA,T\t95,5\tnon_coding1:0:0:A82T:.:foo_bar\tspam eggs', (False, 5.0)),
            ('ref1\t0\t0\t531\t78\tcluster1\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA82T\t1\tA82T\tSNP\t82\t82\tA\t84\t84\tA\t100\tA,T\t90,10\tnon_coding1:0:0:A82T:.:foo_bar\tspam eggs', (True, 10.0)),
            ('ref1\t0\t0\t531\t78\tcluster1\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA82T\t1\tA82T\tSNP\t82\t82\tA\t84\t84\tA\t100\tA,T,C\t90,6,4\tnon_coding1:0:0:A82T:.:foo_bar\tspam eggs', (True, 6.0)),
            ('ref1\t0\t0\t531\t78\tcluster1\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA82T\t1\tA82T\tSNP\t82\t82\tA\t84\t84\tA\t100\tA,T,C\t3,7,90\tnon_coding1:0:0:A82T:.:foo_bar\tspam eggs', (True, 7.0)),
            ('ref1\t0\t0\t531\t78\tcluster1\t120\t100\t98.33\tctg_name\t279\t24.4\t0\tHET\t.\t.\t.\t.\t.\t.\t.\t.\t84\t84\tA\t50\tA,T\t40,10\t.\t.', (True, 20.0)),
            ('23S.rDNA_WHO_F_01358c\t0\t1\t531\t9914\t23S\t2890\t2890\t99.86\t23S.scaffold.1\t3120\t744.8\t1\tSNP\tn\tC2597T\t1\tC2597T\tSNP\t2597\t2597\tC\t2755\t2755\tT\t823\tTC,T\t487,1\t23S.rDNA_WHO_F_01358c:0:1:C2597T:.:E coli C2611T.\tHigh-level resistance to Azithromycin', (False, 100.0)),
            ('23S.rDNA_WHO_F_01358c\t0\t1\t659\t4168\t23S\t2890\t2890\t99.86\t23S.scaffold.1\t3628\t344.0\t1\tSNP\tn\tC2597T\t0\t.\t.\t2597\t2597\tC\t2928\t2928\tC\t410\tC,T\t90,10\t23S.rDNA_WHO_F_01358c:0:1:C2597T:.:E coli C2611T\t.', (True, 10.0)),
            ('23S.rDNA_WHO_F_01358c\t0\t1\t659\t4168\t23S\t2890\t2890\t99.86\t23S.scaffold.1\t3628\t344.0\t1\tSNP\tn\tC2597T\t0\t.\t.\t2597\t2597\tC\t2928\t2928\tC\t410\tC,T\t91,9\t23S.rDNA_WHO_F_01358c:0:1:C2597T:.:E coli C2611T\t.', (False, 9.0)),
            ('23S.rDNA_WHO_F_01358c\t0\t1\t659\t4168\t23S\t2890\t2890\t99.86\t23S.scaffold.1\t3628\t344.0\t1\tSNP\tn\tC2597T\t0\t.\t.\t2597\t2597\tC\t2928\t2928\tC\t410\tC,T\t50,50\t23S.rDNA_WHO_F_01358c:0:1:C2597T:.:E coli C2611T\t.', (True, 50.0)),
            ('23S.rDNA_WHO_F_01358c\t0\t1\t659\t4168\t23S\t2890\t2890\t99.86\t23S.scaffold.1\t3628\t344.0\t1\tSNP\tn\tC2597T\t0\t.\t.\t2597\t2597\tC\t2928\t2928\tC\t410\tC,T\t70,30\t23S.rDNA_WHO_F_01358c:0:1:C2597T:.:E coli C2611T\t.', (True, 30.0)),
            ('23S.rDNA_WHO_F_01358c\t0\t1\t659\t4168\t23S\t2890\t2890\t99.86\t23S.scaffold.1\t3628\t344.0\t1\tSNP\tn\tC2597T\t1\t.\t.\t2597\t2597\tC\t2928\t2928\tT\t410\tT,C\t91,9\t23S.rDNA_WHO_F_01358c:0:1:C2597T:.:E coli C2611T\t.', (False, 91.0)),
            ('23S.rDNA_WHO_F_01358c\t0\t1\t659\t4168\t23S\t2890\t2890\t99.86\t23S.scaffold.1\t3628\t344.0\t1\tSNP\tn\tC2597T\t1\t.\t.\t2597\t2597\tC\t2928\t2928\tT\t410\tT,C\t90,10\t23S.rDNA_WHO_F_01358c:0:1:C2597T:.:E coli C2611T\t.', (True, 90.0)),
            ('23S.rDNA_WHO_F_01358c\t0\t1\t659\t4168\t23S\t2890\t2890\t99.86\t23S.scaffold.1\t3628\t344.0\t1\tSNP\tn\tC2597T\t1\t.\t.\t2597\t2597\tC\t2928\t2928\tT\t410\tT,C\t50,50\t23S.rDNA_WHO_F_01358c:0:1:C2597T:.:E coli C2611T\t.', (True, 50.0)),
            ('23S.rDNA_WHO_F_01358c\t0\t1\t659\t4168\t23S\t2890\t2890\t99.86\t23S.scaffold.1\t3628\t344.0\t1\tSNP\tn\tC2597T\t1\t.\t.\t2597\t2597\tC\t2928\t2928\tT\t410\tT,C\t10,90\t23S.rDNA_WHO_F_01358c:0:1:C2597T:.:E coli C2611T\t.', (True, 10.0)),
            ('23S.rDNA_WHO_F_01358c\t0\t1\t659\t4168\t23S\t2890\t2890\t99.86\t23S.scaffold.1\t3628\t344.0\t1\tSNP\tn\tC2597T\t1\t.\t.\t2597\t2597\tC\t2928\t2928\tT\t410\tT,C\t9,91\t23S.rDNA_WHO_F_01358c:0:1:C2597T:.:E coli C2611T\t.', (True, 9.0)),
        ]

        for line, expected in tests:
            data_dict = summary_cluster.SummaryCluster.line2dict(line)
            got = summary_cluster_variant.SummaryClusterVariant._get_is_het_and_percent(data_dict)
            self.assertEqual(expected, got)


    def test_init(self):
        '''test __init__'''
        lines = [
            'ref1\t1\t0\t531\t78\tcluster1\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tI14L\t1\tI14L\tSNP\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnon_coding1:0:0:I14L:.:foo_bar\tspam eggs',
            'ref1\t0\t0\t531\t78\tcluster1\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnon_coding1:0:0:A14T:id1:foo_bar\tspam eggs',
            'ref1\t0\t0\t531\t78\tcluster1\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t40\tT,A\t10,30\tnon_coding1:0:0:A14T:id1:foo_bar\tspam eggs',
            'ref1\t0\t0\t531\t78\tcluster1\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t40\tT,A,G\t20,10,10\tnon_coding1:0:0:A14T:id1:foo_bar\tspam eggs',
            'ref1\t0\t0\t531\t78\tcluster1\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t0\t.\t.\t13\t13\tA\t84\t84\tA\t100\tA,T\t90,10\tnon_coding1:0:0:A14T:id1:foo_bar\tspam eggs',
        ]

        expected = [
            {'coding': True, 'known': True, 'var_string': 'I14L', 'var_group': '.', 'het_percent': None},
            {'coding': False, 'known': True, 'var_string': 'A14T', 'var_group': 'id1', 'het_percent': 100.0},
            {'coding': False, 'known': True, 'var_string': 'A14T', 'var_group': 'id1', 'het_percent': 25.0},
            {'coding': False, 'known': True, 'var_string': 'A14T', 'var_group': 'id1', 'het_percent': 50.0},
            {'coding': False, 'known': True, 'var_string': 'A14T', 'var_group': 'id1', 'het_percent': 10.0},
        ]
        assert len(lines) == len(expected)

        for i in range(len(lines)):
            data_dict = summary_cluster.SummaryCluster.line2dict(lines[i])
            cluster_var = summary_cluster_variant.SummaryClusterVariant(data_dict)
            for key in expected[i]:
                got_value = eval('cluster_var.' + key)
                print(i, key)
                self.assertEqual(expected[i][key], got_value)

