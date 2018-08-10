import unittest
from ariba import summary_cluster, summary_cluster_variant


class TestSummaryClusterVariant(unittest.TestCase):
    def test_has_nonsynonymous(self):
        '''Test _has_nonsynonymous'''
        lines = [
            'ariba_refname\trefname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSYN\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnoncoding1:0:0:A14T:id1:ref has wild type, foo bar\tsome free text',
            'ariba_refname\trefname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnoncoding1:0:0:A14T:id1:ref has wild type, foo bar\tsome free text',
            'ariba_refname\trefname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t0\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnoncoding1:0:0:A14T:id1:ref has wild type, foo bar\tsome free text',
            'ariba_refname\trefname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t0\tSNP\tn\t.\t.\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnoncoding1:0:0:A14T:id1:ref has wild type, foo bar\tsome free text',
            'ariba_refname\trefname\t1\t0\t528\t2814\tcluster\t1188\t1009\t90.49\tctg_name\t2470\t141.8\t0\t.\tp\t.\t0\t.\tMULTIPLE\t594\t594\tC;T\t1195\t1195\t.;C\t207;204\tC;C\t207;204\t.\t.',
            'ariba_refname\trefname\t1\t0\t528\t2814\tcluster\t1188\t1009\t90.49\tctg_name\t2470\t141.8\t0\t.\tp\t.\t0\t.\tINDELS\t594\t594\tC;T\t1195\t1195\t.;C\t207;204\tC;C\t207;204\t.\t.'
        ]

        dicts = [summary_cluster.SummaryCluster.line2dict(x) for x in lines]
        expected = [False, True, False, True, True, True]
        assert len(dicts) == len(expected)

        for i in range(len(dicts)):
            self.assertEqual(expected[i], summary_cluster_variant.SummaryClusterVariant._has_nonsynonymous(dicts[i]))


    def test_filter_depths(self):
        '''test _filter_depths'''
        tests = [
            ('A', {'A': 1}, {'A': 1}),
            ('A', {'A': 2}, {'A': 2}),
            ('A', {'A': 3}, {'A': 3}),
            ('A', {'A': 4}, {'A': 4}),
            ('A', {'A': 5}, {'A': 5}),
            ('A', {'A': 90, 'C': 9}, {'A': 90}),
            ('C', {'A': 90, 'C': 9}, {'A': 90, 'C': 9}),
            ('C', {'A': 90, 'C': 9, 'G':1}, {'A': 90, 'C': 9}),
            ('A', {'A': 90, 'C': 10}, {'A': 90, 'C': 10}),
            ('A', {'A': 90, 'C': 5, 'G': 5}, {'A': 90}),
            ('A', {'A': 89, 'C': 10, 'G': 1}, {'A': 89, 'C': 10}),
            ('A', {'A': 80, 'C': 10, 'G': 10}, {'A': 80, 'C': 10, 'G': 10}),
        ]

        for ref_base, depths, expected in tests:
            self.assertEqual(expected, summary_cluster_variant.SummaryClusterVariant._filter_depths(ref_base, depths))


    #def test_depths_look_het(self):
    #    '''test _depths_look_het'''
    #    tests = [
    #        ([1], False),
    #        ([2], False),
    #        ([3], False),
    #        ([4], False),
    #        ([5], False),
    #        ([90, 1], False),
    #        ([90, 9], False),
    #        ([90, 10], True),
    #        ([9, 1], False),
    #        ([9, 2], True),
    #        ([1, 2], True),
    #        ([89, 10, 1], True),
    #        ([89, 9, 2], False),
    #        ([90, 2, 1, 1], False),
    #        ([97, 2, 1], False),
    #    ]

    #    for depths, expected in tests:
    #        self.assertEqual(expected, summary_cluster_variant.SummaryClusterVariant._depths_look_het(depths))


    def  test_get_is_het_and_percent(self):
        '''test _get_is_het_and_percent'''
        tests = [
            ('ariba_ref1\tref1\t0\t0\t531\t78\tcluster1\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnon_coding1:0:0:A14T:id1:foo_bar\tspam eggs', (False, 100.0, 'T')),
            ('ariba_ref1\tref1\t0\t0\t531\t78\tcluster1\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA42T\t1\tA42T\tSNP\t42\t42\tA\t84\t84\tT\t40\tT,A\t10,30\tnon_coding1:0:0:A42T:id1:foo_bar\tspam eggs', (True, 25.0, 'AT')),
            ('ariba_ref1\tref1\t0\t0\t531\t78\tcluster1\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA62T\t1\tA62T\tSNP\t62\t62\tA\t84\t84\tA\t40\tA,T\t10,30\tnon_coding1:0:0:A62T:id2:foo_bar\tspam eggs', (True, 75.0, 'AT')),
            ('ariba_ref1\tref1\t0\t0\t531\t78\tcluster1\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA82T\t1\tA82T\tSNP\t82\t82\tA\t84\t84\tA\t100\tA,T,G\t10,40,50\tnon_coding1:0:0:A82T:.:foo_bar\tspam eggs', (True, 40.0, 'AGT')),
            ('ariba_ref1\tref1\t0\t0\t531\t78\tcluster1\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA82T\t1\tA82T\tSNP\t82\t82\tA\t84\t84\tA\t100\tA,T\t95,5\tnon_coding1:0:0:A82T:.:foo_bar\tspam eggs', (False, 5.0, 'A')),
            ('ariba_ref1\tref1\t0\t0\t531\t78\tcluster1\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA82T\t1\tA82T\tSNP\t82\t82\tA\t84\t84\tA\t100\tA,T\t90,10\tnon_coding1:0:0:A82T:.:foo_bar\tspam eggs', (True, 10.0, 'AT')),
            ('ariba_ref1\tref1\t0\t0\t531\t78\tcluster1\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA82T\t1\tA82T\tSNP\t82\t82\tA\t84\t84\tA\t100\tA,T,C\t90,6,4\tnon_coding1:0:0:A82T:.:foo_bar\tspam eggs', (False, 6.0, 'A')),
            ('ariba_ref1\tref1\t0\t0\t531\t78\tcluster1\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA82T\t1\tA82T\tSNP\t82\t82\tA\t84\t84\tA\t100\tA,T,C\t3,7,90\tnon_coding1:0:0:A82T:.:foo_bar\tspam eggs', (True, 7.0, 'ACT')),
            ('ariba_ref1\tref1\t0\t0\t531\t78\tcluster1\t120\t100\t98.33\tctg_name\t279\t24.4\t0\tHET\t.\t.\t.\t.\t.\t.\t.\t.\t84\t84\tA\t50\tA,T\t40,10\t.\t.', (True, 20.0, 'AT')),
            ('ariba_ref1\t23S.rDNA_WHO_F_01358c\t0\t1\t531\t9914\t23S\t2890\t2890\t99.86\t23S.scaffold.1\t3120\t744.8\t1\tSNP\tn\tC2597T\t1\tC2597T\tSNP\t2597\t2597\tC\t2755\t2755\tT\t823\tC,T\t487,1\t23S.rDNA_WHO_F_01358c:0:1:C2597T:.:E coli C2611T.\tHigh-level resistance to Azithromycin', (False, 0.2, 'C')),
            ('ariba\t23S.rDNA_WHO_F_01358c\t0\t1\t659\t4168\t23S\t2890\t2890\t99.86\t23S.scaffold.1\t3628\t344.0\t1\tSNP\tn\tC2597T\t0\t.\t.\t2597\t2597\tC\t2928\t2928\tC\t410\tC,T\t90,10\t23S.rDNA_WHO_F_01358c:0:1:C2597T:.:E coli C2611T\t.', (True, 10.0, 'CT')),
            ('ariba_23S.rDNA_WHO_F_01358c\t23S.rDNA_WHO_F_01358c\t0\t1\t659\t4168\t23S\t2890\t2890\t99.86\t23S.scaffold.1\t3628\t344.0\t1\tSNP\tn\tC2597T\t0\t.\t.\t2597\t2597\tC\t2928\t2928\tC\t410\tC,T\t91,9\t23S.rDNA_WHO_F_01358c:0:1:C2597T:.:E coli C2611T\t.', (False, 9.0, 'C')),
            ('ariba_23S.rDNA_WHO_F_01358c\t23S.rDNA_WHO_F_01358c\t0\t1\t659\t4168\t23S\t2890\t2890\t99.86\t23S.scaffold.1\t3628\t344.0\t1\tSNP\tn\tC2597T\t0\t.\t.\t2597\t2597\tC\t2928\t2928\tC\t410\tC,T\t50,50\t23S.rDNA_WHO_F_01358c:0:1:C2597T:.:E coli C2611T\t.', (True, 50.0, 'CT')),
            ('ariba_23S.rDNA_WHO_F_01358c\t23S.rDNA_WHO_F_01358c\t0\t1\t659\t4168\t23S\t2890\t2890\t99.86\t23S.scaffold.1\t3628\t344.0\t1\tSNP\tn\tC2597T\t0\t.\t.\t2597\t2597\tC\t2928\t2928\tC\t410\tC,T\t70,30\t23S.rDNA_WHO_F_01358c:0:1:C2597T:.:E coli C2611T\t.', (True, 30.0, 'CT')),
            ('ariba_23S.rDNA_WHO_F_01358c\t23S.rDNA_WHO_F_01358c\t0\t1\t659\t4168\t23S\t2890\t2890\t99.86\t23S.scaffold.1\t3628\t344.0\t1\tSNP\tn\tC2597T\t1\t.\t.\t2597\t2597\tC\t2928\t2928\tT\t410\tT,C\t91,9\t23S.rDNA_WHO_F_01358c:0:1:C2597T:.:E coli C2611T\t.', (False, 91.0, 'T')),
            ('ariba_23S.rDNA_WHO_F_01358c\t23S.rDNA_WHO_F_01358c\t0\t1\t659\t4168\t23S\t2890\t2890\t99.86\t23S.scaffold.1\t3628\t344.0\t1\tSNP\tn\tC2597T\t1\t.\t.\t2597\t2597\tC\t2928\t2928\tT\t410\tT,C\t90,10\t23S.rDNA_WHO_F_01358c:0:1:C2597T:.:E coli C2611T\t.', (True, 90.0, 'CT')),
            ('ariba_23S.rDNA_WHO_F_01358c\t23S.rDNA_WHO_F_01358c\t0\t1\t659\t4168\t23S\t2890\t2890\t99.86\t23S.scaffold.1\t3628\t344.0\t1\tSNP\tn\tC2597T\t1\t.\t.\t2597\t2597\tC\t2928\t2928\tT\t410\tT,C\t50,50\t23S.rDNA_WHO_F_01358c:0:1:C2597T:.:E coli C2611T\t.', (True, 50.0, 'CT')),
            ('ariba_23S.rDNA_WHO_F_01358c\t23S.rDNA_WHO_F_01358c\t0\t1\t659\t4168\t23S\t2890\t2890\t99.86\t23S.scaffold.1\t3628\t344.0\t1\tSNP\tn\tC2597T\t1\t.\t.\t2597\t2597\tC\t2928\t2928\tT\t410\tT,C\t10,90\t23S.rDNA_WHO_F_01358c:0:1:C2597T:.:E coli C2611T\t.', (True, 10.0, 'CT')),
            ('ariba_23S.rDNA_WHO_F_01358c\t23S.rDNA_WHO_F_01358c\t0\t1\t659\t4168\t23S\t2890\t2890\t99.86\t23S.scaffold.1\t3628\t344.0\t1\tSNP\tn\tC2597T\t1\t.\t.\t2597\t2597\tC\t2928\t2928\tT\t410\tT,C\t9,91\t23S.rDNA_WHO_F_01358c:0:1:C2597T:.:E coli C2611T\t.', (True, 9.0, 'CT')),
        ]

        for line, expected in tests:
            data_dict = summary_cluster.SummaryCluster.line2dict(line)
            got = summary_cluster_variant.SummaryClusterVariant._get_is_het_and_percent(data_dict)
            self.assertEqual(expected, got)


    def test_init(self):
        '''test __init__'''
        lines = [
            'ariba_ref1\tref1\t1\t0\t531\t78\tcluster1\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tI14L\t1\tI14L\tSNP\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnon_coding1:0:0:I14L:.:foo_bar\tspam eggs',
            'ariba_ref1\tref1\t0\t0\t531\t78\tcluster1\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnon_coding1:0:0:A14T:id1:foo_bar\tspam eggs',
            'ariba_ref1\tref1\t0\t0\t531\t78\tcluster1\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t40\tT,A\t10,30\tnon_coding1:0:0:A14T:id1:foo_bar\tspam eggs',
            'ariba_ref1\tref1\t0\t0\t531\t78\tcluster1\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t40\tT,A,G\t20,10,10\tnon_coding1:0:0:A14T:id1:foo_bar\tspam eggs',
            'ariba_ref1\tref1\t0\t0\t531\t78\tcluster1\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t0\t.\t.\t13\t13\tA\t84\t84\tA\t100\tA,T\t90,10\tnon_coding1:0:0:A14T:id1:foo_bar\tspam eggs',
        ]

        expected = [
            {'coding': True, 'known': True, 'var_string': 'I14L', 'var_group': '.', 'het_percent': None},
            {'coding': False, 'known': True, 'var_string': '14T', 'var_group': 'id1', 'het_percent': 100.0},
            {'coding': False, 'known': True, 'var_string': '14AT', 'var_group': 'id1', 'het_percent': 25.0},
            {'coding': False, 'known': True, 'var_string': '14AGT', 'var_group': 'id1', 'het_percent': 50.0},
            {'coding': False, 'known': True, 'var_string': '14AT', 'var_group': 'id1', 'het_percent': 10.0},
        ]
        assert len(lines) == len(expected)

        for i in range(len(lines)):
            data_dict = summary_cluster.SummaryCluster.line2dict(lines[i])
            cluster_var = summary_cluster_variant.SummaryClusterVariant(data_dict)
            for key in expected[i]:
                got_value = eval('cluster_var.' + key)
                self.assertEqual(expected[i][key], got_value)

