import unittest
import os
from ariba import flag, summary_cluster, summary_cluster_variant

modules_dir = os.path.dirname(os.path.abspath(summary_cluster.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')

class TestSummaryCluster(unittest.TestCase):
    def test_line2dict(self):
        '''Test _line2dict'''
        line = 'ariba_refname\trefname\t1\t0\t19\t78\tcluster\t120\t120\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnoncoding1:1:0:A14T:var_group1:ref has wild type, foo bar\tsome free text'

        expected = {
            'ariba_ref_name': 'ariba_refname',
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
            'ctg_cov': 24.4,
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
            'smtls_nts': 'T',
            'smtls_nts_depth': '17',
            'var_description': 'noncoding1:1:0:A14T:var_group1:ref has wild type, foo bar',
            'var_group': 'var_group1',
            'free_text': 'some free text'
        }

        self.assertEqual(summary_cluster.SummaryCluster.line2dict(line), expected)


    def test_add_data_dict(self):
        '''Test add_data_dict'''
        cluster = summary_cluster.SummaryCluster()
        self.assertTrue(cluster.name is None)
        line1 = 'ariba_refname1\trefname\t1\t0\t19\t78\tcluster\t120\t120\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnoncoding1:1:0:A14T:id1:ref has wild type, foo bar\tsome free text'
        line2 = 'ariba_refname1\trefname\t1\t0\t19\t78\tcluster2\t120\t120\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnoncoding1:1:0:A14T:id2:ref has wild type, foo bar\tsome free text'
        line3 = 'ariba_refname2\trefname2\t1\t0\t19\t78\tcluster\t120\t120\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnoncoding1:1:0:A14T:id3:ref has wild type, foo bar\tsome free text'
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


    def test_has_any_part_of_ref_assembled(self):
        '''Test _has_any_part_of_ref_assembled'''
        line1 = 'ariba_refname\trefname\t1\t0\t19\t78\tcluster\t120\t.\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnoncoding1:1:0:A14T:id1:ref has wild type, foo bar\tsome free text'
        line2 = 'ariba_refname\trefname\t1\t0\t19\t78\tcluster\t120\t0\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnoncoding1:1:0:A14T:id1:ref has wild type, foo bar\tsome free text'
        line3 = 'ariba_refname\trefname\t1\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnoncoding1:1:0:A14T:id1:ref has wild type, foo bar\tsome free text'
        data_dict1 = summary_cluster.SummaryCluster.line2dict(line1)
        data_dict2 = summary_cluster.SummaryCluster.line2dict(line2)
        data_dict3 = summary_cluster.SummaryCluster.line2dict(line3)
        cluster = summary_cluster.SummaryCluster()
        cluster.add_data_dict(data_dict1)
        self.assertFalse(cluster._has_any_part_of_ref_assembled())
        cluster.add_data_dict(data_dict2)
        self.assertFalse(cluster._has_any_part_of_ref_assembled())
        cluster.add_data_dict(data_dict3)
        self.assertTrue(cluster._has_any_part_of_ref_assembled())


    def test_pc_id_and_read_depth_of_longest(self):
        '''Test _pc_id_and_read_depth_of_longest'''
        cluster = summary_cluster.SummaryCluster()
        self.assertTrue(cluster.name is None)
        line1 = 'ariba_refname\trefname\t1\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t42.2\t1\tSNP\tn\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnoncoding1:1:0:A14T:id1:ref has wild type, foo bar\tsome free text'
        line2 = 'ariba_refname\trefname\t1\t0\t19\t78\tcluster\t120\t119\t98.20\tctg_name2\t279\t42.42\t1\tSNP\tn\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnoncoding1:1:0:A14T:id1:ref has wild type, foo bar\tsome free text'
        line3 = 'ariba_refname\trefname\t1\t0\t19\t78\tcluster\t120\t114\t98.32\tctg_name3\t279\t42.4\t1\tSNP\tn\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnoncoding1:1:0:A14T:id1:ref has wild type, foo bar\tsome free text'
        data_dict1 = summary_cluster.SummaryCluster.line2dict(line1)
        data_dict2 = summary_cluster.SummaryCluster.line2dict(line2)
        data_dict3 = summary_cluster.SummaryCluster.line2dict(line3)
        cluster.add_data_dict(data_dict1)
        cluster.add_data_dict(data_dict2)
        cluster.add_data_dict(data_dict3)
        self.assertEqual((98.2, 42.42), cluster._pc_id_and_read_depth_of_longest())


    def test_to_cluster_summary_number(self):
        '''Test _to_cluster_summary_assembled'''
        line = 'ariba_refname\trefname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnoncoding1:0:0:A14T:id1:ref has wild type, foo bar\tsome free text'
        data_dict = summary_cluster.SummaryCluster.line2dict(line)

        tests = [
            ('0', 0, 'partial'),
            ('0', 64, 'no'),
            ('0', 1024, 'no'),
            ('0', 1, 'fragmented'),
            ('0', 3, 'yes_nonunique'),
            ('0', 19, 'yes'),
            ('0', 23, 'yes_nonunique'),
            ('0', 51, 'yes_nonunique'),
            ('0', 147, 'yes_nonunique'),
            ('0', 275, 'yes_nonunique'),
            ('1', 0, 'partial'),
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
            if expected == 'partial':
                original_number = cluster.data[0]['ref_base_assembled']
                cluster.data[0]['ref_base_assembled'] = 0
                self.assertEqual('no', cluster._to_cluster_summary_assembled())
                cluster.data[0]['ref_base_assembled'] = original_number


    def test_has_known_variant(self):
        '''Test _has_known_variant'''
        lines = [
            'ariba_refname\trefname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnoncoding1:0:0:A14T:id1:ref has wild type, foo bar\tsome free text',
            'ariba_refname\trefname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t0\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnoncoding1:0:0:A14T:id1:ref has wild type, foo bar\tsome free text',
            'ariba_refname\trefname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t0\tSNP\tn\t.\t.\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnoncoding1:0:0:A14T:id1:ref has wild type, foo bar\tsome free text',
            'ariba_refname\trefname\t1\t0\t528\t2814\tcluster\t1188\t1009\t90.49\tctg_name\t2470\t141.8\t0\t.\tp\t.\t0\t.\tMULTIPLE\t594\t594\tC;T\t1195\t1195\t.;C\t207;204\t.;.\t207;204\t.\t.',
            'ariba_refname\trefname\t1\t0\t528\t2814\tcluster\t1188\t1009\t90.49\tctg_name\t2470\t141.8\t0\t.\tp\t.\t0\t.\tINDELS\t594\t594\tC;T\t1195\t1195\t.;C\t207;204\t.;.\t207;204\t.\t.',
            'ariba_23S.rDNA_WHO_F_01358c\t23S.rDNA_WHO_F_01358c\t0\t1\t659\t4168\t23S\t2890\t2890\t99.86\t23S.scaffold.1\t3628\t344.0\t1\tSNP\tn\tC2597T\t0\t.\t.\t2597\t2597\tC\t2928\t2928\tC\t410\tC,T\t70,30\t23S.rDNA_WHO_F_01358c:0:1:C2597T:.:E coli C2611T\t.',
        ]

        dicts = [summary_cluster.SummaryCluster.line2dict(x) for x in lines]
        expected = ['yes', 'no', 'no', 'no', 'no', 'het']
        assert len(dicts) == len(expected)

        for i in range(len(dicts)):
            self.assertEqual(expected[i], summary_cluster.SummaryCluster._has_known_variant(dicts[i]))


    def test_has_any_known_variant(self):
        lines = [
            'ariba_refname\trefname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnoncoding1:0:0:A14T:id1:ref has wild type, foo bar\tsome free text',
            'ariba_refname\trefname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t0\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnoncoding1:0:0:A14T:id1:ref has wild type, foo bar\tsome free text',
            'ariba_refname\trefname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t0\tSNP\tn\t.\t.\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnoncoding1:0:0:A14T:id1:ref has wild type, foo bar\tsome free text',
            'ariba_refname\trefname\t1\t0\t528\t2814\tcluster\t1188\t1009\t90.49\tctg_name\t2470\t141.8\t0\t.\tp\t.\t0\t.\tMULTIPLE\t594\t594\tC;T\t1195\t1195\t.;C\t207;204\t.;.\t207;204\t.\t.',
            'ariba_refname\trefname\t1\t0\t528\t2814\tcluster\t1188\t1009\t90.49\tctg_name\t2470\t141.8\t0\t.\tp\t.\t0\t.\tINDELS\t594\t594\tC;T\t1195\t1195\t.;C\t207;204\t.;.\t207;204\t.\t.',
            'ariba_23S.rDNA_WHO_F_01358c\t23S.rDNA_WHO_F_01358c\t0\t1\t659\t4168\t23S\t2890\t2890\t99.86\t23S.scaffold.1\t3628\t344.0\t1\tSNP\tn\tC2597T\t0\t.\t.\t2597\t2597\tC\t2928\t2928\tC\t410\tC,T\t70,30\t23S.rDNA_WHO_F_01358c:0:1:C2597T:.:E coli C2611T\t.',
        ]

        expected = ['yes', 'no', 'no', 'no', 'no', 'het']
        assert len(lines) == len(expected)

        for i in range(len(lines)):
            data_dict = summary_cluster.SummaryCluster.line2dict(lines[i])
            cluster = summary_cluster.SummaryCluster()
            cluster.add_data_dict(data_dict)
            self.assertEqual(expected[i], cluster._has_any_known_variant())


    def test_has_nonsynonymous(self):
        '''Test _has_nonsynonymous'''
        lines = [
            'ariba_refname\trefname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSYN\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnoncoding1:0:0:A14T:id1:ref has wild type, foo bar\tsome free text',
            'ariba_refname\trefname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnoncoding1:0:0:A14T:id1:ref has wild type, foo bar\tsome free text',
            'ariba_refname\trefname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t0\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnoncoding1:0:0:A14T:id1:ref has wild type, foo bar\tsome free text',
            'ariba_refname\trefname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t0\tSNP\tn\t.\t.\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnoncoding1:0:0:A14T:id1:ref has wild type, foo bar\tsome free text',
            'ariba_refname\trefname\t1\t0\t528\t2814\tcluster\t1188\t1009\t90.49\tctg_name\t2470\t141.8\t0\t.\tp\t.\t0\t.\tMULTIPLE\t594\t594\tC;T\t1195\t1195\t.;C\t207;204\t.;.\t207;204\t.\t.',
            'ariba_refname\trefname\t1\t0\t528\t2814\tcluster\t1188\t1009\t90.49\tctg_name\t2470\t141.8\t0\t.\tp\t.\t0\t.\tINDELS\t594\t594\tC;T\t1195\t1195\t.;C\t207;204\t.;.\t207;204\t.\t.',
            'ariba_23S.rDNA_WHO_F_01358c\t23S.rDNA_WHO_F_01358c\t0\t1\t659\t4168\t23S\t2890\t2890\t99.86\t23S.scaffold.1\t3628\t344.0\t0\tHET\t.\t.\t.\t.\t.\t2597\t2597\tC\t2928\t2928\tC\t410\tC,T\t70,30\t.\t.',
            'ariba_ref1\tref1\t0\t0\t531\t78\tcluster1\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA62T\t1\tA62T\tSNP\t62\t62\tA\t84\t84\tA\t40\tA,T\t10,30\tnon_coding1:0:0:A62T:id2:foo_bar\tspam eggs',
        ]

        dicts = [summary_cluster.SummaryCluster.line2dict(x) for x in lines]
        expected = ['no', 'yes', 'no', 'yes', 'yes', 'yes', 'het', 'het']
        assert len(dicts) == len(expected)

        for i in range(len(dicts)):
            self.assertEqual(expected[i], summary_cluster.SummaryCluster._has_nonsynonymous(dicts[i]))


    def test_has_any_nonsynonymous(self):
        '''Test _has_any_nonsynonymous'''
        lines = [
            'ariba_refname\trefname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSYN\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnoncoding1:0:0:A14T:N_ref has wild type, foo bar\tsome free text',
            'ariba_refname\trefname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnoncoding1:0:0:A14T:id1:ref has wild type, foo bar\tsome free text',
            'ariba_refname\trefname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t0\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnoncoding1:0:0:A14T:id1:ref has wild type, foo bar\tsome free text',
            'ariba_refname\trefname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t0\tSNP\tn\t.\t.\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnoncoding1:0:0:A14T:id1:ref has wild type, foo bar\tsome free text',
            'ariba_refname\trefname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t0\tSNP\tn\t.\t.\t.\tMULTIPLE\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnoncoding1:0:0:A14T:id1:ref has wild type, foo bar\tsome free text',
            'ariba_23S.rDNA_WHO_F_01358c\t23S.rDNA_WHO_F_01358c\t0\t1\t659\t4168\t23S\t2890\t2890\t99.86\t23S.scaffold.1\t3628\t344.0\t0\tHET\t.\t.\t.\t.\t.\t2597\t2597\tC\t2928\t2928\tC\t410\tC,T\t70,30\t.\t.',
            'ariba_ref1\tref1\t0\t0\t531\t78\tcluster1\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA62T\t1\tA62T\tSNP\t62\t62\tA\t84\t84\tA\t40\tA,T\t10,30\tnon_coding1:0:0:A62T:id2:foo_bar\tspam eggs',
        ]

        expected = ['no', 'yes', 'no', 'yes', 'yes', 'het', 'het']
        assert len(lines) == len(expected)

        for i in range(len(lines)):
            data_dict = summary_cluster.SummaryCluster.line2dict(lines[i])
            cluster = summary_cluster.SummaryCluster()
            cluster.add_data_dict(data_dict)
            self.assertEqual(expected[i], cluster._has_any_nonsynonymous())


    def test_has_novel_nonsynonymous(self):
        '''Test _has_novel_nonsynonymous'''
        lines = [
            'ariba_refname\trefname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnoncoding1:0:0:A14T:id1:ref has wild type, foo bar\tsome free text',
            'ariba_refname\trefname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t0\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnoncoding1:0:0:A14T:id1:ref has wild type, foo bar\tsome free text',
            'ariba_refname\trefname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t0\tSNP\tn\t.\t.\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnoncoding1:0:0:A14T:id1:ref has wild type, foo bar\tsome free text',
            'ariba_refname\trefname\t1\t0\t528\t2814\tcluster\t1188\t1009\t90.49\tctg_name\t2470\t141.8\t0\t.\tp\t.\t0\t.\tMULTIPLE\t594\t594\tC;T\t1195\t1195\t.;C\t207;204\t.;.\t207;204\t.\t.',
            'ariba_refname\trefname\t1\t0\t528\t2814\tcluster\t1188\t1009\t90.49\tctg_name\t2470\t141.8\t0\t.\tp\t.\t0\t.\tINDELS\t594\t594\tC;T\t1195\t1195\t.;C\t207;204\t.;.\t207;204\t.\t.',
            'ariba_23S.rDNA_WHO_F_01358c\t23S.rDNA_WHO_F_01358c\t0\t1\t659\t4168\t23S\t2890\t2890\t99.86\t23S.scaffold.1\t3628\t344.0\t1\tSNP\tn\tC2597T\t0\t.\t.\t2597\t2597\tC\t2928\t2928\tC\t410\tC,T\t70,30\t23S.rDNA_WHO_F_01358c:0:1:C2597T:.:E coli C2611T\t.',
            'ariba_23S.rDNA_WHO_F_01358c\t23S.rDNA_WHO_F_01358c\t0\t1\t659\t4168\t23S\t2890\t2890\t99.86\t23S.scaffold.1\t3628\t344.0\t0\tHET\t.\t.\t.\t.\t.\t2597\t2597\tC\t2928\t2928\tC\t410\tC,T\t70,30\t.\t.',
        ]

        dicts = [summary_cluster.SummaryCluster.line2dict(x) for x in lines]
        expected = ['no', 'no', 'yes', 'yes', 'yes', 'no', 'het']
        assert len(dicts) == len(expected)

        for i in range(len(dicts)):
            self.assertEqual(expected[i], summary_cluster.SummaryCluster._has_novel_nonsynonymous(dicts[i]))


    def test_has_any_novel_nonsynonymous(self):
        '''Test _has_any_novel_nonsynonymous'''
        lines = [
            'ariba_refname\trefname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnoncoding1:0:0:A14T:id1:ref has wild type, foo bar\tsome free text',
            'ariba_refname\trefname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t0\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnoncoding1:0:0:A14T:id1:ref has wild type, foo bar\tsome free text',
            'ariba_refname\trefname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t0\tSNP\tn\t.\t.\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\tT\t17\t.\tsome free text',
            'ariba_refname\trefname\t1\t0\t528\t2814\tcluster\t1188\t1009\t90.49\tctg_name\t2470\t141.8\t0\t.\tp\t.\t0\t.\tMULTIPLE\t594\t594\tC;T\t1195\t1195\t.;C\t207;204\t.;.\t207;204\t.\t.',
            'ariba_refname\trefname\t1\t0\t528\t2814\tcluster\t1188\t1009\t90.49\tctg_name\t2470\t141.8\t0\t.\tp\t.\t0\t.\tINDELS\t594\t594\tC;T\t1195\t1195\t.;C\t207;204\t.;.\t207;204\t.\t.',
            'ariba_23S.rDNA_WHO_F_01358c\t23S.rDNA_WHO_F_01358c\t0\t1\t659\t4168\t23S\t2890\t2890\t99.86\t23S.scaffold.1\t3628\t344.0\t1\tSNP\tn\tC2597T\t0\t.\t.\t2597\t2597\tC\t2928\t2928\tC\t410\tC,T\t70,30\t23S.rDNA_WHO_F_01358c:0:1:C2597T:.:E coli C2611T\t.',
            'ariba_23S.rDNA_WHO_F_01358c\t23S.rDNA_WHO_F_01358c\t0\t1\t659\t4168\t23S\t2890\t2890\t99.86\t23S.scaffold.1\t3628\t344.0\t0\tHET\t.\t.\t.\t.\t.\t2597\t2597\tC\t2928\t2928\tC\t410\tC,T\t70,30\t.\t.',
        ]

        expected = ['no', 'no', 'yes', 'yes', 'yes', 'no', 'het']
        assert len(lines) == len(expected)

        for i in range(len(lines)):
            data_dict = summary_cluster.SummaryCluster.line2dict(lines[i])
            cluster = summary_cluster.SummaryCluster()
            cluster.add_data_dict(data_dict)
            self.assertEqual(expected[i], cluster._has_any_novel_nonsynonymous())


    def test_to_cluster_summary_has_known_nonsynonymous(self):
        '''Test _to_cluster_summary_has_known_nonsynonymous'''
        lines = [
            'ariba_refname\trefname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSYN\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnoncoding1:0:0:A14T:id1:ref has wild type, foo bar\tsome free text',
            'ariba_refname\trefname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnoncoding1:0:0:A14T:id1:ref has wild type, foo bar\tsome free text',
            'ariba_refname\trefname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t0\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnoncoding1:0:0:A14T:id1:ref has wild type, foo bar\tsome free text',
            'ariba_refname\trefname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t0\tSNP\tn\tA14T\t.\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnoncoding1:0:0:A14T:id1:ref has wild type, foo bar\tsome free text',
            'ariba_refname\trefname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t0\t.\tn\t.\t.\t.\tMULTIPLE\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnoncoding1:0:0:A14T:id1:ref has wild type, foo bar\tsome free text',
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
            'ariba_refname\trefname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSYN\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnoncoding1:0:0:A14T:id1:ref has wild type, foo bar\tsome free text',
            'ariba_refname\trefname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnoncoding1:0:0:A14T:id1:ref has wild type, foo bar\tsome free text',
            'ariba_refname\trefname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t0\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnoncoding1:0:0:A14T:id1:ref has wild type, foo bar\tsome free text',
            'ariba_refname\trefname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t0\tSNP\tn\tA14T\t.\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnoncoding1:0:0:A14T:id1:ref has wild type, foo bar\tsome free text',
            'ariba_refname\trefname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t0\t.\tn\t.\t.\t.\tMULTIPLE\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnoncoding1:0:0:A14T:id1:ref has wild type, foo bar\tsome free text',
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
            'ariba_refname\trefname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSYN\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnoncoding1:0:0:A14T:id1:ref has wild type, foo bar\tsome free text',
            'ariba_refname\trefname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnoncoding1:0:0:A14T:id1:ref has wild type, foo bar\tsome free text',
            'ariba_refname\trefname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t0\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnoncoding1:0:0:A14T:id1:ref has wild type, foo bar\tsome free text',
            'ariba_refname\trefname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t0\tSNP\tn\tA14T\t.\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnoncoding1:0:0:A14T:id1:ref has wild type, foo bar\tsome free text',
            'ariba_refname\trefname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t0\t.\tn\t.\t.\t.\tMULTIPLE\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnoncoding1:0:0:A14T:id1:ref has wild type, foo bar\tsome free text',
        ]

        expected = ['no', 'yes', 'no', 'yes', 'yes']

        for i in range(len(lines)):
            data_dict = summary_cluster.SummaryCluster.line2dict(lines[i])
            cluster = summary_cluster.SummaryCluster()
            cluster.add_data_dict(data_dict)
            for assembled_summary in ['yes', 'fragmented', 'yes_nonunique']:
                self.assertEqual(expected[i], cluster._to_cluster_summary_has_nonsynonymous(assembled_summary))
            self.assertEqual('NA', cluster._to_cluster_summary_has_nonsynonymous('no'))


    def test_get_known_noncoding_het_snp(self):
        '''Test _get_known_noncoding_het_snp'''
        d = {'gene': '1'}
        self.assertEqual(None, summary_cluster.SummaryCluster._get_known_noncoding_het_snp(d))
        d['gene'] = '0'
        d['known_var'] = '0'
        self.assertEqual(None, summary_cluster.SummaryCluster._get_known_noncoding_het_snp(d))

        d['known_var'] = '1'
        d['ref_ctg_effect'] = 'not_a_snp'
        self.assertEqual(None, summary_cluster.SummaryCluster._get_known_noncoding_het_snp(d))

        d['ref_ctg_effect'] = 'SNP'
        d['smtls_nts'] = '.'
        self.assertEqual(None, summary_cluster.SummaryCluster._get_known_noncoding_het_snp(d))

        d['smtls_nts'] = 'A;G;T'
        self.assertEqual(None, summary_cluster.SummaryCluster._get_known_noncoding_het_snp(d))

        d['known_var_change'] = 'A42T'
        d['ctg_nt'] = 'A'
        d['smtls_nts'] = 'A,T'
        d['smtls_nts_depth'] = '52,48'
        self.assertEqual(('A42T', 48.0), summary_cluster.SummaryCluster._get_known_noncoding_het_snp(d))


    def test_get_nonsynonymous_var(self):
        '''Test _get_nonsynonymous_var'''
        d = {
            'ref_name': 'ref',
            'gene': '1',
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
        with self.assertRaises(summary_cluster_variant.Error):
            summary_cluster.SummaryCluster._get_nonsynonymous_var(d)

        d['known_var_change'] = 'I42L'
        self.assertEqual(('ref', 'I42L', 'ungrouped', None), summary_cluster.SummaryCluster._get_nonsynonymous_var(d))

        d['var_group'] = 'vgroup'
        self.assertEqual(('ref', 'I42L', 'grouped', 'vgroup'), summary_cluster.SummaryCluster._get_nonsynonymous_var(d))
        d['var_group'] = '.'

        d['ref_ctg_change'] = 'P43Q'
        with self.assertRaises(summary_cluster_variant.Error):
            summary_cluster.SummaryCluster._get_nonsynonymous_var(d)

        d['known_var_change'] = '.'
        self.assertEqual(('ref', 'P43Q', 'novel', None), summary_cluster.SummaryCluster._get_nonsynonymous_var(d))

        d['ref_ctg_change'] = '.'
        with self.assertRaises(summary_cluster_variant.Error):
            summary_cluster.SummaryCluster._get_nonsynonymous_var(d)

        d['ref_ctg_effect'] = 'MULTIPLE'
        self.assertEqual(('ref', 'MULTIPLE', 'novel', None), summary_cluster.SummaryCluster._get_nonsynonymous_var(d))


    def test_has_match(self):
        '''Test _has_match'''
        lines = [
            'ariba_refname\trefname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSYN\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnoncoding1:0:0:A14T:id1:ref has wild type, foo bar\tsome free text',
            'ariba_refname\trefname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t0\tSNP\tn\tA14T\t.\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnoncoding1:0:0:A14T:id1:ref has wild type, foo bar\tsome free text',
            'ariba_refname\trefname\t0\t1\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSYN\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnoncoding1:0:1:A14T:id1:ref has wild type, foo bar\tsome free text',
            'ariba_refname\trefname\t0\t1\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t0\tSNP\tn\tA14T\t.\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnoncoding1:0:1:A14T:id1:ref has wild type, foo bar\tsome free text',
            'ariba_refname\trefname\t1\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tp\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnoncoding1:1:0:A14T:id1:ref has wild type, foo bar\tsome free text',
            'ariba_refname\trefname\t1\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tp\tA14T\t0\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnoncoding1:1:0:A14T:id1:ref has wild type, foo bar\tsome free text',
            'ariba_refname\trefname\t1\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t0\tSNP\tp\tA14T\t.\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnoncoding1:1:0:A14T:id1:ref has wild type, foo bar\tsome free text',
            'ariba_refname\trefname\t1\t1\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tp\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnoncoding1:1:1:A14T:id1:ref has wild type, foo bar\tsome free text',
            'ariba_refname\trefname\t1\t1\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tp\tA14T\t0\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnoncoding1:1:1:A14T:id1:ref has wild type, foo bar\tsome free text',
            'ariba_refname\trefname\t1\t1\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t0\tSNP\tp\tA14T\t.\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnoncoding1:1:1:A14T:id1:ref has wild type, foo bar\tsome free text',
        ]

        expected = ['yes', 'yes', 'yes', 'no', 'yes', 'yes', 'yes', 'yes', 'no', 'no']

        for i in range(len(lines)):
            data_dict = summary_cluster.SummaryCluster.line2dict(lines[i])
            cluster = summary_cluster.SummaryCluster()
            cluster.add_data_dict(data_dict)
            for assembled_summary in ['yes', 'yes_nonunique']:
                self.assertEqual(expected[i], cluster._has_match(assembled_summary))
            for assembled_summary in ['no', 'fragmented']:
                self.assertEqual('no', cluster._has_match(assembled_summary))


    def test_has_var_groups(self):
        '''Test has_var_groups'''
        lines = [
            'ariba_refname\trefname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSYN\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnoncoding1:0:0:A14T:id1:ref has wild type, foo bar\tsome free text',
            'ariba_refname\trefname\t0\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t0\tSNP\tn\tA14T\t.\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnoncoding1:0:0:A14T:id2:ref has wild type, foo bar\tsome free text',
            'ariba_refname\trefname\t1\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tp\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnoncoding1:1:0:A14T:id3:ref has wild type, foo bar\tsome free text',
            'ariba_refname\trefname\t1\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tp\tA14T\t0\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnoncoding1:1:0:A14T:id4:ref has wild type, foo bar\tsome free text',
            'ariba_refname\trefname\t1\t0\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t0\tSNP\tp\tA14T\t.\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnoncoding1:1:0:A14T:id5:ref has wild type, foo bar\tsome free text',
            'ariba_refname\trefname\t1\t1\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tp\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnoncoding1:1:1:A14T:id6:ref has wild type, foo bar\tsome free text',
            'ariba_refname\trefname\t1\t1\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tp\tA14T\t0\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnoncoding1:1:1:A14T:id7:ref has wild type, foo bar\tsome free text',
            'ariba_refname\trefname\t1\t1\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t0\tSNP\tp\tA14T\t.\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnoncoding1:1:1:A14T:id7:ref has wild type, foo bar\tsome free text',
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
        line1 = 'ariba_ref1\tref1\t0\t0\t531\t78\tcluster1\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnon_coding1:n:A14T:id1:foo_bar\tspam eggs'
        line2 = 'ariba_ref1\tref1\t0\t0\t531\t78\tcluster1\t120\t95\t98.42\tctg_name\t279\t24.4\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\tsome free text'

        data_dict1 = summary_cluster.SummaryCluster.line2dict(line1)
        data_dict2 = summary_cluster.SummaryCluster.line2dict(line2)
        cluster = summary_cluster.SummaryCluster()
        cluster.add_data_dict(data_dict1)
        cluster.add_data_dict(data_dict2)
        expected = {
            'assembled': 'yes',
            'match': 'yes',
            'ref_seq': 'ref1',
            'novel_var': 'no',
            'known_var': 'yes',
            'pct_id': '98.33',
            'ctg_cov': '24.4',
        }
        got = cluster.column_summary_data()
        self.assertEqual(expected, got)


    def test_non_synon_variants(self):
        '''Test non_synon_variants'''
        line1 = 'ariba_ref1\tref1\t0\t0\t531\t78\tcluster1\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnon_coding1:0:0:A14T:id1:foo_bar\tspam eggs'
        line2 = 'ariba_ref1\tref1\t0\t0\t531\t78\tcluster1\t120\t95\t98.42\tctg_name\t279\t24.4\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\tsome free text'

        data_dict1 = summary_cluster.SummaryCluster.line2dict(line1)
        data_dict2 = summary_cluster.SummaryCluster.line2dict(line2)
        cluster = summary_cluster.SummaryCluster()
        cluster.add_data_dict(data_dict1)
        cluster.add_data_dict(data_dict2)
        got = cluster.non_synon_variants()
        expected = {('ref1', 'A14T', 'grouped', 'id1')}
        self.assertEqual(expected, got)


    def test_known_noncoding_het_snps(self):
        '''test known_noncoding_het_snps'''
        lines = [
            'ariba_ref1\tref1\t0\t0\t531\t78\tcluster1\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnon_coding1:0:0:A14T:id1:foo_bar\tspam eggs',
            'ariba_ref1\tref1\t0\t0\t531\t78\tcluster1\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA42T\t1\tA42T\tSNP\t42\t42\tA\t84\t84\tT\t40\tT,A\t10,30\tnon_coding1:0:0:A42T:id1:foo_bar\tspam eggs',
            'ariba_ref1\tref1\t0\t0\t531\t78\tcluster1\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA62T\t1\tA62T\tSNP\t62\t62\tA\t84\t84\tA\t40\tA,T\t10,30\tnon_coding1:0:0:A62T:id2:foo_bar\tspam eggs',
            'ariba_ref1\tref1\t0\t0\t531\t78\tcluster1\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA82T\t1\tA82T\tSNP\t82\t82\tA\t84\t84\tA\t100\tA,T,G\t10,40,50\tnon_coding1:0:0:A82T:.:foo_bar\tspam eggs'
        ]

        cluster = summary_cluster.SummaryCluster()
        for line in lines:
            cluster.add_data_dict(summary_cluster.SummaryCluster.line2dict(line))
        got = cluster.known_noncoding_het_snps()
        expected = {
            '.': {'A82T': 40.0},
            'id1': {'A42T': 25.0, 'A14T': 100.0},
            'id2': {'A62T': 75.0},
        }
        self.assertEqual(expected, got)


    def test_get_all_nonsynon_variants_set(self):
        '''test _get_all_nonsynon_variants_set'''
        lines = [
            'ariba_ref1\tref1\t0\t0\t531\t78\tcluster1\t120\t95\t98.42\tctg_name\t279\t24.4\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\tsome free text',
            'ariba_ref1\tref1\t1\t0\t531\t78\tcluster1\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tI14L\t1\tI14L\tSNP\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnon_coding1:0:0:I14L:.:foo_bar\tspam eggs',
            'ariba_ref1\tref1\t0\t0\t531\t78\tcluster1\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t40\tT,A\t10,30\tnon_coding1:0:0:A14T:id1:foo_bar\tspam eggs',
            'ariba_ref1\tref1\t0\t0\t531\t78\tcluster1\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t40\tT,A,G\t20,10,10\tnon_coding1:0:0:A14T:id1:foo_bar\tspam eggs',
        ]

        data_dicts = [summary_cluster.SummaryCluster.line2dict(x) for x in lines]

        cluster_vars = [summary_cluster_variant.SummaryClusterVariant(x) for x in data_dicts]
        expected = {x for x in cluster_vars if x.has_nonsynon}
        got = summary_cluster.SummaryCluster._get_all_nonsynon_variants_set(data_dicts)
        self.assertEqual(expected, got)


    def test_gather_data(self):
        '''test gather_data'''
        lines = [
            'ariba_ref1\tref1\t0\t0\t531\t78\tcluster1\t120\t95\t98.42\tctg_name\t279\t24.4\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\tsome free text',
            'ariba_ref1\tref1\t1\t0\t531\t78\tcluster1\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tI14L\t1\tI14L\tSNP\t13\t13\tA\t84\t84\tT\t17\tT\t17\tnon_coding1:0:0:I14L:.:foo_bar\tspam eggs',
            'ariba_ref1\tref1\t0\t0\t531\t78\tcluster1\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t40\tT,A\t10,30\tnon_coding1:0:0:A14T:id1:foo_bar\tspam eggs',
            'ariba_ref1\tref1\t0\t0\t531\t78\tcluster1\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t40\tT,A,G\t20,10,10\tnon_coding1:0:0:A14T:id1:foo_bar\tspam eggs',
        ]

        data_dicts = [summary_cluster.SummaryCluster.line2dict(x) for x in lines]
        cluster = summary_cluster.SummaryCluster()
        for data_dict in data_dicts:
            cluster.add_data_dict(data_dict)

        cluster.gather_data()
        expected_summary = {
            'assembled': 'yes',
            'match': 'yes',
            'ref_seq': 'ref1',
            'pct_id': '98.33',
            'ctg_cov': '24.4',
            'known_var': 'yes',
            'novel_var': 'no',
        }
        self.assertEqual(expected_summary, cluster.summary)

        cluster_vars = [summary_cluster_variant.SummaryClusterVariant(x) for x in data_dicts]
        expected_variants = {x for x in cluster_vars if x.has_nonsynon}
        self.assertEqual(expected_variants, cluster.variants)

