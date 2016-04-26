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
        line = 'refname\treftype\t19\t78\tcluster\t120\t120\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1_n_A14T_N_ref has wild type, foo bar\tsome free text'

        expected = {
            'ref_name': 'refname',
            'ref_type': 'reftype',
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
            'var_description': 'noncoding1_n_A14T_N_ref has wild type, foo bar',
            'free_text': 'some free text'
        }

        self.assertEqual(summary_cluster.SummaryCluster.line2dict(line), expected)


    def test_add_data_dict(self):
        '''Test add_data_dict'''
        cluster = summary_cluster.SummaryCluster()
        self.assertTrue(cluster.name is None)
        line1 = 'refname\treftype\t19\t78\tcluster\t120\t120\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1_n_A14T_N_ref has wild type, foo bar\tsome free text'
        line2 = 'refname\treftype\t19\t78\tcluster2\t120\t120\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1_n_A14T_N_ref has wild type, foo bar\tsome free text'
        line3 = 'refname2\treftype\t19\t78\tcluster\t120\t120\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1_n_A14T_N_ref has wild type, foo bar\tsome free text'
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
        line1 = 'refname\treftype\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1_n_A14T_N_ref has wild type, foo bar\tsome free text'
        line2 = 'refname\treftype\t19\t78\tcluster\t120\t119\t98.20\tctg_name2\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1_n_A14T_N_ref has wild type, foo bar\tsome free text'
        line3 = 'refname\treftype\t19\t78\tcluster\t120\t114\t98.32\tctg_name3\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1_n_A14T_N_ref has wild type, foo bar\tsome free text'
        data_dict1 = summary_cluster.SummaryCluster.line2dict(line1)
        data_dict2 = summary_cluster.SummaryCluster.line2dict(line2)
        data_dict3 = summary_cluster.SummaryCluster.line2dict(line3)
        cluster.add_data_dict(data_dict1)
        cluster.add_data_dict(data_dict2)
        cluster.add_data_dict(data_dict3)
        self.assertEqual(98.2, cluster.pc_id_of_longest())


    def test_to_cluster_summary_number(self):
        '''Test _to_cluster_summary_number_assembled'''
        line = 'refname\tnon_coding\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1_n_A14T_N_ref has wild type, foo bar\tsome free text'
        data_dict = summary_cluster.SummaryCluster.line2dict(line)

        tests = [
            ('non_coding', 0, 'no'),
            ('non_coding', 64, 'no'),
            ('non_coding', 1024, 'no'),
            ('non_coding', 1, 'fragmented'),
            ('non_coding', 3, 'yes_nonunique'),
            ('non_coding', 19, 'yes'),
            ('non_coding', 23, 'yes_nonunique'),
            ('non_coding', 51, 'yes_nonunique'),
            ('non_coding', 147, 'yes_nonunique'),
            ('non_coding', 275, 'yes_nonunique'),
            ('presence_absence', 0, 'no'),
            ('presence_absence', 64, 'no'),
            ('presence_absence', 1024, 'no'),
            ('presence_absence', 1, 'fragmented'),
            ('presence_absence', 11, 'yes_nonunique'),
            ('presence_absence', 27, 'yes'),
            ('presence_absence', 29, 'fragmented'),
            ('presence_absence', 59, 'yes_nonunique'),
            ('presence_absence', 155, 'yes_nonunique'),
            ('presence_absence', 283, 'yes_nonunique'),
        ]

        for seq_type, f, expected in tests:
            cluster = summary_cluster.SummaryCluster()
            data_dict['ref_type'] = seq_type
            data_dict['flag'] = flag.Flag(f)
            cluster.add_data_dict(data_dict)
            self.assertEqual(expected, cluster._to_cluster_summary_number_assembled())


    def test_has_nonsynonymous(self):
        '''Test _has_nonsynonymous'''
        lines = [
            'refname\tnon_coding\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSYN\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1_n_A14T_N_ref has wild type, foo bar\tsome free text',
            'refname\tnon_coding\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1_n_A14T_N_ref has wild type, foo bar\tsome free text',
            'refname\tnon_coding\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t0\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1_n_A14T_N_ref has wild type, foo bar\tsome free text',
            'refname\tnon_coding\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t0\tSNP\tn\tA14T\t.\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1_n_A14T_N_ref has wild type, foo bar\tsome free text',
            'refname\tnon_coding\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t0\tSNP\tn\tA14T\t.\t.\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1_n_A14T_N_ref has wild type, foo bar\tsome free text',
        ]

        dicts = [summary_cluster.SummaryCluster.line2dict(x) for x in lines]
        expected = [False, True, False, True, False]
        assert len(dicts) == len(expected)

        for i in range(len(dicts)):
            self.assertEqual(expected[i], summary_cluster.SummaryCluster._has_nonsynonymous(dicts[i]))


    def test_has_any_nonsynonymous(self):
        '''Test _has_any_nonsynonymous'''
        lines = [
            'refname\tnon_coding\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSYN\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1_n_A14T_N_ref has wild type, foo bar\tsome free text',
            'refname\tnon_coding\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1_n_A14T_N_ref has wild type, foo bar\tsome free text',
            'refname\tnon_coding\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t0\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1_n_A14T_N_ref has wild type, foo bar\tsome free text',
            'refname\tnon_coding\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t0\tSNP\tn\tA14T\t.\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1_n_A14T_N_ref has wild type, foo bar\tsome free text',
            'refname\tnon_coding\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t0\tSNP\tn\tA14T\t.\t.\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1_n_A14T_N_ref has wild type, foo bar\tsome free text',
        ]

        expected = ['no', 'yes', 'no', 'yes', 'no']
        assert len(lines) == len(expected)

        for i in range(len(lines)):
            data_dict = summary_cluster.SummaryCluster.line2dict(lines[i])
            cluster = summary_cluster.SummaryCluster()
            cluster.add_data_dict(data_dict)
            self.assertEqual(expected[i], cluster._has_any_nonsynonymous())


    def test_to_cluster_summary_number_has_nonsynonymous(self):
        '''Test _to_cluster_summary_number_has_nonsynonymous'''
        lines = [
            'refname\tnon_coding\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSYN\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1_n_A14T_N_ref has wild type, foo bar\tsome free text',
            'refname\tnon_coding\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1_n_A14T_N_ref has wild type, foo bar\tsome free text',
            'refname\tnon_coding\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t0\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1_n_A14T_N_ref has wild type, foo bar\tsome free text',
            'refname\tnon_coding\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t0\tSNP\tn\tA14T\t.\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1_n_A14T_N_ref has wild type, foo bar\tsome free text',
            'refname\tnon_coding\t19\t78\tcluster\t120\t100\t98.33\tctg_name\t279\t24.4\t0\tSNP\tn\tA14T\t.\t.\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tnoncoding1_n_A14T_N_ref has wild type, foo bar\tsome free text',
        ]

        expected = ['no', 'yes', 'no', 'yes', 'no']

        for i in range(len(lines)):
            data_dict = summary_cluster.SummaryCluster.line2dict(lines[i])
            cluster = summary_cluster.SummaryCluster()
            cluster.add_data_dict(data_dict)
            for assembled_summary in ['yes', 'fragmented', 'yes_nonunique']:
                self.assertEqual(expected[i], cluster._to_cluster_summary_number_has_nonsynonymous(assembled_summary))
            self.assertEqual('NA', cluster._to_cluster_summary_number_has_nonsynonymous('no'))


    def test_get_nonsynonymous_var(self):
        '''Test _get_nonsynonymous_var'''
        d = {
            'ref_name': 'ref',
            'var_type': '.',
            'known_var_change': '.',
            'has_known_var': '.',
            'known_var': '0',
            'ref_ctg_change': '.',
            'ref_ctg_effect': 'SNP',
            'var_seq_type': '.'
        }

        self.assertEqual(None, summary_cluster.SummaryCluster._get_nonsynonymous_var(d))

        d['var_type'] = 'p'
        d['has_known_var'] = '1'
        with self.assertRaises(summary_cluster.Error):
            summary_cluster.SummaryCluster._get_nonsynonymous_var(d)

        d['known_var_change'] = 'I42L'
        self.assertEqual('ref.I42L', summary_cluster.SummaryCluster._get_nonsynonymous_var(d))

        d['ref_ctg_change'] = 'P43Q'
        with self.assertRaises(summary_cluster.Error):
            summary_cluster.SummaryCluster._get_nonsynonymous_var(d)

        d['known_var_change'] = '.'
        self.assertEqual('ref.P43Q', summary_cluster.SummaryCluster._get_nonsynonymous_var(d))


    def test_column_summary_data(self):
        '''Test column_summary_data'''
        line1 = 'ref1\tnon_coding\t531\t78\tcluster1\t120\t100\t98.33\tctg_name\t279\t24.4\t1\tSNP\tn\tA14T\t1\tA14T\tSNP\t13\t13\tA\t84\t84\tT\t17\t.\t17\tfoo bar\tspam eggs'
        line2 = 'ref1\tnon_coding\t531\t78\tcluster1\t120\t95\t98.42\tctg_name\t279\t24.4\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\tsome free text'

        data_dict1 = summary_cluster.SummaryCluster.line2dict(line1)
        data_dict2 = summary_cluster.SummaryCluster.line2dict(line2)
        cluster = summary_cluster.SummaryCluster()
        cluster.add_data_dict(data_dict1)
        cluster.add_data_dict(data_dict2)
        expected = {
            'cluster1': 'yes',
            'cluster1.allele': 'ref1',
            'cluster1.any_var': 'yes',
            'cluster1.pct_id': '98.33',
        }
        got = cluster.column_summary_data()
        self.assertEqual(expected, got)

