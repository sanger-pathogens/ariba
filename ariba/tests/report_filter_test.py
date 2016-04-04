import unittest
import os
import filecmp
import pyfastaq
from ariba import flag, report_filter

modules_dir = os.path.dirname(os.path.abspath(report_filter.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


class TestReportFilter(unittest.TestCase):
    def test_init_good_file(self):
        '''test __init__ on good input file'''
        infile = os.path.join(data_dir, 'report_filter_test_init_good.tsv')
        rf = report_filter.ReportFilter(infile=infile)
        line1 = '\t'.join(['cluster1', 'non_coding', '27', '10000', 'cluster1', '1000', '999', '99.42', 'cluster1.scaffold.1', '1300', '1', 'SNP', 'n', 'C42T', '0', '.', '.', '42', '42', 'C', '142', '142', 'C', '500', '.', '500', 'Description_of_variant.C42T', 'free_text'])
        line2 = '\t'.join(['cluster1', 'non_coding', '27', '10000', 'cluster1', '1000', '999', '99.42', 'cluster1.scaffold.1', '1300', '1', 'SNP', 'n', 'A51G', '0', '.', '.', '51', '51', 'C', '151', '151', 'C', '542', '.', '542', 'Description_of_variant.A51G', 'free_text2'])
        line3 = '\t'.join(['cluster2', 'variants_only', '179', '20000', 'cluster2', '1042', '1042', '42.42', 'cluster2.scaffold.1', '1442', '1', 'SNP', 'p', 'I42L', '1', 'I42L', 'NONSYN', '112', '112', 'C', '442', '442', 'T', '300', '.', '290', 'Description_of_variant.I42L', 'free_text3'])

        expected = {
            'cluster1': [
                report_filter.ReportFilter._report_line_to_dict(line1),
                report_filter.ReportFilter._report_line_to_dict(line2),
            ],
            'cluster2': [
                report_filter.ReportFilter._report_line_to_dict(line3),
            ]
        }

        self.assertEqual(expected, rf.report)


    def test_init_bad_file(self):
        '''test __init__ on bad input file'''
        infile = os.path.join(data_dir, 'report_filter_test_init_bad.tsv')
        with self.assertRaises(report_filter.Error):
            rf = report_filter.ReportFilter(infile=infile)


    def test_report_line_to_dict(self):
        line = 'cluster1\tnon_coding\t27\t10000\tcluster1\t1000\t999\t99.42\tcluster1.scaffold.1\t1300\t1\tSNP\tn\tC42T\t0\t.\t.\t42\t42\tC\t142\t142\tC\t500\t.\t500\tDescription_of_variant C42T\tfree text'
        expected = {
            'ref_name':           'cluster1',
            'ref_type':           'non_coding',
            'flag':               flag.Flag(27),
            'reads':              10000,
            'cluster_rep':        'cluster1',
            'ref_len':            1000,
            'ref_base_assembled': 999,
            'pc_ident':           99.42,
            'ctg':                'cluster1.scaffold.1',
            'ctg_len':            1300,
            'known_var':          '1',
            'var_type':           'SNP',
            'var_seq_type':       'n',
            'known_var_change':   'C42T',
            'has_known_var':      '0',
            'ref_ctg_change':     '.',
            'ref_ctg_effect':     '.',
            'ref_start':          42,
            'ref_end':            42,
            'ref_nt':             'C',
            'ctg_start':          142,
            'ctg_end':            142,
            'ctg_nt':             'C',
            'smtls_total_depth':  '500',
            'smtls_alt_nt':       '.',
            'smtls_alt_depth':    '500',
            'var_description':    'Description_of_variant C42T',
            'free_text':          'free text',
        }
        self.assertEqual(expected, report_filter.ReportFilter._report_line_to_dict(line))

        bad_line = '\t'.join(line.split('\t')[:-1])
        self.assertEqual(None, report_filter.ReportFilter._report_line_to_dict(bad_line))


    def test_load_report(self):
        good_infile = os.path.join(data_dir, 'report_filter_test_load_report_good.tsv')
        bad_infile = os.path.join(data_dir, 'report_filter_test_load_report_bad.tsv')

        line1 = '\t'.join(['cluster1', 'non_coding', '27', '10000', 'cluster1', '1000', '999', '99.42', 'cluster1.scaffold.1', '1300', '1', 'SNP', 'n', 'C42T', '0', '.', '.', '42', '42', 'C', '142', '142', 'C', '500', '.', '500', 'Description_of_variant.C42T', 'free_text'])
        line2 = '\t'.join(['cluster1', 'non_coding', '27', '10000', 'cluster1', '1000', '999', '99.42', 'cluster1.scaffold.1', '1300', '1', 'SNP', 'n', 'A51G', '0', '.', '.', '51', '51', 'C', '151', '151', 'C', '542', '.', '542', 'Description_of_variant.A51G', 'free_text2'])
        line3 = '\t'.join(['cluster2', 'variants_only', '179', '20000', 'cluster2', '1042', '1042', '42.42', 'cluster2.scaffold.1', '1442', '1', 'SNP', 'p', 'I42L', '1', 'I42L', 'NONSYN', '112', '112', 'C', '442', '442', 'T', '300', '.', '290', 'Description_of_variant.I42L', 'free_text3'])

        expected = {
            'cluster1': [
                report_filter.ReportFilter._report_line_to_dict(line1),
                report_filter.ReportFilter._report_line_to_dict(line2),
            ],
            'cluster2': [
                report_filter.ReportFilter._report_line_to_dict(line3),
            ]
        }

        got = report_filter.ReportFilter._load_report(good_infile)
        self.assertEqual(expected, got)
        with self.assertRaises(report_filter.Error):
            report_filter.ReportFilter._load_report(bad_infile)


    def test_report_dict_passes_known_variant_filter(self):
        tests = [
            ('.', True, False),
            ('0', True, False),
            ('1', True, True),
            ('.', False, True),
            ('0', False, True),
            ('1', False, True),
        ]

        for has_known_var, ignore_not_has_known_variant, expected in tests:
            d = {'has_known_var': has_known_var}
            self.assertEqual(expected, report_filter.ReportFilter._report_dict_passes_known_variant_filter(ignore_not_has_known_variant, d))


    def test_report_dict_passes_filters_known_variants(self):
        '''Test _report_dict_passes_filters with different known variants options'''
        rf = report_filter.ReportFilter()

        test_dict = {
            'pc_ident': 95.0,
            'ref_base_assembled': 10,
            'has_known_var': '1'
        }

        self.assertTrue(rf._report_dict_passes_filters(test_dict))
        tests = [
            ('.', True, False),
            ('0', True, False),
            ('1', True, True),
            ('.', False, True),
            ('0', False, True),
            ('1', False, True),
        ]

        for has_known_var, ignore_not_has_known_variant, expected in tests:
            test_dict['has_known_var'] = has_known_var
            rf = report_filter.ReportFilter(ignore_not_has_known_variant=ignore_not_has_known_variant)
            self.assertEqual(expected, rf._report_dict_passes_filters(test_dict))


    def test_report_dict_passes_filters_pc_ident(self):
        '''Test _report_dict_passes_filters with different percent identities'''
        test_dict = {
            'pc_ident': 95.0,
            'ref_base_assembled': 10,
            'has_known_var': '1'
        }

        tests = [
            (94.9, True),
            (95.0, True),
            (95.1, False)
        ]

        for cutoff, expected in tests:
            rf = report_filter.ReportFilter(min_pc_ident=cutoff)
            self.assertEqual(expected, rf._report_dict_passes_filters(test_dict))


    def test_report_dict_passes_filters_ref_base_assembled(self):
        '''Test _report_dict_passes_filters with different ref bases assembled'''
        test_dict = {
            'pc_ident': 95.0,
            'ref_base_assembled': 10,
            'has_known_var': '1'
        }

        tests = [
            (9, True),
            (10, True),
            (11, False)
        ]

        for cutoff, expected in tests:
            rf = report_filter.ReportFilter(min_ref_base_assembled=cutoff)
            self.assertEqual(expected, rf._report_dict_passes_filters(test_dict))

