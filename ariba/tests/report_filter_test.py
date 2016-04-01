import unittest
import os
import filecmp
import pyfastaq
from ariba import report_filter

modules_dir = os.path.dirname(os.path.abspath(report_filter.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


class TestReportFilter(unittest.TestCase):
    def test_init_good_file(self):
        '''test __init__ on good input file'''
        infile = os.path.join(data_dir, 'report_filter_test_init_good.tsv')
        rf = report_filter.ReportFilter(infile)
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
            rf = report_filter.ReportFilter(infile)


    def test_report_line_to_dict(self):
        line = 'cluster1\tnon_coding\t27\t10000\tcluster1\t1000\t999\t99.42\tcluster1.scaffold.1\t1300\t1\tSNP\tn\tC42T\t0\t.\t.\t42\t42\tC\t142\t142\tC\t500\t.\t500\tDescription_of_variant C42T\tfree text'
        expected = {
            'ref_name':          'cluster1',
            'ref_type':          'non_coding',
            'flag':              '27',
            'reads':             '10000',
            'cluster_rep':       'cluster1',
            'ref_len':           '1000',
            'ref_base_assembled':'999',
            'pc_ident':          '99.42',
            'ctg':               'cluster1.scaffold.1',
            'ctg_len':           '1300',
            'known_var':         '1',
            'var_type':          'SNP',
            'var_seq_type':      'n',
            'known_var_change':  'C42T',
            'has_known_var':     '0',
            'ref_ctg_change':    '.',
            'ref_ctg_effect':    '.',
            'ref_start':         '42',
            'ref_end':           '42',
            'ref_nt':            'C',
            'ctg_start':         '142',
            'ctg_end':           '142',
            'ctg_nt':            'C',
            'smtls_total_depth': '500',
            'smtls_alt_nt':      '.',
            'smtls_alt_depth':   '500',
            'var_description':   'Description_of_variant C42T',
            'free_text':         'free text',
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
