import unittest
import os
import filecmp
from ariba import flag, report_filter, report

modules_dir = os.path.dirname(os.path.abspath(report_filter.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


class TestReportFilter(unittest.TestCase):
    def test_init_good_file(self):
        '''test __init__ on good input file'''
        infile = os.path.join(data_dir, 'report_filter_test_init_good.tsv')
        rf = report_filter.ReportFilter(infile=infile)
        line1 = '\t'.join(['ariba_cluster1', 'cluster1', '0', '0', '27', '10000', 'cluster1', '1000', '999', '99.42', 'cluster1.scaffold.1', '1300', '10.5', '1', 'SNP', 'n', 'C42T', '0', '.', '.', '42', '42', 'C', '142', '142', 'C', '500', 'C', '500', 'a:n:C42T:id1:foo', 'free_text'])
        line2 = '\t'.join(['ariba_cluster1', 'cluster1', '0', '0', '27', '10000', 'cluster1', '1000', '999', '99.42', 'cluster1.scaffold.1', '1300', '10.5', '1', 'SNP', 'n', 'A51G', '0', '.', '.', '51', '51', 'C', '151', '151', 'C', '542', 'C', '542', 'a:n:A51G:id2:bar', 'free_text2'])
        line3 = '\t'.join(['ariba_cluster1', 'cluster1', '0', '0', '27', '10000', 'cluster1', '1000', '999', '99.42', 'cluster1.scaffold.2', '1300', '12.4', '1', 'SNP', 'n', 'A51G', '0', '.', '.', '51', '51', 'C', '151', '151', 'C', '542', 'C', '542', 'a:n:A51G:id3:spam', 'free_text3'])
        line4 = '\t'.join(['ariba_cluster2', 'cluster2', '1', '0', '179', '20000', 'cluster2', '1042', '1042', '42.42', 'cluster2.scaffold.1', '1442', '20.2', '1', 'SNP', 'p', 'I42L', '1', 'I42L', 'NONSYN', '112', '112', 'C', '442', '442', 'T', '300', 'T', '290', 'a:v:I42L:id4:eggs', 'free_text3'])

        expected = {
            'cluster1': {
                'cluster1.scaffold.1': [report_filter.ReportFilter._report_line_to_dict(line1), report_filter.ReportFilter._report_line_to_dict(line2)],
                'cluster1.scaffold.2': [report_filter.ReportFilter._report_line_to_dict(line3)],
            },
            'cluster2': {
                'cluster2.scaffold.1': [report_filter.ReportFilter._report_line_to_dict(line4)]
            }
        }

        self.assertEqual(expected, rf.report)


    def test_init_bad_file(self):
        '''test __init__ on bad input file'''
        infile = os.path.join(data_dir, 'report_filter_test_init_bad.tsv')
        with self.assertRaises(report_filter.Error):
            report_filter.ReportFilter(infile=infile)


    def test_report_line_to_dict(self):
        line = 'ariba_cluster1\tcluster1\t0\t0\t27\t10000\tcluster1\t1000\t999\t99.42\tcluster1.scaffold.1\t999\t23.2\t1\tSNP\tn\tC42T\t0\t.\t.\t42\t42\tC\t142\t142\tC\t500\t.\t500\ta:n:C42T:id1:foo\tfree text'
        expected = {
            'ariba_ref_name':     'ariba_cluster1',
            'ref_name':           'cluster1',
            'gene':               '0',
            'var_only':           '0',
            'flag':               flag.Flag(27),
            'reads':              10000,
            'cluster':        'cluster1',
            'ref_len':            1000,
            'ref_base_assembled': 999,
            'pc_ident':           99.42,
            'ctg':                'cluster1.scaffold.1',
            'ctg_len':            999,
            'ctg_cov':            23.2,
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
            'smtls_nts':       '.',
            'smtls_nts_depth':    '500',
            'var_description':    'a:n:C42T:id1:foo',
            'free_text':          'free text',
        }

        self.assertEqual(expected, report_filter.ReportFilter._report_line_to_dict(line))

        bad_line = '\t'.join(line.split('\t')[:-1])
        self.assertEqual(None, report_filter.ReportFilter._report_line_to_dict(bad_line))


    def test_dict_to_report_line(self):
        '''Test _dict_to_report_line'''
        report_dict = {
            'ariba_ref_name':     'ariba_cluster1',
            'ref_name':           'cluster1',
            'gene':               '0',
            'var_only':           '0',
            'flag':               flag.Flag(27),
            'reads':              10000,
            'cluster':        'cluster1',
            'ref_len':            1000,
            'ref_base_assembled': 999,
            'pc_ident':           99.42,
            'ctg':                'cluster1.scaffold.1',
            'ctg_len':            1300,
            'ctg_cov':            42.4,
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
            'smtls_nts':       '.',
            'smtls_nts_depth':    '500',
            'var_description':    'a:n:C42T:id1:foo',
            'free_text':          'free text',
        }

        expected = 'ariba_cluster1\tcluster1\t0\t0\t27\t10000\tcluster1\t1000\t999\t99.42\tcluster1.scaffold.1\t1300\t42.4\t1\tSNP\tn\tC42T\t0\t.\t.\t42\t42\tC\t142\t142\tC\t500\t.\t500\ta:n:C42T:id1:foo\tfree text'
        self.assertEqual(expected, report_filter.ReportFilter._dict_to_report_line(report_dict))


    def test_load_report(self):
        good_infile = os.path.join(data_dir, 'report_filter_test_load_report_good.tsv')
        bad_infile = os.path.join(data_dir, 'report_filter_test_load_report_bad.tsv')

        line1 = '\t'.join(['ariba_cluster1', 'cluster1', '0', '0', '27', '10000', 'cluster1', '1000', '999', '99.42', 'cluster1.scaffold.1', '1300', '12.2', '1', 'SNP', 'n', 'C42T', '0', '.', '.', '42', '42', 'C', '142', '142', 'C', '500', 'C', '500', 'a:n:C42T:id1:foo', 'free_text'])
        line2 = '\t'.join(['ariba_cluster1', 'cluster1', '0', '0', '27', '10000', 'cluster1', '1000', '999', '99.42', 'cluster1.scaffold.1', '1300', '12.2', '1', 'SNP', 'n', 'A51G', '0', '.', '.', '51', '51', 'C', '151', '151', 'C', '542', 'C', '542', 'a:n:A51G:id2:bar', 'free_text2'])
        line3 = '\t'.join(['ariba_cluster1', 'cluster1', '0', '0', '27', '10000', 'cluster1', '1000', '999', '99.42', 'cluster1.scaffold.2', '1300', '22.2', '1', 'SNP', 'n', 'A51G', '0', '.', '.', '51', '51', 'C', '151', '151', 'C', '542', 'C', '542', 'a:n:A51G:id3:spam', 'free_text3'])
        line4 = '\t'.join(['ariba_cluster2', 'cluster2', '1', '1', '179', '20000', 'cluster2', '1042', '1042', '42.42', 'cluster2.scaffold.1', '1442', '33.3', '1', 'SNP', 'p', 'I42L', '1', 'I42L', 'NONSYN', '112', '112', 'C', '442', '442', 'T', '300', 'T', '290', 'a:v:I42L:id4:eggs', 'free_text3'])

        expected = {
            'cluster1': {
                'cluster1.scaffold.1': [report_filter.ReportFilter._report_line_to_dict(line1), report_filter.ReportFilter._report_line_to_dict(line2)],
                'cluster1.scaffold.2': [report_filter.ReportFilter._report_line_to_dict(line3)],
            },
            'cluster2': {
                'cluster2.scaffold.1': [report_filter.ReportFilter._report_line_to_dict(line4)]
            }
        }

        got = report_filter.ReportFilter._load_report(good_infile)
        self.maxDiff = None
        self.assertEqual(expected, got)
        with self.assertRaises(report_filter.Error):
            report_filter.ReportFilter._load_report(bad_infile)


    def test_report_dict_passes_non_essential_filters_known_vars(self):
        '''Test _report_dict_passes_non_essential_filters with known vars'''
        tests = [
            ('.', '.', True, True),
            ('.', '.', False, True),
            ('0', '0', True, True),
            ('0', '0', False, True),
            ('1', '0', True, False),
            ('1', '1', True, True),
            ('1', '0', False, True),
            ('1', '1', False, True),
        ]

        for known_var, has_known_var, ignore_not_has_known_variant, expected in tests:
            d = {'known_var': known_var, 'has_known_var': has_known_var}
            rf = report_filter.ReportFilter(ignore_not_has_known_variant=ignore_not_has_known_variant)
            self.assertEqual(expected, rf._report_dict_passes_non_essential_filters(d))


    def test_report_dict_passes_non_essential_filters_synonymous(self):
        '''Test _report_dict_passes_non_essential_filters with synonymous AA changes'''
        tests = [
             ('.', True, True),
             ('.', False, True),
             ('SNP', True, True),
             ('SNP', False, True),
             ('SYN', True, False),
             ('SYN', False, True),
        ]

        for var, remove_synonymous_snps, expected in tests:
            d = {'known_var': '1', 'ref_ctg_effect': var, 'has_known_var': '1'}
            rf = report_filter.ReportFilter(remove_synonymous_snps=remove_synonymous_snps)
            self.assertEqual(expected, rf._report_dict_passes_non_essential_filters(d))


    def test_report_dict_passes_essential_filters(self):
        '''Test _report_dict_passes_essential_filters'''
        line1 = 'ariba_cluster1\tcluster1\t0\t0\t27\t10000\tcluster1\t1000\t999\t98.42\tcluster1.scaffold.1\t400\t12.2\t1\tSNP\tn\tC42T\t0\t.\t.\t42\t42\tC\t142\t142\tC\t500\t.\t500\tDescription_of_variant C42T\tfree text'
        line2 = 'ariba_cluster1\tcluster1\t0\t0\t27\t10000\tcluster1\t1000\t0\t98.42\tcluster1.scaffold.1\t400\t12.2\t1\tSNP\tn\tC42T\t0\t.\t.\t42\t42\tC\t142\t142\tC\t500\t.\t500\tDescription_of_variant C42T\tfree text'
        line3 = 'ariba_cluster1\tcluster1\t0\t0\t27\t10000\tcluster1\t1000\t999\t78.42\tcluster1.scaffold.1\t400\t12.2\t1\tSNP\tn\tC42T\t0\t.\t.\t42\t42\tC\t142\t142\tC\t500\t.\t500\tDescription_of_variant C42T\tfree text'
        tests = [
            (report_filter.ReportFilter._report_line_to_dict(line1), True),
            (report_filter.ReportFilter._report_line_to_dict(line2), False),
            (report_filter.ReportFilter._report_line_to_dict(line3), False),
        ]

        for test_dict, expected in tests:
            rf = report_filter.ReportFilter()
            self.assertEqual(expected,  rf._report_dict_passes_essential_filters(test_dict))


    def test_flag_passes_filter(self):
        '''Test _flag_passes_filter'''
        rf = report_filter.ReportFilter()
        exclude_flags = ['assembly_fail', 'ref_seq_choose_fail']
        f = flag.Flag()
        self.assertTrue(rf._flag_passes_filter(f, exclude_flags))
        f.add('assembled')
        self.assertTrue(rf._flag_passes_filter(f, exclude_flags))
        f = flag.Flag()
        f.add('assembly_fail')
        self.assertFalse(rf._flag_passes_filter(f, exclude_flags))
        f = flag.Flag()
        f.add('ref_seq_choose_fail')
        self.assertFalse(rf._flag_passes_filter(f, exclude_flags))


    def test_filter_list_of_dicts_all_fail(self):
        '''Test _filter_list_of_dicts where all fail'''
        rf = report_filter.ReportFilter()
        line1 = 'ariba_cluster1\tcluster1\t0\t0\t27\t10000\tcluster1\t1000\t999\t88.42\tcluster1.scaffold.1\t400\t12.2\t1\tSNP\tn\tC42T\t0\t.\t.\t42\t42\tC\t142\t142\tC\t500\t.\t500\tDescription_of_variant C42T\tfree text'
        line2 = 'ariba_cluster1\tcluster1\t0\t0\t27\t10000\tcluster1\t1000\t999\t78.42\tcluster1.scaffold.1\t400\t12.2\t1\tSNP\tn\tC42T\t0\t.\t.\t42\t42\tC\t142\t142\tC\t500\t.\t500\tDescription_of_variant C42T\tfree text'
        dict1 = report_filter.ReportFilter._report_line_to_dict(line1)
        dict2 = report_filter.ReportFilter._report_line_to_dict(line2)
        got = rf._filter_list_of_dicts([dict1, dict2])
        self.assertEqual([], got)


    def test_filter_list_of_dicts_with_essential(self):
        '''Test _filter_list_of_dicts with an essential line but all others fail'''
        rf = report_filter.ReportFilter(ignore_not_has_known_variant=True)
        line1 = 'ariba_cluster1\tcluster1\t0\t0\t27\t10000\tcluster1\t1000\t999\t98.42\tcluster1.scaffold.1\t400\t12.2\t1\tSNP\tn\tC42T\t0\t.\t.\t42\t42\tC\t142\t142\tC\t500\t.\t500\tDescription_of_variant C42T\tfree text'
        line2 = 'ariba_cluster1\tcluster1\t0\t0\t27\t10000\tcluster1\t1000\t999\t78.42\tcluster1.scaffold.1\t400\t12.2\t1\tSNP\tn\tC42T\t0\t.\t.\t42\t42\tC\t142\t142\tC\t500\t.\t500\tDescription_of_variant C42T\tfree text'
        dict1 = report_filter.ReportFilter._report_line_to_dict(line1)
        dict2 = report_filter.ReportFilter._report_line_to_dict(line2)
        expected_line = 'ariba_cluster1\tcluster1\t0\t0\t27\t10000\tcluster1\t1000\t999\t98.42\tcluster1.scaffold.1\t400\t12.2\t' + '\t'.join(['.'] * 17) + '\tfree text'
        expected = [report_filter.ReportFilter._report_line_to_dict(expected_line)]
        assert expected != [None]
        got = rf._filter_list_of_dicts([dict1, dict2])
        self.assertEqual(expected, got)


    def test_filter_list_of_dicts_with_pass(self):
        '''Test _filter_list_of_dicts with a line that passes'''
        rf = report_filter.ReportFilter(ignore_not_has_known_variant=True)
        line1 = 'ariba_cluster1\tcluster1\t0\t0\t27\t10000\tcluster1\t1000\t999\t98.42\tcluster1.scaffold.1\t500\t12.1\t1\tSNP\tn\tC42T\t0\t.\t.\t42\t42\tC\t142\t142\tC\t500\t.\t500\tDescription_of_variant C42T\tfree text'
        line2 = 'ariba_cluster1\tcluster1\t0\t0\t27\t10000\tcluster1\t1000\t999\t98.42\tcluster1.scaffold.1\t500\t12.1\t1\tSNP\tn\tC46T\t1\t.\t.\t42\t42\tC\t142\t142\tC\t500\t.\t500\tDescription_of_variant C46T\tfree text'
        line3 = 'ariba_cluster1\tcluster1\t0\t0\t27\t10000\tcluster1\t1000\t999\t78.42\tcluster1.scaffold.1\t500\t12.1\t1\tSNP\tn\tC42T\t0\t.\t.\t42\t42\tC\t142\t142\tC\t500\t.\t500\tDescription_of_variant C42T\tfree text'
        dict1 = report_filter.ReportFilter._report_line_to_dict(line1)
        dict2 = report_filter.ReportFilter._report_line_to_dict(line2)
        dict3 = report_filter.ReportFilter._report_line_to_dict(line3)
        got = rf._filter_list_of_dicts([dict1, dict2, dict3])
        self.assertEqual([dict2], got)


    def test_remove_all_after_first_frameshift(self):
        '''Test _remove_all_after_first_frameshift'''
        self.assertEqual([], report_filter.ReportFilter._remove_all_after_first_frameshift([]))
        line1 = 'ariba_cluster1\tcluster1\t1\t0\t528\t1874\tcluster1\t1188\t1097\t92.43\tcluster1.scaffold.1\t2218\t42.42\t0\t.\tp\t.\t0\tE89G\tNONSYN\t65\t265\tA;A\t766\t766\tG;C\t88;90\t.;.\t87;90\t.\t.'
        line2 = 'ariba_cluster1\tcluster1\t1\t0\t528\t1874\tcluster1\t1188\t1097\t92.43\tcluster1.scaffold.1\t2218\t42.42\t0\t.\tp\t.\t0\tQ37fs\tFSHIFT\t109\t109\tA\t634\t634\t.\t67\t.\t67\t.\t.'
        line3 = 'ariba_cluster1\tcluster1\t1\t0\t528\t1874\tcluster1\t1188\t1097\t92.43\tcluster1.scaffold.1\t2218\t42.42\t0\t.\tp\t.\t0\tE89G\tNONSYN\t265\t265\tA;A\t766\t766\tG;C\t88;90\t.;.\t87;90\t.\t.'
        dict1 = report_filter.ReportFilter._report_line_to_dict(line1)
        dict2 = report_filter.ReportFilter._report_line_to_dict(line2)
        dict3 = report_filter.ReportFilter._report_line_to_dict(line3)
        self.assertEqual([dict1], report_filter.ReportFilter._remove_all_after_first_frameshift([dict1]))
        self.assertEqual([dict1, dict2], report_filter.ReportFilter._remove_all_after_first_frameshift([dict1, dict2]))
        self.assertEqual([dict2], report_filter.ReportFilter._remove_all_after_first_frameshift([dict2, dict3]))
        self.assertEqual([dict1, dict2], report_filter.ReportFilter._remove_all_after_first_frameshift([dict1, dict2, dict3]))


    def test_remove_all_after_first_frameshift_complete_gene(self):
        '''Test _remove_all_after_first_frameshift when gene is complete'''
        self.assertEqual([], report_filter.ReportFilter._remove_all_after_first_frameshift([]))
        line1 = 'ariba_cluster1\tcluster1\t1\t0\t538\t1874\tcluster1\t1188\t1097\t92.43\tcluster1.scaffold.1\t2218\t42.42\t0\t.\tp\t.\t0\tE89G\tNONSYN\t65\t265\tA;A\t766\t766\tG;C\t88;90\t.;.\t87;90\t.\t.'
        line2 = 'ariba_cluster1\tcluster1\t1\t0\t538\t1874\tcluster1\t1188\t1097\t92.43\tcluster1.scaffold.1\t2218\t42.42\t0\t.\tp\t.\t0\tQ37fs\tFSHIFT\t109\t109\tA\t634\t634\t.\t67\t.\t67\t.\t.'
        line3 = 'ariba_cluster1\tcluster1\t1\t0\t538\t1874\tcluster1\t1188\t1097\t92.43\tcluster1.scaffold.1\t2218\t42.42\t0\t.\tp\t.\t0\tE89G\tNONSYN\t265\t265\tA;A\t766\t766\tG;C\t88;90\t.;.\t87;90\t.\t.'
        dict1 = report_filter.ReportFilter._report_line_to_dict(line1)
        dict2 = report_filter.ReportFilter._report_line_to_dict(line2)
        dict3 = report_filter.ReportFilter._report_line_to_dict(line3)
        self.assertEqual([dict1], report_filter.ReportFilter._remove_all_after_first_frameshift([dict1]))
        self.assertEqual([dict1, dict2], report_filter.ReportFilter._remove_all_after_first_frameshift([dict1, dict2]))
        self.assertEqual([dict2, dict3], report_filter.ReportFilter._remove_all_after_first_frameshift([dict2, dict3]))
        self.assertEqual([dict1, dict2, dict3], report_filter.ReportFilter._remove_all_after_first_frameshift([dict1, dict2, dict3]))


    def test_filter_dicts(self):
        '''Test _filter_dicts'''
        rf = report_filter.ReportFilter(min_ref_base_assembled=10, ignore_not_has_known_variant=True)
        ref_2_dict = {x: '.' for x in report.columns}
        ref_2_dict['pc_ident'] = 91.0
        ref_2_dict['ref_base_assembled'] = 10
        ref_2_dict['has_known_var'] = '0'
        ref_2_dict['flag'] = flag.Flag(27)
        ref_2_dict['var_type'] = '.'

        rf.report = {
            'ref1': {
                'ref1.scaff1': [
                    {'flag': flag.Flag(27), 'pc_ident': 91.0, 'ref_base_assembled': 9, 'known_var': '1', 'has_known_var': '1', 'var_type': 'SNP'},
                    {'flag': flag.Flag(27), 'pc_ident': 91.5, 'ref_base_assembled': 11, 'known_var': '1', 'has_known_var': '1', 'var_type': 'HET'},
                    {'flag': flag.Flag(27), 'pc_ident': 89.0, 'ref_base_assembled': 10, 'known_var': '1', 'has_known_var': '1', 'var_type': 'SNP'},
                    {'flag': flag.Flag(27), 'pc_ident': 90.0, 'ref_base_assembled': 11, 'known_var': '1', 'has_known_var': '0', 'var_type': 'SNP'},
                    {'flag': flag.Flag(27), 'pc_ident': 90.0, 'ref_base_assembled': 11, 'known_var': '1', 'has_known_var': '1', 'var_type': 'SNP'},
                ]
            },
            'ref2': {
                'ref2.scaff1': [
                    ref_2_dict
                ]
            },
            'ref3': {
                'ref3.scaff1': [
                    {'flag': flag.Flag(27), 'pc_ident': 84.0, 'ref_base_assembled': 10, 'known_var': '1', 'has_known_var': '0', 'var_type': 'SNP'},
                ]
            },
            'ref4': {
                'ref4.scaff1': [
                    {'flag': flag.Flag(64), 'pc_ident': '.', 'ref_base_assembled': '.', 'known_var': '.', 'has_known_var': '.', 'var_type': '.'},
                ]
            }
        }

        expected = {
            'ref1': {
                'ref1.scaff1': [
                    {'flag': flag.Flag(27), 'pc_ident': 91.5, 'ref_base_assembled': 11, 'known_var': '1', 'has_known_var': '1', 'var_type': 'HET'},
                    {'flag': flag.Flag(27), 'pc_ident': 90.0, 'ref_base_assembled': 11, 'known_var': '1', 'has_known_var': '1', 'var_type': 'SNP'},
                ]
            },
            'ref2': {
                'ref2.scaff1': [ref_2_dict]
            }
        }

        rf._filter_dicts()
        self.assertEqual(expected, rf.report)


    def test_write_report_tsv(self):
        '''Test write_report_tsv'''
        infile = os.path.join(data_dir, 'report_filter_test_write_report.tsv')
        tmpfile = 'tmp.test.report_filter.write_report.tsv'
        rf = report_filter.ReportFilter(infile=infile)
        rf._write_report_tsv(tmpfile)
        self.assertTrue(filecmp.cmp(tmpfile, infile, shallow=False))
        os.unlink(tmpfile)


    def test_run(self):
        '''Test run'''
        infile = os.path.join(data_dir, 'report_filter_test_run.in.tsv')
        expected_file = os.path.join(data_dir, 'report_filter_test_run.expected.tsv')
        tmpfile = 'tmp.test.report_filter.run.out.tsv'
        rf = report_filter.ReportFilter(infile=infile)
        rf.run(tmpfile)
        self.assertTrue(filecmp.cmp(expected_file, tmpfile, shallow=False))
        os.unlink(tmpfile)

