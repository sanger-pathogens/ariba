import unittest
import filecmp
import os
from ariba import mic_plotter

modules_dir = os.path.dirname(os.path.abspath(mic_plotter.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


class TestMicPlotter(unittest.TestCase):
    def test_mic_string_to_float(self):
        '''Test _mic_string_to_float'''
        self.assertEqual(mic_plotter.MicPlotter._mic_string_to_float('42.42'), 42.42)
        self.assertEqual(mic_plotter.MicPlotter._mic_string_to_float('42'), 42.0)
        self.assertEqual(mic_plotter.MicPlotter._mic_string_to_float('>42'), 84.0)
        self.assertEqual(mic_plotter.MicPlotter._mic_string_to_float('> 42'), 84.0)
        self.assertEqual(mic_plotter.MicPlotter._mic_string_to_float('>=42'), 42.0)
        self.assertEqual(mic_plotter.MicPlotter._mic_string_to_float('>= 42'), 42.0)
        self.assertEqual(mic_plotter.MicPlotter._mic_string_to_float('<42'), 21.0)
        self.assertEqual(mic_plotter.MicPlotter._mic_string_to_float('< 42'), 21.0)
        self.assertEqual(mic_plotter.MicPlotter._mic_string_to_float('<=42'), 42.0)
        self.assertEqual(mic_plotter.MicPlotter._mic_string_to_float('<= 42'), 42.0)
        self.assertEqual(mic_plotter.MicPlotter._mic_string_to_float('   <=  42.0   '), 42.0)
        self.assertEqual(mic_plotter.MicPlotter._mic_string_to_float('na'), 'NA')
        self.assertEqual(mic_plotter.MicPlotter._mic_string_to_float('NA'), 'NA')
        self.assertEqual(mic_plotter.MicPlotter._mic_string_to_float('.'), 'NA')
        self.assertEqual(mic_plotter.MicPlotter._mic_string_to_float(' '), 'NA')
        self.assertEqual(mic_plotter.MicPlotter._mic_string_to_float(''), 'NA')


    def test_load_mic_file(self):
        '''Test _load_mic_file'''
        infile = os.path.join(data_dir, 'mic_plotter_load_mic_file.tsv')
        got = mic_plotter.MicPlotter._load_mic_file(infile)
        expected = {
            'sample1': {'antibio1': 0.25, 'antibio2': 0.004},
            'sample2': {'antibio1': 0.125, 'antibio2': 0.004},
            'sample3': {'antibio1': 0.125, 'antibio2': 0.004},
            'sample4': {'antibio1': 512.0, 'antibio2': 256.0},
            'sample5': {'antibio1': 512.0, 'antibio2': 256.0},
            'sample6': {'antibio1': 'NA', 'antibio2': 1.0},
        }

        self.assertEqual(got, expected)


    def test_load_summary_file(self):
        '''Test _load_summary_file'''
        infile = os.path.join(data_dir, 'mic_plotter_load_summary_file.tsv')
        got = mic_plotter.MicPlotter._load_summary_file(infile)
        expected = {
            'name1': {
                'cluster1': {'assembled': 'yes', 'match': 'yes', 'ref_seq': 'ref1', 'pct_id': 100.0, 'known_var': 'no', 'novel_var': 'no', 'group1.A42T': 'no', 'group1.A42T.%': 'NA'},
                'cluster2': {'assembled': 'yes', 'match': 'yes', 'ref_seq': 'ref2', 'pct_id': 99.0, 'known_var': 'yes', 'novel_var': 'no', 'group1.A42T': 'yes', 'group1.A42T.%': 95.42},
            },
            'name2': {
                'cluster1': {'assembled': 'yes', 'match': 'yes_nonunique', 'ref_seq': 'ref3', 'pct_id': 99.0, 'known_var': 'yes', 'novel_var': 'no', 'group1.A42T': 'yes', 'group1.A42T.%': 90.90},
                'cluster2': {'assembled': 'no', 'match': 'no', 'ref_seq': 'NA', 'pct_id': 'NA', 'known_var': 'NA', 'novel_var': 'NA', 'group1.A42T': 'NA', 'group1.A42T.%': 'NA'},
            },
        }
        self.maxDiff = None
        self.assertEqual(got, expected)


    def test_to_boxplot_tsv(self):
        '''Test _to_boxplot_tsv'''
        mic_data = {
            'name1': {'antibio1': 0.25, 'antibio2': 0.004},
            'name2': {'antibio1': 0.125, 'antibio2': 'NA'},
            'name3': {'antibio1': 'NA', 'antibio2': 0.002},
        }

        summary_data = {
            'name1': {
                'cluster1': {'assembled': 'yes', 'match': 'yes', 'ref_seq': 'ref1', 'pct_id': 100.0, 'known_var': 'no', 'novel_var': 'no', 'group1.A42T': 'no', 'group1.A42T.%': 'NA'},
                'cluster2': {'assembled': 'yes', 'match': 'yes', 'ref_seq': 'ref2', 'pct_id': 99.0, 'known_var': 'yes', 'novel_var': 'no', 'group2.A43T': 'yes', 'group2.A43T.%': 95.42},
                'cluster3': {'assembled': 'interrupted', 'match': 'no', 'ref_seq': 'ref3', 'pct_id': 99.0, 'known_var': 'no', 'novel_var': 'yes', 'A42T': 'no', 'A44T.%': 'NA'},
            },
            'name2': {
                'cluster1': {'assembled': 'yes', 'match': 'yes_nonunique', 'ref_seq': 'ref3', 'pct_id': 99.0, 'known_var': 'yes', 'novel_var': 'no', 'group1.A42T': 'yes', 'group1.A42T.%': 90.90},
                'cluster2': {'assembled': 'no', 'match': 'no', 'ref_seq': 'NA', 'pct_id': 'NA', 'known_var': 'NA', 'novel_var': 'NA', 'group2.A43T': 'NA', 'group2.A43T.%': 'NA'},
                'cluster3': {'assembled': 'no', 'match': 'no', 'ref_seq': 'NA', 'pct_id': 'NA', 'known_var': 'no', 'novel_var': 'no', 'A42T': 'no', 'A44T.%': 'NA'},
            },
            'name3': {
                'cluster1': {'assembled': 'yes', 'match': 'yes', 'ref_seq': 'ref_seq42', 'pct_id': 100.0, 'known_var': 'no', 'novel_var': 'no', 'group1.A42T': 'no', 'group1.A42T.%': 'NA'},
                'cluster2': {'assembled': 'no', 'match': 'no', 'ref_seq': 'NA', 'pct_id': 'NA', 'known_var': 'NA', 'novel_var': 'NA', 'group2.A43T': 'NA', 'group2.A43T.%': 'NA'},
                'cluster3': {'assembled': 'no', 'match': 'no', 'ref_seq': 'NA', 'pct_id': 'NA', 'known_var': 'no', 'novel_var': 'no', 'A42T': 'no', 'A44T.%': 'NA'},
            },
        } 

        tmp_tsv = 'tmp.mic_plotter_test.to_boxplot.tsv'
        expected_mutations = [
            {'cluster1.group1.A42T', 'cluster2.group2.A43T', 'cluster3.interrupted'},
            {'without_mutation', 'cluster2.group2.A43T', 'cluster3.interrupted'},
        ]

        expected_combs = [
            {('cluster2.group2.A43T', 'cluster3.interrupted'), ('cluster1.group1.A42T',)},
            {('without_mutation',), ('cluster2.group2.A43T', 'cluster3.interrupted')}
        ]

        for i in [1, 2]:
            got_muts, got_combs = mic_plotter.MicPlotter._to_boxplot_tsv(summary_data, mic_data, 'antibio' + str(i), tmp_tsv)
            expected_tsv = os.path.join(data_dir, 'mic_plotter_to_boxplot_tsv.antibio' + str(i) + '.tsv')
            self.assertTrue(filecmp.cmp(tmp_tsv, expected_tsv, shallow=False))
            self.assertEqual(got_muts, expected_mutations[i-1])
            self.assertEqual(got_combs, expected_combs[i-1])
            os.unlink(tmp_tsv)


    def test_to_dots_tsv(self):
        '''test _to_dots_tsv'''
        all_mutations = {'m1', 'm2', 'm3'}
        combinations = {
            ('m1',),
            ('m1', 'm3'),
            ('m2', 'm3'),
        }

        tmp_tsv = 'tmp.test.mic_plotter_to_dots.tsv'
        expected_tsv = os.path.join(data_dir, 'mic_plotter_to_dots.tsv')
        mic_plotter.MicPlotter._to_dots_tsv(all_mutations, combinations, tmp_tsv)
        self.assertTrue(filecmp.cmp(tmp_tsv, expected_tsv, shallow=False))
        os.unlink(tmp_tsv)

        all_mutations.update({'without_mutation', 'z1'})
        combinations.add(('without_mutation',))
        combinations.add(('m1', 'z1'))
        expected_tsv = os.path.join(data_dir, 'mic_plotter_to_dots_without_mutation.tsv')
        mic_plotter.MicPlotter._to_dots_tsv(all_mutations, combinations, tmp_tsv)
        self.assertTrue(filecmp.cmp(tmp_tsv, expected_tsv, shallow=False))
        os.unlink(tmp_tsv)

