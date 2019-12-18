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


    def test_get_colours(self):
        '''test _get_colours'''
        col1 = (0.0, 0.0, 0.5, 1.0)
        col2 = (0.0, 0.0, 0.517825311942959, 1.0)

        tests = [
            (1, 1, 'jet', ["black"]),
            (2, 1, 'jet', ["black", "black"]),
            (3, 1, 'jet', ["black", "black", "black"]),
            (2, 2, 'jet', [col1, col2]),
            (3, 2, 'jet', [col1, col2, col1]),
            (4, 2, 'jet', [col1, col2, col1, col2]),
            (3, 0, 'jet', [(0.0, 0.0, 0.5, 1.0), (0.49019607843137247, 1.0, 0.47754585705249841, 1.0), (0.5, 0.0, 0.0, 1.0)])
        ]

        for total_length, number_of_colours, colormap, expected in tests:
            self.assertEqual(expected, mic_plotter.MicPlotter._get_colours(total_length, number_of_colours, colormap))


    def test_get_top_plot_data(self):
        '''Test _get_top_plot_data'''
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
                'cluster4': {'assembled': 'yes', 'match': 'yes', 'ref_seq': 'ref4', 'pct_id': 100.0, 'known_var': 'no', 'novel_var': 'no', 'group4.A44T': 'no', 'group4.A44T.%': 'NA'},
            },
            'name2': {
                'cluster1': {'assembled': 'yes', 'match': 'yes_nonunique', 'ref_seq': 'ref3', 'pct_id': 99.0, 'known_var': 'yes', 'novel_var': 'no', 'group1.A42T': 'yes', 'group1.A42T.%': 90.90},
                'cluster2': {'assembled': 'no', 'match': 'no', 'ref_seq': 'NA', 'pct_id': 'NA', 'known_var': 'NA', 'novel_var': 'NA', 'group2.A43T': 'NA', 'group2.A43T.%': 'NA'},
                'cluster3': {'assembled': 'no', 'match': 'no', 'ref_seq': 'NA', 'pct_id': 'NA', 'known_var': 'no', 'novel_var': 'no', 'A42T': 'no', 'A44T.%': 'NA'},
                'cluster4': {'assembled': 'yes', 'match': 'yes', 'ref_seq': 'ref4', 'pct_id': 99.0, 'known_var': 'yes', 'novel_var': 'no', 'group4.A44T': 'het', 'group4.A44T.%': '50.0'},
            },
            'name3': {
                'cluster1': {'assembled': 'yes', 'match': 'yes', 'ref_seq': 'ref_seq42', 'pct_id': 100.0, 'known_var': 'no', 'novel_var': 'no', 'group1.A42T': 'no', 'group1.A42T.%': 'NA'},
                'cluster2': {'assembled': 'no', 'match': 'no', 'ref_seq': 'NA', 'pct_id': 'NA', 'known_var': 'NA', 'novel_var': 'NA', 'group2.A43T': 'NA', 'group2.A43T.%': 'NA'},
                'cluster3': {'assembled': 'no', 'match': 'no', 'ref_seq': 'NA', 'pct_id': 'NA', 'known_var': 'no', 'novel_var': 'no', 'A42T': 'no', 'A44T.%': 'NA'},
                'cluster4': {'assembled': 'yes', 'match': 'yes', 'ref_seq': 'ref4', 'pct_id': 100.0, 'known_var': 'yes', 'novel_var': 'no', 'group4.A44T': 'no', 'group4.A44T.%': 'NA'},
            },
        }

        expected_top_plot_data = {
            'antibio1': {
                'yes': {'cluster1.group1.A42T.cluster4.group4.A44T': [0.125], 'cluster2.group2.A43T.cluster3.interrupted': [0.25]},
                'no': {'cluster1.group1.A42T': [0.125], 'cluster2.group2.A43T.cluster3.interrupted': [0.25]},
                'exclude': {'cluster2.group2.A43T.cluster3.interrupted': [0.25]},
            },
            'antibio2': {
                'yes': {'without_mutation': [0.002], 'cluster2.group2.A43T.cluster3.interrupted': [0.004]},
                'no': {'without_mutation': [0.002], 'cluster2.group2.A43T.cluster3.interrupted': [0.004]},
                'exclude': {'without_mutation': [0.002], 'cluster2.group2.A43T.cluster3.interrupted': [0.004]},
            }
        }

        expected_top_plot_data_no_combs = {
            'antibio1': {
                'yes': {'cluster2.group2.A43T': [0.25], 'cluster4.group4.A44T': [0.125], 'cluster3.interrupted': [0.25], 'cluster1.group1.A42T': [0.125]},
                'no': {'cluster2.group2.A43T': [0.25], 'cluster3.interrupted': [0.25], 'cluster1.group1.A42T': [0.125]},
                'exclude': {'cluster2.group2.A43T': [0.25], 'cluster3.interrupted': [0.25]},
            },
            'antibio2': {
                'yes': {'cluster2.group2.A43T': [0.004], 'without_mutation': [0.002], 'cluster3.interrupted': [0.004]},
                'no': {'cluster2.group2.A43T': [0.004], 'without_mutation': [0.002], 'cluster3.interrupted': [0.004]},
                'exclude': {'cluster2.group2.A43T': [0.004], 'without_mutation': [0.002], 'cluster3.interrupted': [0.004]},
            }
        }

        expected_mutations = {
            'antibio1': {
                'yes': {'cluster1.group1.A42T', 'cluster2.group2.A43T', 'cluster3.interrupted', 'cluster4.group4.A44T'},
                'no': {'cluster1.group1.A42T', 'cluster2.group2.A43T', 'cluster3.interrupted'},
                'exclude': {'cluster2.group2.A43T', 'cluster3.interrupted'},
            },
            'antibio2': {
                'yes': {'without_mutation', 'cluster2.group2.A43T', 'cluster3.interrupted'},
                'no': {'without_mutation', 'cluster2.group2.A43T', 'cluster3.interrupted'},
                'exclude': {'without_mutation', 'cluster2.group2.A43T', 'cluster3.interrupted'},
            }
        }

        expected_combs = {
            'antibio1': {
                'yes': {('cluster2.group2.A43T', 'cluster3.interrupted'), ('cluster1.group1.A42T', 'cluster4.group4.A44T')},
                'no': {('cluster2.group2.A43T', 'cluster3.interrupted'), ('cluster1.group1.A42T',)},
                'exclude': {('cluster2.group2.A43T', 'cluster3.interrupted')}
            },
            'antibio2': {
                'yes': {('without_mutation',), ('cluster2.group2.A43T', 'cluster3.interrupted')},
                'no': {('without_mutation',), ('cluster2.group2.A43T', 'cluster3.interrupted')},
                'exclude': {('without_mutation',), ('cluster2.group2.A43T', 'cluster3.interrupted')}
            }
        }


        expected_no_combs = {
            'antibio1': {
                'yes': {('cluster2.group2.A43T',),  ('cluster3.interrupted',), ('cluster1.group1.A42T',), ('cluster4.group4.A44T',)},
                'no': {('cluster2.group2.A43T',),  ('cluster3.interrupted',), ('cluster1.group1.A42T',)},
                'exclude': {('cluster2.group2.A43T',),  ('cluster3.interrupted',)}
            },
            'antibio2': {
                'yes': {('without_mutation',), ('cluster2.group2.A43T', ), ('cluster3.interrupted',)},
                'no': {('without_mutation',), ('cluster2.group2.A43T', ), ('cluster3.interrupted',)},
                'exclude': {('without_mutation',), ('cluster2.group2.A43T', ), ('cluster3.interrupted',)}
            }
        }


        tmp_tsv = 'tmp.mic_plotter_test.to_boxplot.tsv'

        for antibio in ['antibio1', 'antibio2']:
            for use_het in ['no', 'yes', 'exclude']:
                got_data, got_muts, got_combs = mic_plotter.MicPlotter._get_top_plot_data(summary_data, mic_data, antibio, use_het, interrupted=True, outfile=tmp_tsv)
                expected_tsv = os.path.join(data_dir, 'mic_plotter_to_boxplot_tsv.' + antibio + '.' + use_het + '.tsv')
                self.assertTrue(filecmp.cmp(tmp_tsv, expected_tsv, shallow=False))
                self.assertEqual(got_muts, expected_mutations[antibio][use_het])
                self.assertEqual(got_combs, expected_combs[antibio][use_het])
                self.assertEqual(got_data, expected_top_plot_data[antibio][use_het])
                os.unlink(tmp_tsv)

                got_data, got_muts, got_combs = mic_plotter.MicPlotter._get_top_plot_data(summary_data, mic_data, antibio, use_het, no_combinations=True, interrupted=True, outfile=tmp_tsv)
                expected_tsv = os.path.join(data_dir, 'mic_plotter_to_boxplot_tsv.' + antibio + '.' + use_het +  '.no_combinations.tsv')
                self.assertTrue(filecmp.cmp(tmp_tsv, expected_tsv, shallow=False))
                self.assertEqual(got_muts, expected_mutations[antibio][use_het])
                self.assertEqual(got_combs, expected_no_combs[antibio][use_het])
                self.assertEqual(got_data, expected_top_plot_data_no_combs[antibio][use_het])
                os.unlink(tmp_tsv)


    def test_filter_top_plot_data(self):
        '''test _filter_top_plot_data'''
        top_plot_data = {
            'var1': [1, 2, 3],
            'var2.var3': [1],
            'var1.var3': [1, 2],
        }

        all_mutations = {'var1', 'var2', 'var3'}
        seen_combinations = {('var1',), ('var1', 'var3'), ('var2', 'var3')}

        got_top, got_all, got_seen = mic_plotter.MicPlotter._filter_top_plot_data(top_plot_data, all_mutations, seen_combinations, 1)
        self.assertEqual(got_top, top_plot_data)
        self.assertEqual(got_all, all_mutations)
        self.assertEqual(got_seen, seen_combinations)


        got_top, got_all, got_seen = mic_plotter.MicPlotter._filter_top_plot_data(top_plot_data, all_mutations, seen_combinations, 2)
        expected_top_plot_data = {
            'var1': [1, 2, 3],
            'var1.var3': [1, 2],
        }
        expected_all_mutations = {'var1', 'var3'}
        expected_seen_combinations = {('var1',), ('var1', 'var3')}
        self.assertEqual(got_top, expected_top_plot_data)
        self.assertEqual(got_all, expected_all_mutations)
        self.assertEqual(got_seen, expected_seen_combinations)

        got_top, got_all, got_seen = mic_plotter.MicPlotter._filter_top_plot_data(top_plot_data, all_mutations, seen_combinations, 3)
        expected_top_plot_data = {
            'var1': [1, 2, 3],
        }
        expected_all_mutations = {'var1'}
        expected_seen_combinations = {('var1',),}
        self.assertEqual(got_top, expected_top_plot_data)
        self.assertEqual(got_all, expected_all_mutations)
        self.assertEqual(got_seen, expected_seen_combinations)


    def test_ordered_bottom_plot_rows(self):
        '''test _ordered_bottom_plot_rows'''
        to_order = {'clust1.grp1.42T', 'clust1.grp1.47G', 'clust0.10T', 'abcdefg'}
        got = mic_plotter.MicPlotter._ordered_bottom_plot_rows(to_order)
        expected = ['abcdefg', 'clust0.10T', 'clust1.grp1.42T', 'clust1.grp1.47G']
        self.assertEqual(expected, got)


    def test_ordered_columns(self):
        '''test _ordered_colunns'''
        top_plot_data = {}
        # FIXME


    def test_bottom_scatter_data(self):
        '''test _bottom_scatter_data'''
        #FIXME
        pass


    def test_top_plot_y_ticks(self):
        '''test _top_plot_y_ticks'''
        # FIXME
        pass


    def test_top_plot_scatter_counts(self):
        '''test _top_plot_scatter_counts'''
        top_plot_data = {}
        # FIXME


    def test_top_plot_scatter_data(self):
        '''test _top_plot_scatter_data'''
        top_plot_data = {}
        # FIXME


    def test_top_plot_violin_data(self):
        '''test _top_plot_violin_data'''
        top_plot_data = {}
        # FIXME

