import unittest
import os
from ariba import tb

modules_dir = os.path.dirname(os.path.abspath(tb.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


class TestTb(unittest.TestCase):
    def test_report_to_resistance_dict(self):
        '''test report_to_resistance_dict'''
        infile = os.path.join(data_dir, 'tb_report_to_resistance_dict.tsv')
        got = tb.report_to_resistance_dict(infile)
        expect = {
            'Ethambutol': [('ref1', 'I42J'), ('ref2', 'R10S')],
            'Rifampicin': [('ref1', 'I42J')],
            'Isoniazid': [('katG', 'Incomplete_gene')],
        }
        self.assertEqual(expect, got)


    def test_genbank_to_gene_coords(self):
        '''test genbank_to_gene_coords'''
        infile = os.path.join(data_dir, 'tb_genbank_to_gene_coords.gb')
        expect = {
            'gene1': {'start': 42, 'end': 45},
            'gene2': {'start': 10, 'end': 5},
        }
        genes = {'gene1', 'gene2'}

        got = tb.genbank_to_gene_coords(infile, genes)
        self.assertEqual(expect, got)

