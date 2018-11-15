import filecmp
import json
import os
import unittest

import pyfastaq

from ariba import common, tb

modules_dir = os.path.dirname(os.path.abspath(tb.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


class TestTb(unittest.TestCase):
    def test_report_to_resistance_dict(self):
        '''test report_to_resistance_dict'''
        infile = os.path.join(data_dir, 'tb_report_to_resistance_dict.tsv')
        got = tb.report_to_resistance_dict(infile)
        expect = {
            'Ethambutol': [('ref1', 'I42J'), ('ref2_upstream', 'A-3G')],
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



    def test_load_mutations(self):
        '''test load_mutations'''
        mutation_to_drug_json = 'tmp.tb.load_mutations.json'
        variants_txt = 'tmp.tb.load_mutations.txt'

        gene_coords = {
            'gene1': {'start': 42, 'end': 45},
            'gene2': {'start': 100, 'end': 50},
            'gene3': {'start': 200, 'end': 300},
        }

        mutation_to_drug = {
            'gene1_A-5T': ['drug1'],
            'gene1_A5C': ['drug1', 'drug2'],
            'gene2_A-3G': ['drug4'],
            'gene2_A10C': ['drug5'],
            'gene3_AC10A': ['drug5'],
        }
        with open(mutation_to_drug_json, 'w') as f:
            json.dump(mutation_to_drug, f, indent=4, sort_keys=True)

        with open(variants_txt, 'w') as f:
            print('gene1', 'A-5T', 'DNA', sep='\t', file=f)
            print('gene1', 'A5C', 'PROT', sep='\t', file=f)
            print('gene2', 'A-3G', 'DNA', sep='\t', file=f)
            print('gene2', 'A10C', 'PROT', sep='\t', file=f)
            print('gene3', 'AC10A', 'DNA', sep='\t', file=f)
        got_mutations, got_genes_with_indels, got_genes_need_upstream, got_genes_non_upstream = tb.load_mutations(gene_coords, mutation_to_drug_json, variants_txt, upstream_before=10)
        expect_mutations = [
            {'gene': 'gene1', 'var': 'A6T', 'coding': 1, 'upstream': True, 'drugs': 'drug1', 'original_mutation': 'A-5T'},
            {'gene': 'gene1', 'var': 'A5C', 'coding': 1, 'upstream': False, 'drugs': 'drug1,drug2'},
            {'gene': 'gene2', 'var': 'A8G', 'coding': 1, 'upstream': True, 'drugs': 'drug4', 'original_mutation': 'A-3G'},
            {'gene': 'gene2', 'var': 'A10C', 'coding': 1, 'upstream': False, 'drugs': 'drug5'}
        ]
        expect_genes_with_indels = {'gene3'}
        expect_genes_need_upstream = {'gene1', 'gene2'}
        expect_genes_non_upstream = {'gene1', 'gene2'}
        self.assertEqual(expect_mutations, got_mutations)
        self.assertEqual(expect_genes_with_indels, got_genes_with_indels)
        self.assertEqual(expect_genes_need_upstream, got_genes_need_upstream)
        self.assertEqual(expect_genes_non_upstream, got_genes_non_upstream)
        os.unlink(mutation_to_drug_json)
        os.unlink(variants_txt)


    def test_write_prepareref_fasta_file(self):
        '''test write_prepareref_fasta_file'''
        outfile = 'tmp.write_prepareref_fasta_file.fa'
        # start of the ref seq is:
        # 0         10        20
        # 0123456789012345678901234567890
        # TTGACCGATGACCCCGGTTCAGGCTTCACCA
        gene_coords = {
            'gene1': {'start': 0, 'end': 6},
            'gene2': {'start': 12, 'end': 10},
            'gene3': {'start': 13, 'end': 16},

        }
        expect_seqs = {
            'gene1': pyfastaq.sequences.Fasta('gene1', 'TTGACCG'),
            'gene2': pyfastaq.sequences.Fasta('gene2', 'GGT'),
            'gene2_upstream': pyfastaq.sequences.Fasta('gene2_upstream', 'CGGGGTCA'),
            'gene3_upstream': pyfastaq.sequences.Fasta('gene3_upstream', 'ACCCCGGT'),
        }
        genes_need_upstream = {'gene2', 'gene3'}
        genes_non_upstream = {'gene1', 'gene2'}
        tb.write_prepareref_fasta_file(outfile, gene_coords, genes_need_upstream, genes_non_upstream, upstream_before=3, upstream_after=5)
        got_seqs = {}
        pyfastaq.tasks.file_to_dict(outfile, got_seqs)
        os.unlink(outfile)
        self.assertEqual(expect_seqs, got_seqs)


    def test_write_prepareref_metadata_file(self):
        '''test write_prepareref_metadata_file'''
        mutations = [
            {'gene': 'gene1', 'var': 'A6T', 'coding': 1, 'upstream': True, 'drugs': 'drug1', 'original_mutation': 'A-5T'},
            {'gene': 'gene1', 'var': 'A5C', 'coding': 1, 'upstream': False, 'drugs': 'drug1,drug2'},
            {'gene': 'gene2', 'var': 'A8G', 'coding': 1, 'upstream': True, 'drugs': 'drug4', 'original_mutation': 'A-3G'},
            {'gene': 'gene2', 'var': 'A10C', 'coding': 1, 'upstream': False, 'drugs': 'drug5'}
        ]
        outfile = 'tmp.write_prepareref_metadata_file.tsv'
        tb.write_prepareref_metadata_file(mutations, outfile)
        expected = os.path.join(data_dir, 'tb_write_prepareref_metadata_file.tsv')
        self.assertTrue(filecmp.cmp(expected, outfile, shallow=False))
        os.unlink(outfile)


    def test_make_prepareref_files(self):
        '''test make_prepareref_files'''
        outprefix = 'tmp.make_prepareref_files'
        expected_files = [outprefix + '.fa', outprefix + '.tsv']
        for fname in expected_files:
            if os.path.exists(fname):
                os.unlink(fname)
        tb.make_prepareref_files(outprefix)

        for fname in expected_files:
            self.assertTrue(os.path.exists(fname))
            os.unlink(fname)


    def test_make_prepareref_dir(self):
        '''test make_prepareref_dir'''
        outdir = 'tmp.make_prepareref_dir'
        common.rmtree(outdir)
        tb.make_prepareref_dir(outdir)
        self.assertTrue(os.path.exists(outdir))
        json_file = os.path.join(outdir, '00.params.json')
        common.rmtree(outdir)

