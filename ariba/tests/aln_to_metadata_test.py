import unittest
import os
import copy
import shutil
import filecmp
import pyfastaq
from ariba import aln_to_metadata, sequence_variant

modules_dir = os.path.dirname(os.path.abspath(aln_to_metadata.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


class TestAlnToMetadata(unittest.TestCase):
    def test_load_aln_file(self):
        '''test _load_aln_file'''
        aln_file = os.path.join(data_dir, 'aln_to_metadata_load_aln_file.in.fa')
        expected = {
            'seq1': pyfastaq.sequences.Fasta('seq1', 'ABC-DE'),
            'seq2': pyfastaq.sequences.Fasta('seq2', 'ABCQDE'),
        }
        got = aln_to_metadata.AlnToMetadata._load_aln_file(aln_file)
        self.assertEqual(expected, got)


    def test_load_vars_file_good_file(self):
        '''test _load_vars_file good input file'''
        infile = os.path.join(data_dir, 'aln_to_metadata_load_vars_file_good.tsv')
        variant1 = sequence_variant.Variant('p', 'A42B', 'id1')
        variant2 = sequence_variant.Variant('p', 'C43D', 'id2')
        variant3 = sequence_variant.Variant('p', 'E100F', 'id3')
        expected = {
            'seq1': [(variant1, 'description 1')],
            'seq2': [(variant2, 'description 2'), (variant3, 'description 3')]
        }
        got = aln_to_metadata.AlnToMetadata._load_vars_file(infile, True)
        self.assertEqual(expected, got)


    def test_load_vars_bad_files(self):
        '''test _load_vars_file bad input files'''
        infiles = [
            os.path.join(data_dir, 'aln_to_metadata_load_vars_file_bad.1.tsv'),
            os.path.join(data_dir, 'aln_to_metadata_load_vars_file_bad.2.tsv')
        ]

        for infile in infiles:
            with self.assertRaises(aln_to_metadata.Error):
                aln_to_metadata.AlnToMetadata._load_vars_file(infile, True)


    def test_make_unpadded_seqs(self):
        '''test _make_unpadded_seqs'''
        padded = {
            'seq1': pyfastaq.sequences.Fasta('seq1', 'acg---t'),
            'seq2': pyfastaq.sequences.Fasta('seq2', '---a-cgt-'),
        }
        expected = {
            'seq1': pyfastaq.sequences.Fasta('seq1', 'acgt'),
            'seq2': pyfastaq.sequences.Fasta('seq2', 'acgt'),
        }
        got = aln_to_metadata.AlnToMetadata._make_unpadded_seqs(padded)
        self.assertEqual(expected, got)


    def test_check_seq_lengths_same(self):
        '''test _check_seq_lengths_same'''
        seqs = {
            'seq1': pyfastaq.sequences.Fasta('seq1', 'acgt'),
            'seq2': pyfastaq.sequences.Fasta('seq2', 'acgt'),
        }

        self.assertTrue(aln_to_metadata.AlnToMetadata._check_seq_lengths_same(seqs))
        seqs['seq1'].seq = 'a'
        with self.assertRaises(aln_to_metadata.Error):
            aln_to_metadata.AlnToMetadata._check_seq_lengths_same(seqs)


    def test_insertion_coords(self):
        '''test _insertion_coords'''
        ivl = pyfastaq.intervals.Interval
        tests = [
            ('acgt', []),
            ('-a', [pyfastaq.intervals.Interval(0, 0)]),
            ('a---cgt--', [pyfastaq.intervals.Interval(1, 3), pyfastaq.intervals.Interval(7, 8)]),
        ]

        for seq, expected in tests:
            fa = pyfastaq.sequences.Fasta('x', seq)
            got = aln_to_metadata.AlnToMetadata._insertion_coords(fa)
            self.assertEqual(expected, got)
