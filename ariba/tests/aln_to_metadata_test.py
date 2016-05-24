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


    def test_check_insertion_coords(self):
        '''test _check_insertion_coords'''
        seq = pyfastaq.sequences.Fasta('name', 'AAA---GGG------TTT---')
        self.assertTrue(aln_to_metadata.AlnToMetadata._check_insertion_coords(seq))

        bad_seqs = [
            pyfastaq.sequences.Fasta('name', 'AAA--GGG'),  # bad length
            pyfastaq.sequences.Fasta('name', 'A---AA'),  # bad start position
            pyfastaq.sequences.Fasta('name', 'AA---AA'), # bad start position
        ]

        for seq in bad_seqs:
            with self.assertRaises(aln_to_metadata.Error):
                aln_to_metadata.AlnToMetadata._check_insertion_coords(seq)


    def test_check_coding_seq(self):
        '''test _check_coding_seq'''
        seq = pyfastaq.sequences.Fasta('name', 'ATGCTTTAG')
        self.assertTrue(aln_to_metadata.AlnToMetadata._check_coding_seq(seq))

        bad_seqs = [
            pyfastaq.sequences.Fasta('name', 'TTGCTTAG'), # length not a mutliple of 3
            pyfastaq.sequences.Fasta('name', 'TTTCTTTAG'), # no start codon
            pyfastaq.sequences.Fasta('name', 'ATGTAGCTTTAG'), # stop codon in middle
            pyfastaq.sequences.Fasta('name', 'TTGCTTTTT'), # no stop at end
        ]

        for seq in bad_seqs:
            with self.assertRaises(aln_to_metadata.Error):
                aln_to_metadata.AlnToMetadata._check_coding_seq(seq)


    def test_check_sequences_non_coding(self):
        '''test _check_sequences with noncoding seqs'''
        padded_sequences = {
            'seq1': pyfastaq.sequences.Fasta('seq1', 'AC-T')
        }

        unpadded_sequences = aln_to_metadata.AlnToMetadata._make_unpadded_seqs(padded_sequences)
        self.assertTrue(aln_to_metadata.AlnToMetadata._check_sequences(padded_sequences, unpadded_sequences, False))
        padded_sequences['seq2'] = pyfastaq.sequences.Fasta('seq2', 'AC-')
        unpadded_sequences = aln_to_metadata.AlnToMetadata._make_unpadded_seqs(padded_sequences)
        with self.assertRaises(aln_to_metadata.Error):
            aln_to_metadata.AlnToMetadata._check_sequences(padded_sequences, unpadded_sequences, False)


    def test_check_sequences_coding(self):
        '''test _check_sequences with coding seqs'''
        padded_sequences = {
            'seq1': pyfastaq.sequences.Fasta('seq1', 'ATGCTTTAG'),
            'seq2': pyfastaq.sequences.Fasta('seq2', 'ATG---TAG')
        }

        unpadded_sequences = aln_to_metadata.AlnToMetadata._make_unpadded_seqs(padded_sequences)

        self.assertTrue(aln_to_metadata.AlnToMetadata._check_sequences(padded_sequences, unpadded_sequences, True))

        bad_seqs = [
            'ATGCTTAG', # length not a mutliple of 3
            'TTTCTTTAG', # no start codon
            'ATGTAGCTTTAG', # stop codon in middle
            'ATGTTTTTT', # no stop at end
            'ATGC---TTTAG', # bad insertion
            'ATGCT---TTAG', # bad insertion
            'ATG-CTTTAG', # bad insertion
            'ATG--CTTTAG', # bad insertion
            'ATG----CTTTAG', # bad insertion
        ]

        for seq in bad_seqs:
            padded_sequences['seq2'] = pyfastaq.sequences.Fasta('seq2', seq)
            unpadded_sequences = aln_to_metadata.AlnToMetadata._make_unpadded_seqs(padded_sequences)
            with self.assertRaises(aln_to_metadata.Error):
                aln_to_metadata.AlnToMetadata._check_sequences(padded_sequences, unpadded_sequences, True)
