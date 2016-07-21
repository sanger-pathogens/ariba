import unittest
import os
import pyfastaq
from ariba import sequence_variant

modules_dir = os.path.dirname(os.path.abspath(sequence_variant.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


class TestSequenceVariant(unittest.TestCase):
    def test_init_fails_on_bad_variant_strings(self):
        '''Test init fails on bad variant strings'''
        bad_variants = [
            'x',
            'x1',
            '1x',
            '1x1',
            'I42K43',
            'I-1K',
        ]

        for var in bad_variants:
            with self.assertRaises(sequence_variant.Error):
                sequence_variant.Variant('p', var, '.')


    def test_init_ok(self):
        '''Test init ok'''
        variants = [('I42K', '.'), ('i42k', 'id1'), ('I42k', 'id2'), ('i42K', 'id3')]

        for var, identifier in variants:
            aa_var = sequence_variant.Variant('p', var, identifier)
            self.assertEqual(41, aa_var.position)
            self.assertEqual('I', aa_var.wild_value)
            self.assertEqual('K', aa_var.variant_value)
            if identifier == '.':
                self.assertIsNone(aa_var.identifier)
            else:
                self.assertEqual(identifier, aa_var.identifier)


    def test_init_str(self):
        '''Test init ok and str'''
        variants = ['I42K', 'i42k', 'I42k', 'i42K']
        expected = 'I42K'

        for var in variants:
            self.assertEqual(expected, str(sequence_variant.Variant('p', var, '.')))


    def test_sanity_check_against_seq_no_translate(self):
        '''test sanity_check_against_seq with translate False'''
        seq = 'BrissSpecialStvff'
        tests = [
            ('I3K', True),
            ('K3I', True),
            ('A2b', False),
            ('x1000y', False)
        ]

        for var, expected in tests:
            variant = sequence_variant.Variant('p', var, '.')
            self.assertEqual(expected, variant.sanity_check_against_seq(seq))


    def test_sanity_check_against_seq_translate(self):
        '''test sanity_check_against_seq with translate True'''
        seq = 'AGTACGACGTAC'  # translates to STTY
        tests = [
            ('S1X', True),
            ('x1s', True),
            ('a1y', False),
            ('x5y', False)
        ]

        for var, expected in tests:
            variant = sequence_variant.Variant('p', var, '.')
            self.assertEqual(expected, variant.sanity_check_against_seq(seq, translate_seq=True))


    def test_has_variant(self):
        '''test has_variant'''
        seq = pyfastaq.sequences.Fasta('name', 'ATGTATTGCTGA') # translation: MYC*
        tests = [
            (sequence_variant.Variant('n', 'A2T', '.'), True),
            (sequence_variant.Variant('n', 'T2A', '.'), False),
            (sequence_variant.Variant('p', 'I2Y', '.'), True),
            (sequence_variant.Variant('p', 'Y2I', '.'), False),
        ]

        for var, expected in tests:
            self.assertEqual(expected, var.has_variant(seq))


    def test_nucleotide_range(self):
        '''test nucleotide_range'''
        sv = sequence_variant.Variant('n', 'A2T', '.')
        self.assertEqual((1, 1), sv.nucleotide_range())
        sv = sequence_variant.Variant('p', 'I42L', '.')
        self.assertEqual((123, 125), sv.nucleotide_range())
