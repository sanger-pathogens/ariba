import unittest
import os
from ariba import sequence_variant

modules_dir = os.path.dirname(os.path.abspath(sequence_variant.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


class TestSequenceVariant(unittest.TestCase):
    def test_init_fails(self):
        '''Test init fails'''
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
                v = sequence_variant.Variant(var)


    def test_init_ok(self):
        '''Test init ok'''
        variants = ['I42K', 'i42k', 'I42k', 'i42K']

        for var in variants:
            aa_var = sequence_variant.Variant(var)
            self.assertEqual(41, aa_var.position)
            self.assertEqual('I', aa_var.wild_aa)
            self.assertEqual('K', aa_var.variant_aa)


    def test_init_str(self):
        '''Test init ok and str'''
        variants = ['I42K', 'i42k', 'I42k', 'i42K']
        expected = 'I42K'

        for var in variants:
            self.assertEqual(expected, str(sequence_variant.Variant(var)))


    def test_agrees_with_protein_seq(self):
        '''test agrees_with_protein_seq'''
        protein_seq = 'BrissSpecialStvff'
        tests = [
            ('I3K', True),
            ('K3I', True),
            ('A2b', False),
            ('x1000y', False)
        ]

        for var, expected in tests:
            variant = sequence_variant.Variant(var)
            self.assertEqual(expected, variant.agrees_with_protein_seq(protein_seq))

