import unittest
import os
from ariba import sequence_metadata, sequence_variant

modules_dir = os.path.dirname(os.path.abspath(sequence_metadata.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


class TestSequenceMetadata(unittest.TestCase):
    def test_init_fails_on_bad_lines(self):
        '''Test init fails on bad lines'''
        lines = [
            'only one column\n',
            'two\tcolumns is not enough\n',
            'five\tcolumns\tis\ttoo\tmany\n',
        ]

        for line in lines:
            with self.assertRaises(sequence_metadata.Error):
                sequence_metadata.SequenceMetadata(line)

        with self.assertRaises(sequence_variant.Error):
            sequence_metadata.SequenceMetadata('gene\tx\tI42L\n')


    def test_init_on_good_input(self):
        '''test init ok on good input'''
        data = sequence_metadata.SequenceMetadata('gene\tn\tI42L\tspam spam wonderful spam')
        self.assertEqual(data.name, 'gene')
        self.assertEqual(data.variant.variant_type, sequence_variant.NUCLEOTIDE_SNP)
        self.assertEqual(data.variant.wild_value, 'I')
        self.assertEqual(data.variant.variant_value, 'L')
        self.assertEqual(data.free_text, 'spam spam wonderful spam')
