import unittest
import os
import pyfastaq
from ariba import sequence_metadata, sequence_variant

modules_dir = os.path.dirname(os.path.abspath(sequence_metadata.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


class TestSequenceMetadata(unittest.TestCase):
    def test_init_fails_on_bad_lines(self):
        '''Test init fails on bad lines'''
        lines = [
            'only one column\n',
            'two\tcolumns is not enough\n',
            'two\tcolumns\tis not enough\n',
            'six\tcolumns\tis\tone\ttoo\tmany\n',
            'name\tp\tI42L\tx', # column 4 should be "Y" or "N"
        ]

        for line in lines:
            with self.assertRaises(sequence_metadata.Error):
                sequence_metadata.SequenceMetadata(line)

        with self.assertRaises(sequence_variant.Error):
            sequence_metadata.SequenceMetadata('gene\tx\tI42L\tN\n')


    def test_init_on_good_input(self):
        '''test init ok on good input'''
        data = sequence_metadata.SequenceMetadata('gene\tn\tI42L\tN\tspam spam wonderful spam')
        self.assertEqual(data.name, 'gene')
        self.assertEqual(data.variant_type, 'n')
        self.assertEqual(data.variant.wild_value, 'I')
        self.assertEqual(data.variant.variant_value, 'L')
        self.assertFalse(data.always_report)
        self.assertEqual(data.free_text, 'spam spam wonderful spam')


    def test_str(self):
        '''test __str__'''
        lines = [
            'gene1\tn\tA42G\tY\tspam',
            'gene2\t.\t.\tN',
            'gene3\t.\t.\tY\teggs',
            'gene4\tp\tI42K\tN\tthis mutation kills tardigrades',
        ]

        for line in lines:
            self.assertEqual(line, str(sequence_metadata.SequenceMetadata(line)))


    def test_has_variant(self):
        '''test has_variant'''
        tests = [
            ('gene1\t.\t.\tN', False),
            ('gene1\tn\tA2T\tN', True),
            ('gene1\tn\tT2A\tN', False),
            ('gene1\tp\tI2Y\tN', True),
            ('gene1\tp\tY2I\tN', False),
        ]

        seq = pyfastaq.sequences.Fasta('name', 'ATGTATTGCTGA') # translation: MYC*

        for line, expected in tests:
            metadata = sequence_metadata.SequenceMetadata(line)
            self.assertEqual(expected, metadata.has_variant(seq))
