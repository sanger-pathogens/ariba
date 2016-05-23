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
            'only one column. There can NOT be only one\n',
            'two\tcolumns is not enough\n',
            'three\tcolumns\tis still not enough\n',
            'six\tcolumns\tis\tone\ttoo\tmany\n',
        ]

        for line in lines:
            with self.assertRaises(sequence_metadata.Error):
                sequence_metadata.SequenceMetadata(line)

        lines = [
            'gene\tx\tI42L\t.\n',
        ]

        for line in lines:
            with self.assertRaises(sequence_variant.Error):
                sequence_metadata.SequenceMetadata(line)


    def test_init_on_good_input(self):
        '''test init ok on good input'''
        data = sequence_metadata.SequenceMetadata('gene\tn\tI42L\tid\tspam spam wonderful spam')
        self.assertEqual(data.name, 'gene')
        self.assertEqual(data.variant_type, 'n')
        self.assertEqual(data.variant.wild_value, 'I')
        self.assertEqual(data.variant.variant_value, 'L')
        self.assertEqual(data.variant.identifier, 'id')
        self.assertEqual(data.free_text, 'spam spam wonderful spam')


    def test_str(self):
        '''test __str__'''
        lines = [
            'gene1\tn\tA42G\tid1\tspam',
            'gene2\t.\t.\t.\t.',
            'gene3\t.\t.\t.\teggs',
            'gene4\tp\tI42K\tid\tthis mutation kills tardigrades',
        ]

        for line in lines:
            self.assertEqual(line, str(sequence_metadata.SequenceMetadata(line)))


    def test_has_variant(self):
        '''test has_variant'''
        tests = [
            ('gene1\t.\t.\t.', False),
            ('gene1\tn\tA2T\t.', True),
            ('gene1\tn\tT2A\t.', False),
            ('gene1\tp\tI2Y\t.', True),
            ('gene1\tp\tY2I\t.', False),
        ]

        seq = pyfastaq.sequences.Fasta('name', 'ATGTATTGCTGA') # translation: MYC*

        for line, expected in tests:
            metadata = sequence_metadata.SequenceMetadata(line)
            self.assertEqual(expected, metadata.has_variant(seq))
