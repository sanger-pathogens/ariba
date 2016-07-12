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
            'four\tcolumns\tis\talso not enough\n',
            'five\tcolumns\tis\talso\tnot enough\n',
            'seven\tcolumns\tis\tone\tmore\tthan\nwe want',
        ]

        for line in lines:
            with self.assertRaises(sequence_metadata.Error):
                sequence_metadata.SequenceMetadata(line)

        tests = [
            ('gene\tx\t0\t.\t.\tfoo\n', sequence_metadata.Error),
            ('gene\t1\t2\t.\t.\tfoo\n', sequence_metadata.Error),
            ('gene\t1\t1\tI42\t.\tfoo\n', sequence_variant.Error),
        ]

        for line, err in tests:
            with self.assertRaises(err):
                sequence_metadata.SequenceMetadata(line)


    def test_init_on_good_input(self):
        '''test init ok on good input'''
        data = sequence_metadata.SequenceMetadata('gene\t1\t0\tI42L\tid\tspam spam wonderful spam')
        self.assertEqual(data.name, 'gene')
        self.assertEqual(data.seq_type, 'p')
        self.assertEqual(data.variant_only, False)
        self.assertEqual(data.variant.wild_value, 'I')
        self.assertEqual(data.variant.variant_value, 'L')
        self.assertEqual(data.variant.identifier, 'id')
        self.assertEqual(data.free_text, 'spam spam wonderful spam')


    def test_str(self):
        '''test __str__'''
        lines = [
            'gene1\t1\t1\tA42G\tid1\tspam',
            'gene2\t0\t0\t.\t.\t.',
            'gene3\t0\t0\t.\t.\teggs',
            'gene4\t1\t0\tI42K\tid\tthis mutation kills tardigrades',
        ]

        for line in lines:
            self.assertEqual(line, str(sequence_metadata.SequenceMetadata(line)))


    def test_has_variant(self):
        '''test has_variant'''
        tests = [
            ('gene1\t0\t0\t.\t.\t.', False),
            ('gene1\t0\t0\tA2T\t.\t,', True),
            ('gene1\t0\t0\tT2A\t.\t.', False),
            ('gene1\t1\t0\tI2Y\t.\t.', True),
            ('gene1\t1\t0\tY2I\t.\t.', False),
        ]

        seq = pyfastaq.sequences.Fasta('name', 'ATGTATTGCTGA') # translation: MYC*

        for line, expected in tests:
            metadata = sequence_metadata.SequenceMetadata(line)
            self.assertEqual(expected, metadata.has_variant(seq))


    def test_to_string(self):
        '''test to_string'''
        lines = [
            ('gene1', '0', '0', 'A42G', 'id1', 'spam'),
            ('gene2', '0', '0', '.', '.', '.'),
            ('gene3', '0', '0', '.', '.', 'eggs'),
            ('gene4', '1', '0', 'I42K', 'id', 'this mutation kills tardigrades'),
        ]

        for line in lines:
            m = sequence_metadata.SequenceMetadata('\t'.join(line))
            for separator in ('_', '\t'):
                expected = separator.join(line)
                self.assertEqual(expected, m.to_string(separator=separator))

