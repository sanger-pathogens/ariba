import unittest
import filecmp
import os
from ariba import vfdb_parser

modules_dir = os.path.dirname(os.path.abspath(vfdb_parser.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


class TestVfdbParser(unittest.TestCase):
    def test_fa_header_to_name_pieces(self):
        '''test _fa_header_to_name_pieces'''
        tests = [
            ('VF123(gi:1234) (abcD) foobar description [abc] [genus species]', ('VF123(gi:1234)', 'abcD', 'foobar description [abc]', 'genus species')),
            ('F123(gi:1234) (abcD) foobar description [abc] [genus species]', None), # no V at start
            ('VF123(gi:1234) (abcD) foobar description [abc]', None),  # missing genus species
            ('VF123(gi:1234) abcD foobar description [abc] [genus species]', None), # brackets missing around abcD
        ]

        for header, expected in tests:
            got = vfdb_parser.VfdbParser._fa_header_to_name_pieces(header)
            self.assertEqual(expected, got)


    def test_fa_header_to_name_and_metadata(self):
        '''test _fa_header_to_name_and_metadata'''
        headers = [
            'VF123(gi:1234) (abcD) foobar description [abc] [genus species]',
            'F123(gi:1234) (abcD) foobar description [abc] [genus species]', # no V at start
            'VF123(gi:1234) (abcD) foobar description [abc]', # missing genus species
            'VF123(gi:1234) abcD foobar description [abc] [genus species]', # brackets missing around abcD
        ]

        expected = [
            ('abcD.VF123(gi:1234).genus_species', 'foobar description [abc]'),
            (headers[1], '.'),
            (headers[2], '.'),
            (headers[3], '.'),
        ]

        assert len(headers) == len(expected)
        for i in range(len(headers)):
            got = vfdb_parser.VfdbParser._fa_header_to_name_and_metadata(headers[i])
            self.assertEqual(expected[i], got)


    def test_run(self):
        '''test run'''
        infile = os.path.join(data_dir, 'vfdb_parser_test_run.in.fa')
        expected_tsv = os.path.join(data_dir, 'vfdb_parser_test_run.out.tsv')
        expected_fa = os.path.join(data_dir, 'vfdb_parser_test_run.out.fa')
        outprefix = 'tmp.vfdb_parser_test_run'
        got_tsv = outprefix + '.tsv'
        got_fa = outprefix + '.fa'
        vp = vfdb_parser.VfdbParser(infile, outprefix)
        vp.run()
        self.assertTrue(filecmp.cmp(expected_tsv, got_tsv, shallow=False))
        self.assertTrue(filecmp.cmp(expected_fa, got_fa, shallow=False))
        os.unlink(got_tsv)
        os.unlink(got_fa)
