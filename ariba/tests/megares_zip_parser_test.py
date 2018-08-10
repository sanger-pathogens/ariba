import unittest
import filecmp
import os
import pyfastaq
from ariba import common, megares_zip_parser

modules_dir = os.path.dirname(os.path.abspath(megares_zip_parser.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


class TestMegaresZipParser(unittest.TestCase):
    def test_extract_files_ok(self):
        '''test _extract_files when all ok'''
        zip_file = os.path.join(data_dir, 'megares_zip_parse_extract_files_ok.zip')
        tmp_dir = 'tmp.test_megares_extract_files_ok'
        got = megares_zip_parser.MegaresZipParser._extract_files(zip_file, tmp_dir)
        common_dir = os.path.join('megares_zip_parse_extract_files_ok', 'megares_v1.01')
        expected = {
            'annotations': os.path.join(common_dir, 'megares_annotations_v1.01.csv'),
            'fasta': os.path.join(common_dir, 'megares_database_v1.01.fasta'),
            'header_mappings': os.path.join(common_dir, 'megares_to_external_header_mappings_v1.01.tsv')
        }

        self.assertEqual(expected, got)

        for filename in expected.values():
            self.assertTrue(os.path.exists(os.path.join(tmp_dir, filename)))

        common.rmtree(tmp_dir)


    def test_extract_files_one_missing(self):
        '''test _extract_files when one missing'''
        zip_file = os.path.join(data_dir, 'megares_zip_parse_extract_files_one_missing.zip')
        tmp_dir = 'tmp.test_megares_extract_files_one_missing'
        with self.assertRaises(megares_zip_parser.Error):
            got = megares_zip_parser.MegaresZipParser._extract_files(zip_file, tmp_dir)


    def test_load_annotations_file(self):
        '''test _load_annotations_file'''
        infile = os.path.join(data_dir, 'megares_zip_parser_load_annotations.csv')
        expected = {
            'Bla|OXA-1|JN123456|42-141|100|betalactams|Class_A_betalactamases|OXA': {'class': 'betalactams', 'mechanism': 'Class A betalactamases', 'group': 'OXA'},
            'Foo|Bar-1|JN42|1-11|10|foobar|Class_foobar|Bar': {'class': 'Class', 'mechanism': 'foobar', 'group': 'Bar'}
        }
        got = megares_zip_parser.MegaresZipParser._load_annotations_file(infile)
        self.maxDiff = None
        self.assertEqual(expected, got)


    def test_load_header_mappings_file(self):
        '''test _load_header_mappings_file'''
        infile = os.path.join(data_dir, 'megares_zip_parser_load_header_mappings.tsv')
        expected = {
            'Bla|OXA-1|JN123456|42-141|100|betalactams|Class_A_betalactamases|OXA': {'Source_Database': 'SOURCE1', 'Source_Headers(space_separated)': 'source header 1'},
            'Foo|Bar-1|JN42|1-11|10|foobar|Class_foobar|Bar': {'Source_Database': 'SOURCE2', 'Source_Headers(space_separated)': 'source header 2'},
        }
        got = megares_zip_parser.MegaresZipParser._load_header_mappings_file(infile)
        self.maxDiff = None
        self.assertEqual(expected, got)


    def test_write_files(self):
        '''test _write_files'''
        fasta_file = os.path.join(data_dir, 'megares_zip_parser_write_files', 'megares_database_v1.01.fasta')
        annotations_file = os.path.join(data_dir, 'megares_zip_parser_write_files', 'megares_annotations_v1.01.csv')
        mappings_file = os.path.join(data_dir, 'megares_zip_parser_write_files', 'megares_to_external_header_mappings_v1.01.tsv')
        sequences = {}
        pyfastaq.tasks.file_to_dict(fasta_file, sequences)
        annotation_data = megares_zip_parser.MegaresZipParser._load_annotations_file(annotations_file)
        mappings_data = megares_zip_parser.MegaresZipParser._load_header_mappings_file(mappings_file)

        tmp_prefix = 'tmp.test_megares_zip_parser_write_files'
        megares_zip_parser.MegaresZipParser._write_files(tmp_prefix, sequences, annotation_data, mappings_data)

        expected_fasta = os.path.join(data_dir, 'megares_zip_parser_write_files.expect.fa')
        expected_tsv = os.path.join(data_dir, 'megares_zip_parser_write_files.expect.tsv')
        self.assertTrue(filecmp.cmp(expected_fasta, tmp_prefix + '.fa', shallow=False))
        self.assertTrue(filecmp.cmp(expected_tsv, tmp_prefix + '.tsv', shallow=False))
        os.unlink(tmp_prefix + '.fa')
        os.unlink(tmp_prefix + '.tsv')
