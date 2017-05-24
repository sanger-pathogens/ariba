import unittest
import os
from ariba import card_record

modules_dir = os.path.dirname(os.path.abspath(card_record.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


class TestCardRecord(unittest.TestCase):
    def test_ARO_id(self):
        '''test _ARO_id'''
        d = {'spam': 'eggs'}
        self.assertEqual(None, card_record.CardRecord._ARO_id(d))
        d['ARO_id'] = '123'
        self.assertEqual('123', card_record.CardRecord._ARO_id(d))


    def test_ARO_accession(self):
        '''test _ARO_accession'''
        d = {'spam': 'eggs'}
        self.assertEqual(None, card_record.CardRecord._ARO_accession(d))
        d['ARO_accession'] = '321'
        self.assertEqual('321', card_record.CardRecord._ARO_accession(d))


    def test_ARO_name(self):
        '''test _ARO_name'''
        d = {'spam': 'eggs'}
        self.assertEqual(None, card_record.CardRecord._ARO_name(d))
        d['ARO_name'] = 'Dave Lister'
        self.assertEqual('Dave Lister', card_record.CardRecord._ARO_name(d))


    def test_ARO_description(self):
        '''test _ARO_description'''
        d = {'spam': 'eggs'}
        self.assertEqual(None, card_record.CardRecord._ARO_description(d))
        d['ARO_description'] = 'Technician, Third Class'
        self.assertEqual('Technician, Third Class', card_record.CardRecord._ARO_description(d))


    def test_ARO_name_to_fasta_name(self):
        '''test _ARO_name_to_fasta_name'''
        tests = [
            ('x', 'x'),
            ('abcD foo bar match at the start', 'abcD'),
            ('foo bar abcD match in the middle', 'abcD'),
            ('match at the end abcD', 'abcD'),
            ('use first three foo bar', 'use_first_three'),
            ('remove.any.dots first three', 'remove_any_dots_first_three')
        ]

        for aro_name, expected in tests:
            got = card_record.CardRecord._ARO_name_to_fasta_name(aro_name)
            self.assertEqual(expected, got)


    def test_dna_seqs_and_genbank_ids(self):
        '''test _dna_seqs_and_genbank_ids'''
        d = {'spam': 'eggs'}
        self.assertEqual([], card_record.CardRecord._dna_seqs_and_genbank_ids(d))
        d['model_sequences'] = {}
        self.assertEqual([], card_record.CardRecord._dna_seqs_and_genbank_ids(d))
        d['model_sequences']['sequence'] = {}
        self.assertEqual([], card_record.CardRecord._dna_seqs_and_genbank_ids(d))
        d['model_sequences']['sequence']['foo'] = {}
        self.assertEqual([], card_record.CardRecord._dna_seqs_and_genbank_ids(d))
        d['model_sequences']['sequence']['foo']['dna_sequence'] = {}
        self.assertEqual([], card_record.CardRecord._dna_seqs_and_genbank_ids(d))
        d['model_sequences']['sequence']['foo']['dna_sequence']['sequence'] = 'ACGT'
        self.assertEqual([], card_record.CardRecord._dna_seqs_and_genbank_ids(d))
        d['model_sequences']['sequence']['foo']['dna_sequence']['accession'] = 'ABC123'
        self.assertEqual([], card_record.CardRecord._dna_seqs_and_genbank_ids(d))
        d['model_sequences']['sequence']['foo']['dna_sequence']['fmin'] = '42'
        self.assertEqual([], card_record.CardRecord._dna_seqs_and_genbank_ids(d))
        d['model_sequences']['sequence']['foo']['dna_sequence']['fmax'] = '4242'
        self.assertEqual([], card_record.CardRecord._dna_seqs_and_genbank_ids(d))
        d['model_sequences']['sequence']['foo']['protein_sequence'] = {}
        d['model_sequences']['sequence']['foo']['protein_sequence']['GI'] = '123456789'
        self.assertEqual([], card_record.CardRecord._dna_seqs_and_genbank_ids(d))
        d['model_sequences']['sequence']['foo']['protein_sequence']['sequence'] = 'III'
        expected = [('foo', '123456789', 'ABC123', '42', '4242', 'ACGT', 'III')]
        self.assertEqual(expected, card_record.CardRecord._dna_seqs_and_genbank_ids(d))
        del d['model_sequences']['sequence']['foo']['protein_sequence']['GI']
        d['model_sequences']['sequence']['foo']['protein_sequence']['accession'] = '123456789'
        self.assertEqual(expected, card_record.CardRecord._dna_seqs_and_genbank_ids(d))

        d['model_sequences']['sequence']['bar'] = {
            'dna_sequence': {
                'sequence': 'TTTT',
                'fmin': '42',
                'fmax': '4242',
                'accession': 'BCD234',
            },
            'protein_sequence': {
                'GI': '',
                'sequence': ''
            }
        }

        expected = [('bar', 'NA', 'BCD234', '42', '4242', 'TTTT', '')] + expected
        self.assertEqual(expected, card_record.CardRecord._dna_seqs_and_genbank_ids(d))


    def test_snps(self):
        '''test snps'''
        d = {'spam': 'eggs'}
        self.assertEqual(set(), card_record.CardRecord._snps(d))
        d['model_param'] = {}
        self.assertEqual(set(), card_record.CardRecord._snps(d))
        d['model_param']['snp'] = {}
        self.assertEqual(set(), card_record.CardRecord._snps(d))
        d['model_param']['snp']['param_value'] = {}
        self.assertEqual(set(), card_record.CardRecord._snps(d))
        d['model_param']['snp']['param_value'] = {
                '1': 'I42L',
                '2': 'S100T',
       }
        self.assertEqual({'I42L', 'S100T'}, card_record.CardRecord._snps(d))


    def test_get_data(self):
        d = {
            'ARO_id': '123',
            'ARO_accession': '1234567',
            'ARO_name': 'ARO_name1',
            'ARO_description': 'ARO description that we want.',
            'model_id': '1',
            'model_name': 'Model_name1',
            'model_type': 'protein homolog model',
            'model_type_id': '12345',
            'model_description': 'Models to detect proteins conferring antibiotic resistance, which include a reference protein sequence and a curated BLASTP cut-off.',
            'model_sequences': {
                'sequence': {
                    '1234': {
                        'protein_sequence': {
                            'sequence': 'MCDE*',
                            'GI': '229597524'
                        },
                        'dna_sequence': {
                            'sequence': 'ATGTGCGATGAATAA',
                            'strand': '+',
                            'fmax': '1194',
                            'fmin': '0',
                            'accession': 'XX0000001'
                        },
                        'NCBI_taxonomy': {
                            'NCBI_taxonomy_cvterm_id': '234567',
                            'NCBI_taxonomy_id': '42',
                            'NCBI_taxonomy_name': 'Genus1 species1'
                        }
                    }
                }
            },
            'model_param': {
                'blastp_evalue': {} # we're ignoring this, so make it empty for tests to save a few lines
            },
            'ARO_category': {
                '36696': {
                    'category_aro_description': 'Enzyme that catalyzes the inactivation of an antibiotic resulting in resistance.  Inactivation includes chemical modification, destruction, etc.',
                    'category_aro_cvterm_id': '36696',
                    'category_aro_accession': '3000557',
                    'category_aro_name': 'antibiotic inactivation enzyme'
                },
                '36268': {
                    'category_aro_description': 'Genes conferring resistance to beta-lactams.',
                    'category_aro_cvterm_id': '36268',
                    'category_aro_accession': '3000129',
                    'category_aro_name': 'beta-lactam resistance gene'
                }
            },
        }

        expected = {
            'ARO_id': '123',
            'ARO_accession': '1234567',
            'ARO_name': 'ARO_name1',
            'ARO_description': 'ARO description that we want.',
            'dna_seqs_and_ids': [(
                '1234',
                '229597524',
                'XX0000001',
                '0',
                '1194',
                'ATGTGCGATGAATAA',
                'MCDE*'
            )],
            'snps': set(),
        }

        record = card_record.CardRecord(d)
        got = record.get_data()
        self.assertEqual(expected, got)

        d['model_param'] = {
            'snp': {
                'param_value': {
                    '1': 'I42L',
                    '2': 'S100T',
                }
            }
        }

        expected['snps'] = {'I42L', 'S100T'}
        record = card_record.CardRecord(d)
        got = record.get_data()
        self.assertEqual(expected, got)

