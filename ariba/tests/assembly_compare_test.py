import unittest
import os
import copy
import shutil
import filecmp
import pyfastaq
import pysam
import pymummer
from ariba import assembly_compare

modules_dir = os.path.dirname(os.path.abspath(assembly_compare.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


class TestAssemblyCompare(unittest.TestCase):
    def test_parse_nucmer_coords_file(self):
        '''test _parse_nucmer_coords_file'''
        coords_file = os.path.join(data_dir, 'assembly_compare_parse_nucmer_coords_file.coords')
        ref_name = 'ref'
        got = assembly_compare.AssemblyCompare._parse_nucmer_coords_file(coords_file, ref_name)
        line1 = ['1', '1000', '1', '1000', '1000', '1000', '100.00', '1000', '1000', '1', '1', 'ref', 'contig1']
        line2 = ['1', '240', '1', '240', '240', '240', '100.00', '1000', '580', '1', '1', 'ref', 'contig2']
        line3 = ['661', '1000', '241', '580', '340', '340', '100.00', '1000', '580', '1', '1', 'ref', 'contig2']
        expected = {
            'contig1': [pymummer.alignment.Alignment('\t'.join(line1))],
            'contig2': [pymummer.alignment.Alignment('\t'.join(line2)), pymummer.alignment.Alignment('\t'.join(line3))],
        }
        self.assertEqual(expected, got)


    def test_nucmer_hits_to_percent_identity(self):
        '''test _nucmer_hits_to_percent_identity'''
        hits = [
            ['1', '10', '1', '10', '10', '10', '90.00', '1000', '1000', '1', '1', 'ref', 'scaff1'],
            ['9', '42', '9', '42', '34', '34', '100.00', '1000', '1000', '1', '1', 'ref', 'scaff1'],
            ['1', '42', '1', '42', '42', '42', '42.42', '1000', '1000', '1', '1', 'ref', 'scaff2'],
        ]
        nucmer_hits = {
            'scaff1': [
                pymummer.alignment.Alignment('\t'.join(hits[0])),
                pymummer.alignment.Alignment('\t'.join(hits[1])),
            ],
            'scaff2': [
                pymummer.alignment.Alignment('\t'.join(hits[2])),
            ]
        }
        expected = {'scaff1': round((90*10 + 100*34) / (10+34), 2), 'scaff2': 42.42}
        got = assembly_compare.AssemblyCompare._nucmer_hits_to_percent_identity(nucmer_hits)
        self.assertEqual(expected, got)


    def test_nucmer_hits_to_assembly_coords(self):
        '''test _nucmer_hits_to_assembly_coords'''
        hits = [
            ['1', '10', '1', '10', '10', '10', '100.00', '1000', '1000', '1', '1', 'ref', 'scaff1'],
            ['9', '42', '9', '42', '34', '34', '100.00', '1000', '1000', '1', '1', 'ref', 'scaff1'],
            ['50', '52', '50', '52', '3', '3', '100.00', '1000', '1000', '1', '1', 'ref', 'scaff1'],
            ['1', '42', '1', '42', '42', '42', '100.00', '1000', '1000', '1', '1', 'ref', 'scaff2'],
        ]
        nucmer_hits = {
            'scaff1': [
                pymummer.alignment.Alignment('\t'.join(hits[0])),
                pymummer.alignment.Alignment('\t'.join(hits[1])),
                pymummer.alignment.Alignment('\t'.join(hits[2])),
            ],
            'scaff2': [
                pymummer.alignment.Alignment('\t'.join(hits[3])),
            ]
        }
        expected = {
            'scaff1': [
                pyfastaq.intervals.Interval(0, 41),
                pyfastaq.intervals.Interval(49, 51)
            ],
            'scaff2': [
                pyfastaq.intervals.Interval(0, 41),
            ]
        }
        got = assembly_compare.AssemblyCompare._nucmer_hits_to_assembly_coords(nucmer_hits)
        self.assertEqual(expected, got)


    def test_nucmer_hits_to_ref_coords(self):
        '''test nucmer_hits_to_ref_coords'''
        hits = [
            ['1', '42', '1', '42', '42', '42', '100.00', '1000', '1000', '1', '1', 'ref', 'contig1'],
            ['31', '52', '1', '22', '22', '22', '100.00', '1000', '1000', '1', '1', 'ref', 'contig1'],
            ['100', '142', '200', '242', '42', '42', '99.42', '1000', '1000', '1', '1', 'ref', 'contig1'],
            ['100', '110', '200', '210', '11', '11', '99.42', '1000', '1000', '1', '1', 'ref', 'contig2'],
        ]
        nucmer_hits = {
            'contig1': [
                pymummer.alignment.Alignment('\t'.join(hits[0])),
                pymummer.alignment.Alignment('\t'.join(hits[1])),
                pymummer.alignment.Alignment('\t'.join(hits[2])),
            ],
            'contig2': [
                pymummer.alignment.Alignment('\t'.join(hits[3])),
            ]
        }
        got = assembly_compare.AssemblyCompare.nucmer_hits_to_ref_coords(nucmer_hits)
        expected = {
            'contig1': [pyfastaq.intervals.Interval(0,51), pyfastaq.intervals.Interval(99, 141)],
            'contig2': [pyfastaq.intervals.Interval(99, 109)]
        }

        self.assertEqual(expected, got)

        got = assembly_compare.AssemblyCompare.nucmer_hits_to_ref_coords(nucmer_hits, contig='contig2')
        del expected['contig1']
        self.assertEqual(expected, got)


    def test_ref_cov_per_contig(self):
        '''test ref_cov_per_contig'''
        hits = [
            ['1', '42', '1', '42', '42', '42', '100.00', '1000', '1000', '1', '1', 'ref', 'contig1'],
            ['100', '142', '200', '242', '42', '42', '99.42', '1000', '1000', '1', '1', 'ref', 'contig1'],
            ['100', '110', '200', '210', '11', '11', '99.42', '1000', '1000', '1', '1', 'ref', 'contig2'],
        ]
        nucmer_hits = {
            'contig1': [
                pymummer.alignment.Alignment('\t'.join(hits[0])),
                pymummer.alignment.Alignment('\t'.join(hits[1])),
            ],
            'contig2': [
                pymummer.alignment.Alignment('\t'.join(hits[2])),
            ]
        }

        expected = {'contig1': 85, 'contig2': 11}
        got = assembly_compare.AssemblyCompare.ref_cov_per_contig(nucmer_hits)
        self.assertEqual(expected, got)


    def test_get_assembled_reference_sequences(self):
        '''test _get_assembled_reference_sequences'''
        ref_sequence = pyfastaq.sequences.Fasta('ref_seq', 'ATGGTACAAGACGGCCCTTTGCAGTCCTGTGTACTTGCGGGTCGCTCCTTTGCATTGAATTATCGAACATCGTCGCGTTCAAGATCCCGCGAAAAAAATTATAGATCGCAGGATATCACTGCCAGTGGCATCTGTGTAAGCGCTTAG')
        assembly = {
            'contig1': pyfastaq.sequences.Fasta('contig1', 'CATCTATGCTGCATCGATCACTGACGTATCATCATCAGCGTACTGACGTATTAGTTTGTAATGGTACAAGACGGCCCTTTGCAGTCCTGTGTACTTGCGGGTCGCTCCTTTGCATTGAATTATCGAACATCGTCGCGTTCAAGATCCCGCGAAAAAAATTATAGATCGCAGGATATCACTGCCAGTGGCATCTGTGTAAGCGCTTAGACGTCGTACTACTGTATATGCATCGATCTGAA'),
            'contig2': pyfastaq.sequences.Fasta('contig2', 'AGTGATATCCTGCGATCTATAATTTTTTTCGCGGGATCTTGAACGCGACGATGTTCGATAATTCAATGCAAAGGAGCGACCCGCAAGTACACAGGACTGCAAA')
        }

        hits = [
            ['1', '147', '61', '207', '147', '147', '100.00', '147', '239', '1', '1', 'ref_seq', 'contig1'],
            ['18', '120', '103', '1', '103', '103', '100.00', '147', '103', '1', '-1', 'ref_seq', 'contig2']
        ]
        nucmer_hits = {
            'contig1': [
                pymummer.alignment.Alignment('\t'.join(hits[0])),
            ],
            'contig2': [
                pymummer.alignment.Alignment('\t'.join(hits[1])),
            ]
        }

        expected = {'ref_seq.1.147.contig1.61.207.+.complete': pyfastaq.sequences.Fasta('ref_seq.1.147.contig1.61.207.+.complete', 'ATGGTACAAGACGGCCCTTTGCAGTCCTGTGTACTTGCGGGTCGCTCCTTTGCATTGAATTATCGAACATCGTCGCGTTCAAGATCCCGCGAAAAAAATTATAGATCGCAGGATATCACTGCCAGTGGCATCTGTGTAAGCGCTTAG'),
            'ref_seq.18.120.contig2.1.103.-': pyfastaq.sequences.Fasta('ref_seq.18.120.contig2.1.103.-', 'TTTGCAGTCCTGTGTACTTGCGGGTCGCTCCTTTGCATTGAATTATCGAACATCGTCGCGTTCAAGATCCCGCGAAAAAAATTATAGATCGCAGGATATCACT')
        }

        got = assembly_compare.AssemblyCompare._get_assembled_reference_sequences(nucmer_hits, ref_sequence, assembly)
        self.assertEqual(expected, got)


    def test_whole_gene_covered_by_nucmer_hits(self):
        '''test _whole_gene_covered_by_nucmer_hits'''
        ref_seq = pyfastaq.sequences.Fasta('ref', 'ACGTGTGCAT')
        hit1 = ['1', '10', '1', '10', '10', '10', '100.00', '10', '10', '1', '1', 'ref', 'contig1']
        hit2 = ['1', '5', '1', '5', '5', '5', '100.00', '10', '10', '1', '1', 'ref', 'contig2']
        hit3 = ['6', '10', '6', '10', '5', '5', '100.00', '10', '10', '1', '1', 'ref', 'contig2']
        nucmer_hits = [
            {'contig1': [pymummer.alignment.Alignment('\t'.join(hit1))]},
            {'contig2': [pymummer.alignment.Alignment('\t'.join(hit2))]},
            {'contig2': [pymummer.alignment.Alignment('\t'.join(hit2)), pymummer.alignment.Alignment('\t'.join(hit3))]}
        ]
        expected = [True, False, True]
        for i in range(len(nucmer_hits)):
            got = assembly_compare.AssemblyCompare._whole_gene_covered_by_nucmer_hits(nucmer_hits[i], ref_seq, 0.95)
            self.assertEqual(expected[i], got)


    def test_ref_has_region_assembled_twice(self):
        '''test _ref_has_region_assembled_twice'''
        ref_seq = pyfastaq.sequences.Fasta('gene', 'ACGTGTGCAT')
        hit1 = ['1', '10', '1', '10', '10', '10', '100.00', '10', '10', '1', '1', 'gene', 'contig1']
        hit2 = ['1', '5', '1', '5', '5', '5', '100.00', '10', '10', '1', '1', 'gene', 'contig2']
        nucmer_hits = { 'contig1': [pymummer.alignment.Alignment('\t'.join(hit1))] }
        self.assertFalse(assembly_compare.AssemblyCompare._ref_has_region_assembled_twice(nucmer_hits, ref_seq, 0.03))
        nucmer_hits['contig2'] = [pymummer.alignment.Alignment('\t'.join(hit2))]
        self.assertTrue(assembly_compare.AssemblyCompare._ref_has_region_assembled_twice(nucmer_hits, ref_seq, 0.03))


    def test_longest_nucmer_hit_in_ref(self):
        '''test _longest_nucmer_hit_in_ref'''
        hits = [
            ['1', '39', '1', '39', '39', '39', '100.00', '39', '39', '1', '1', 'gene', 'contig1'],
            ['1', '20', '1', '20', '20', '20', '100.00', '39', '39', '1', '1', 'gene', 'contig1'],
            ['21', '39', '21', '39', '19', '19', '100.00', '39', '39', '1', '1', 'gene', 'contig2'],
        ]
        alignments = [pymummer.alignment.Alignment('\t'.join(x)) for x in hits]
        nucmer_hits = {
            'contig1': [alignments[0]],
            'contig2': [alignments[1], alignments[2]],
        }
        got = assembly_compare.AssemblyCompare._longest_nucmer_hit_in_ref(nucmer_hits)
        self.assertEqual(alignments[0], got)


    def test_find_previous_start_codon(self):
        '''test _find_previous_start_codon'''
        tests = [
            ('ATGTTTAAA', 0, 10, 0),
            ('TATGTTTAAA', 0, 10, None),
            ('TATGTTTAAA', 1, 10, 1),
            ('ATATGTTTAAA', 2, 10, 2),
            ('AATGTTTAAA', 7, 1, None),
            ('AATGTTTAAA', 7, 2, None),
            ('AATGTTTAAA', 7, 3, None),
            ('AATGTTTAAA', 7, 4, None),
            ('AATGTTTAAA', 7, 5, None),
            ('AATGTTTAAA', 7, 6, 1),
            ('AATGTTTAAA', 7, 7, 1),
            ('AGTGTTTAAA', 7, 7, None),
        ]

        for seq, start_coord, max_nt_to_extend, expected in tests:
            fa = pyfastaq.sequences.Fasta('x', seq)
            got = assembly_compare.AssemblyCompare._find_previous_start_codon(fa, start_coord, max_nt_to_extend)
            self.assertEqual(expected, got)


    def test_ref_covered_by_complete_contig_with_orf(self):
        '''test _ref_covered_by_complete_contig_with_orf'''
        gene = pyfastaq.sequences.Fasta('gene', 'GATCGCGAAGCGATGACCCATGAAGCGACCGAACGCTGA')
        gene_no_orf = pyfastaq.sequences.Fasta('gene', 'GATTGAGAAGCGATGACCCATGAAGCGACCGAACGCTGA')
        hit1 = ['1', '39', '1', '39', '39', '39', '100.00', '39', '39', '1', '1', 'gene', 'contig1']
        hit2 = ['1', '20', '1', '20', '20', '20', '100.00', '39', '39', '1', '1', 'gene', 'contig1']
        hit3 = ['21', '39', '21', '39', '19', '19', '100.00', '39', '39', '1', '1', 'gene', 'contig2']
        nucmer_hits = [
            {'contig1': [pymummer.alignment.Alignment('\t'.join(hit1))]},
            {'contig1': [pymummer.alignment.Alignment('\t'.join(hit1))]},
            {'contig2': [pymummer.alignment.Alignment('\t'.join(hit2))]},
            {'contig2': [pymummer.alignment.Alignment('\t'.join(hit2)), pymummer.alignment.Alignment('\t'.join(hit3))]},
        ]
        expected = [True, False, False, False]
        assemblies = [
            {'contig1': gene},
            {'contig1': gene_no_orf},
            {'contig1': gene},
            {'contig1': gene, 'contig2': pyfastaq.sequences.Fasta('contig2', 'ACGT')}
        ]
        assert len(expected) == len(nucmer_hits) == len(assemblies)
        for i in range(len(expected)):
            self.assertEqual(expected[i], assembly_compare.AssemblyCompare._ref_covered_by_complete_contig_with_orf(nucmer_hits[i], assemblies[i]))


    def test_ref_covered_by_at_least_one_full_length_contig(self):
        '''test _ref_covered_by_at_least_one_full_length_contig'''
        ref = pyfastaq.sequences.Fasta('gene', 'GATCGCGAAGCGATGACCCATGAAGCGACCGAACGCTGA')
        hit1 = ['1', '39', '1', '39', '39', '39', '100.00', '39', '39', '1', '1', 'ref', 'contig1']
        hit2 = ['1', '20', '1', '20', '20', '20', '100.00', '39', '39', '1', '1', 'ref', 'contig1']
        nucmer_hits = [
            {'contig1': [pymummer.alignment.Alignment('\t'.join(hit1))]},
            {'contig1': [pymummer.alignment.Alignment('\t'.join(hit2))]},
        ]
        expected = [True, False]
        assert len(expected) == len(nucmer_hits)
        for i in range(len(expected)):
            self.assertEqual(expected[i], assembly_compare.AssemblyCompare._ref_covered_by_at_least_one_full_length_contig(nucmer_hits[i]))

    def test_nucmer_hit_containing_reference_position(self):
        '''test nucmer_hit_containing_reference_position'''
        listhit1 = ['100', '200', '300', '400', '100', '100', '100.00', '600', '500', '1', '1', 'ref', 'contig1']
        listhit2 = ['400', '500', '500', '600', '100', '100', '100.00', '600', '600', '1', '1', 'ref', 'contig2']
        hit1 = pymummer.alignment.Alignment('\t'.join(listhit1))
        hit2 = pymummer.alignment.Alignment('\t'.join(listhit2))
        nucmer_hits = {
            'contig1': [hit1],
            'contig2': [hit2],
        }

        tests = [
            ('ref2', 150, None),
            ('ref', 42, None),
            ('ref', 98, None),
            ('ref', 200, None),
            ('ref', 99, hit1),
            ('ref', 142, hit1),
            ('ref', 199, hit1),
            ('ref', 200, None),
            ('ref', 398, None),
            ('ref', 399, hit2),
            ('ref', 442, hit2),
            ('ref', 499, hit2),
            ('ref', 500, None),
        ]

        for ref_name, ref_pos, expected in tests:
            got = assembly_compare.AssemblyCompare.nucmer_hit_containing_reference_position(nucmer_hits, ref_name, ref_pos)
            self.assertEqual(expected, got)
