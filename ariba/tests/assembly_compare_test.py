import unittest
import os
import pyfastaq
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
        expected = {'scaff1': 100.0, 'scaff2': 42.42}
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


    def test_nucmer_hits_to_ref_and_qry_coords(self):
        '''test _nucmer_hits_to_ref_and_qry_coords'''
        hits = [
            ['1', '42', '1', '42', '42', '42', '100.00', '1000', '1000', '1', '1', 'ref', 'contig1'],
            ['31', '52', '1', '22', '22', '22', '100.00', '1000', '1000', '1', '1', 'ref', 'contig1'],
            ['11', '32', '1000', '1022', '22', '22', '100.00', '1000', '1000', '1', '1', 'ref', 'contig1'],
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

        got_ctg, got_ref = assembly_compare.AssemblyCompare.nucmer_hits_to_ref_and_qry_coords(nucmer_hits)
        expected_ctg = {
            'contig1': [pyfastaq.intervals.Interval(0,51), pyfastaq.intervals.Interval(99, 141)],
            'contig2': [pyfastaq.intervals.Interval(99, 109)]
        }
        expected_ref = {
            'contig1': [pyfastaq.intervals.Interval(0,51), pyfastaq.intervals.Interval(99, 141)],
            'contig2': [pyfastaq.intervals.Interval(99, 109)]
        }


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
            got = assembly_compare.AssemblyCompare._whole_gene_covered_by_nucmer_hits(nucmer_hits[i], ref_seq, 0.95, 0)
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
            ('ATGTTTAAA', 0, 0, 0),
            ('TATGTTTAAA', 0, 0, None),
            ('TATGTTTAAA', 1, 0, 1),
            ('ATATGTTTAAA', 2, 0, 2),
            ('AATGTTTAAA', 7, 6, None),
            ('AATGTTTAAA', 7, 5, None),
            ('AATGTTTAAA', 7, 4, None),
            ('AATGTTTAAA', 7, 3, None),
            ('AATGTTTAAA', 7, 2, None),
            ('AATGTTTAAA', 7, 1, 1),
            ('AATGTTTAAA', 7, 0, 1),
            ('AGTGTTTAAA', 7, 0, None),
            ('AATGTAGAAA', 7, 0, None),
        ]

        for seq, start_coord, min_coord, expected in tests:
            fa = pyfastaq.sequences.Fasta('x', seq)
            got = assembly_compare.AssemblyCompare._find_previous_start_codon(fa, start_coord, min_coord)
            self.assertEqual(expected, got)


    def test_find_next_stop_codon(self):
        '''test _find_next_stop_codon'''
        tests = [
            ('ATGTTTAGA', 2, 5, None),
            ('ATGTTTAGA', 2, 6, None),
            ('ATGTTTAGA', 2, 7, 5),
            ('ATGTTTGGA', 2, 4, None),
            ('ATGTTTAGA', 5, 6, None),
            ('ATGTTTAGA', 5, 7, 5),
            ('ATGTTTGGA', 5, 7, None),
            ('ATGTTTAGA', 4, 7, None),
            ('ATGTTTAGA', 3, 7, None),
        ]

        for seq, end_coord, max_coord, expected in tests:
            fa = pyfastaq.sequences.Fasta('x', seq)
            got = assembly_compare.AssemblyCompare._find_next_stop_codon(fa, end_coord, max_coord)
            self.assertEqual(expected, got)


    def test_gene_from_nucmer_match(self):
        '''test _gene_from_nucmer_match'''
        tests = [
            (
             ['1', '15', '2', '16', '11', '11', '100.00', '20', '20', '1', '1', 'ref', 'contig'], 'AATGTAGTTTCCCTAGATAT', 5,
             (pyfastaq.sequences.Fasta('contig.2-16', 'ATGTAGTTTCCCTAG'), 'HAS_STOP', None, None)
            ),
            (
             ['1', '15', '2', '16', '11', '11', '100.00', '20', '20', '1', '1', 'ref', 'contig'], 'AATGAAATTTCCCTAGATAT', 5,
             (pyfastaq.sequences.Fasta('contig.2-16', 'ATGAAATTTCCCTAG'), 'GENE_FOUND', 0, 0)
            ),
            (
             ['1', '15', '2', '15', '11', '11', '100.00', '20', '20', '1', '1', 'ref', 'contig'], 'AATGAAATTTCCCTAGATAT', 5,
             (pyfastaq.sequences.Fasta('contig.2-16', 'ATGAAATTTCCCTAG'), 'GENE_FOUND', 0, 1)
            ),
            (
             ['1', '15', '2', '14', '11', '11', '100.00', '20', '20', '1', '1', 'ref', 'contig'], 'AATGAAATTTCCCTAGATAT', 5,
             (pyfastaq.sequences.Fasta('contig.2-16', 'ATGAAATTTCCCTAG'), 'GENE_FOUND', 0, 2)
            ),
            (
             ['1', '15', '2', '13', '11', '11', '100.00', '20', '20', '1', '1', 'ref', 'contig'], 'AATGAAATTTCCCTAGATAT', 5,
             (pyfastaq.sequences.Fasta('contig.2-16', 'ATGAAATTTCCCTAG'), 'GENE_FOUND', 0, 3)
            ),
            (
             ['2', '15', '3', '16', '11', '11', '100.00', '20', '20', '1', '1', 'ref', 'contig'], 'AATGAAATTTCCCTAGATAT', 5,
             (pyfastaq.sequences.Fasta('contig.2-16', 'ATGAAATTTCCCTAG'), 'GENE_FOUND', 1, 0)
            ),
            (
             ['3', '15', '4', '16', '11', '11', '100.00', '20', '20', '1', '1', 'ref', 'contig'], 'AATGAAATTTCCCTAGATAT', 5,
             (pyfastaq.sequences.Fasta('contig.2-16', 'ATGAAATTTCCCTAG'), 'GENE_FOUND', 2, 0)
            ),
            (
             ['4', '15', '5', '16', '11', '11', '100.00', '20', '20', '1', '1', 'ref', 'contig'], 'AATGAAATTTCCCTAGATAT', 5,
             (pyfastaq.sequences.Fasta('contig.2-16', 'ATGAAATTTCCCTAG'), 'GENE_FOUND', 3, 0)
            ),
            (
             ['1', '15', '2', '14', '11', '11', '100.00', '20', '20', '1', '1', 'ref', 'contig'], 'AATGAAATTTCCCTAGATAT', 0,
             (pyfastaq.sequences.Fasta('contig.2-13', 'ATGAAATTTCCC'), 'START_OR_END_FAIL', 0, None)
            ),
            (
             ['1', '15', '2', '14', '11', '11', '100.00', '20', '20', '1', '1', 'ref', 'contig'], 'AATGAAATTTCCCTAGATAT', 1,
             (pyfastaq.sequences.Fasta('contig.2-13', 'ATGAAATTTCCC'), 'START_OR_END_FAIL', 0, None)
            ),
            (
             ['1', '15', '2', '14', '11', '11', '100.00', '20', '20', '1', '1', 'ref', 'contig'], 'AATGAAATTTCCCTAGATAT', 2,
             (pyfastaq.sequences.Fasta('contig.2-16', 'ATGAAATTTCCCTAG'), 'GENE_FOUND', 0, 2)
            ),
            (
             ['2', '15', '3', '14', '11', '11', '100.00', '20', '20', '1', '1', 'ref', 'contig'], 'AATGAAATTTCCCTAGATAT', 2,
             (pyfastaq.sequences.Fasta('contig.2-16', 'ATGAAATTTCCCTAG'), 'GENE_FOUND', 1, 2)
            ),
            (
             ['2', '15', '18', '7', '11', '11', '100.00', '20', '20', '1', '1', 'ref', 'contig'], 'ATATCTAGGGAAATTTCATT', 2,
             (pyfastaq.sequences.Fasta('contig.16-2', 'ATGAAATTTCCCTAG'), 'GENE_FOUND', 1, 2)
            ),
        ]

        for hit, seq, max_extend, expected in tests:
            nucmer_match = pymummer.alignment.Alignment('\t'.join(hit))
            contig = pyfastaq.sequences.Fasta('contig', seq)
            got = assembly_compare.AssemblyCompare._gene_from_nucmer_match(nucmer_match, contig, max_extend)
            self.assertEqual(expected, got)


    def test_get_gene_matching_ref(self):
        '''test _get_gene_matching_ref'''
        hit1 = ['2', '15', '3', '14', '11', '11', '100.00', '20', '20', '1', '1', 'ref', 'contig']
        hit2 = ['2', '7', '3', '8', '6', '6', '100.00', '20', '20', '1', '1', 'ref', 'contig2']
        contigs = {
            'contig': pyfastaq.sequences.Fasta('contig', 'AATGAAATTTCCCTAGATAT'),
            'contig2': pyfastaq.sequences.Fasta('contig2', 'AATGAAATTTCCCTAGATAT')
        }
        nucmer_hits = {
            'contig': [pymummer.alignment.Alignment('\t'.join(hit1))],
            'contig2': [pymummer.alignment.Alignment('\t'.join(hit2))],
        }

        got = assembly_compare.AssemblyCompare._get_gene_matching_ref(nucmer_hits, contigs, 10)
        expected = ('contig', pyfastaq.sequences.Fasta('contig.2-16', 'ATGAAATTTCCCTAG'), 'GENE_FOUND', 1, 2)
        self.assertEqual(expected, got)


    def test_ref_covered_by_at_least_one_full_length_contig(self):
        '''test _ref_covered_by_at_least_one_full_length_contig'''
        hits = [
            ['1', '100', '1', '100', '100', '100', '100.00', '100', '100', '1', '1', 'ref', 'contig1'],
            ['1', '99', '1', '99', '99', '99', '100.00', '100', '100', '1', '1', 'ref', 'contig1'],
            ['1', '96', '1', '96', '96', '96', '100.00', '100', '100', '1', '1', 'ref', 'contig1'],
            ['1', '95', '1', '95', '95', '95', '100.00', '100', '100', '1', '1', 'ref', 'contig1'],
            ['1', '94', '1', '94', '94', '94', '100.00', '100', '100', '1', '1', 'ref', 'contig1'],
        ]
        nucmer_hits = [{'contig1': [pymummer.alignment.Alignment('\t'.join(hit))]} for hit in hits]
        expected = [True, True, True, True, False]
        assert len(expected) == len(nucmer_hits)
        for i in range(len(expected)):
            self.assertEqual(expected[i], assembly_compare.AssemblyCompare._ref_covered_by_at_least_one_full_length_contig(nucmer_hits[i], 0.95, 0))


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
