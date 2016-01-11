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


class TestAssemblyCompare(unittest.TestCase):
    def test_parse_nucmer_coords_file(self):
        '''test _parse_nucmer_coords_file'''
        cluster_dir = os.path.join(data_dir, 'cluster_test_parse_assembly_vs_gene_coords')
        coords_file = os.path.join(data_dir, 'cluster_test_parse_assembly_vs_gene_coords.coords')
        clean_cluster_dir(cluster_dir)
        c = cluster.Cluster(cluster_dir, 'name')
        shutil.copyfile(coords_file, c.assembly_vs_gene_coords)
        c.gene = pyfastaq.sequences.Fasta('gene', 'AAACCCGGGTTT')
        c._parse_assembly_vs_gene_coords()
        line1 = ['1', '1000', '1', '1000', '1000', '1000', '100.00', '1000', '1000', '1', '1', 'gene', 'contig1']
        line2 = ['1', '240', '1', '240', '240', '240', '100.00', '1000', '580', '1', '1', 'gene', 'contig2']
        line3 = ['661', '1000', '241', '580', '340', '340', '100.00', '1000', '580', '1', '1', 'gene', 'contig2']
        expected = {
            'contig1': [pymummer.alignment.Alignment('\t'.join(line1))],
            'contig2': [pymummer.alignment.Alignment('\t'.join(line2)), pymummer.alignment.Alignment('\t'.join(line3))],
        }
        self.assertEqual(expected, c.nucmer_hits)
        clean_cluster_dir(cluster_dir)


    def test_nucmer_hits_to_percent_identity(self):
        '''test _nucmer_hits_to_percent_identity'''
        cluster_dir = os.path.join(data_dir, 'cluster_test_generic')
        clean_cluster_dir(cluster_dir)
        c = cluster.Cluster(cluster_dir, 'name')
        hits = [
            ['1', '10', '1', '10', '10', '10', '90.00', '1000', '1000', '1', '1', 'gene', 'scaff1'],
            ['9', '42', '9', '42', '34', '34', '100.00', '1000', '1000', '1', '1', 'gene', 'scaff1'],
            ['1', '42', '1', '42', '42', '42', '42.42', '1000', '1000', '1', '1', 'gene', 'scaff2'],
        ]
        c.nucmer_hits = {
            'scaff1': [
                pymummer.alignment.Alignment('\t'.join(hits[0])),
                pymummer.alignment.Alignment('\t'.join(hits[1])),
            ],
            'scaff2': [
                pymummer.alignment.Alignment('\t'.join(hits[2])),
            ]
        }
        expected = {'scaff1': round((90*10 + 100*34) / (10+34), 2), 'scaff2': 42.42}
        c._nucmer_hits_to_percent_identity()
        self.assertEqual(expected, c.percent_identities)


    def test_nucmer_hits_to_assembly_coords(self):
        '''test _nucmer_hits_to_assembly_coords'''
        cluster_dir = os.path.join(data_dir, 'cluster_test_generic')
        clean_cluster_dir(cluster_dir)
        c = cluster.Cluster(cluster_dir, 'name')
        hits = [
            ['1', '10', '1', '10', '10', '10', '100.00', '1000', '1000', '1', '1', 'gene', 'scaff1'],
            ['9', '42', '9', '42', '34', '34', '100.00', '1000', '1000', '1', '1', 'gene', 'scaff1'],
            ['50', '52', '50', '52', '3', '3', '100.00', '1000', '1000', '1', '1', 'gene', 'scaff1'],
            ['1', '42', '1', '42', '42', '42', '100.00', '1000', '1000', '1', '1', 'gene', 'scaff2'],
        ]
        c.nucmer_hits = {
            'scaff1': [
                pymummer.alignment.Alignment('\t'.join(hits[0])),
                pymummer.alignment.Alignment('\t'.join(hits[1])),
                pymummer.alignment.Alignment('\t'.join(hits[2])),
            ],
            'scaff2': [
                pymummer.alignment.Alignment('\t'.join(hits[3])),
            ]
        }
        got = c._nucmer_hits_to_scaff_coords()
        expected = {
            'scaff1': [
                pyfastaq.intervals.Interval(0, 41),
                pyfastaq.intervals.Interval(49, 51)
            ],
            'scaff2': [
                pyfastaq.intervals.Interval(0, 41),
            ]
        }
        self.assertEqual(got, expected)
        clean_cluster_dir(cluster_dir)


    def test_nucmer_hits_to_ref_coords(self):
        '''test _nucmer_hits_to_ref_coords'''
        cluster_dir = os.path.join(data_dir, 'cluster_test_generic')
        clean_cluster_dir(cluster_dir)
        c = cluster.Cluster(cluster_dir, 'name')
        hits = [
            ['1', '42', '1', '42', '42', '42', '100.00', '1000', '1000', '1', '1', 'gene', 'contig1'],
            ['100', '142', '200', '242', '42', '42', '99.42', '1000', '1000', '1', '1', 'gene', 'contig1'],
            ['100', '110', '200', '210', '11', '11', '99.42', '1000', '1000', '1', '1', 'gene', 'contig2'],
        ]
        c.nucmer_hits = {
            'contig1': [
                pymummer.alignment.Alignment('\t'.join(hits[0])),
                pymummer.alignment.Alignment('\t'.join(hits[1])),
            ],
            'contig2': [
                pymummer.alignment.Alignment('\t'.join(hits[2])),
            ]
        }
        got_coords = c._nucmer_hits_to_ref_coords()
        expected = [
            pyfastaq.intervals.Interval(0,41),
            pyfastaq.intervals.Interval(99, 109),
            pyfastaq.intervals.Interval(99, 141),
        ]
        self.assertEqual(got_coords, expected)

        got_coords = c._nucmer_hits_to_ref_coords(contig='contig2')
        expected = [
            pyfastaq.intervals.Interval(99, 109),
        ]
        self.assertEqual(got_coords, expected)
        clean_cluster_dir(cluster_dir)


    def test_ref_cov_per_contig(self):
        '''test __ref_cov_per_contig'''
        cluster_dir = os.path.join(data_dir, 'cluster_test_generic')
        clean_cluster_dir(cluster_dir)
        c = cluster.Cluster(cluster_dir, 'name')
        hits = [
            ['1', '42', '1', '42', '42', '42', '100.00', '1000', '1000', '1', '1', 'gene', 'contig1'],
            ['100', '142', '200', '242', '42', '42', '99.42', '1000', '1000', '1', '1', 'gene', 'contig1'],
            ['100', '110', '200', '210', '11', '11', '99.42', '1000', '1000', '1', '1', 'gene', 'contig2'],
        ]
        c.nucmer_hits = {
            'contig1': [
                pymummer.alignment.Alignment('\t'.join(hits[0])),
                pymummer.alignment.Alignment('\t'.join(hits[1])),
            ],
            'contig2': [
                pymummer.alignment.Alignment('\t'.join(hits[2])),
            ]
        }

        expected = {'contig1': 85, 'contig2': 11}
        self.assertEqual(expected, c._nucmer_hits_to_gene_cov_per_contig())
        clean_cluster_dir(cluster_dir)


    def test_write_assembled_reference_sequences(self):
        '''test _write_assembled_reference_sequences'''
        ref_gene = pyfastaq.sequences.Fasta('ref_gene', 'ATGGTACAAGACGGCCCTTTGCAGTCCTGTGTACTTGCGGGTCGCTCCTTTGCATTGAATTATCGAACATCGTCGCGTTCAAGATCCCGCGAAAAAAATTATAGATCGCAGGATATCACTGCCAGTGGCATCTGTGTAAGCGCTTAG')
        assembly = {
            'contig1': pyfastaq.sequences.Fasta('contig1', 'CATCTATGCTGCATCGATCACTGACGTATCATCATCAGCGTACTGACGTATTAGTTTGTAATGGTACAAGACGGCCCTTTGCAGTCCTGTGTACTTGCGGGTCGCTCCTTTGCATTGAATTATCGAACATCGTCGCGTTCAAGATCCCGCGAAAAAAATTATAGATCGCAGGATATCACTGCCAGTGGCATCTGTGTAAGCGCTTAGACGTCGTACTACTGTATATGCATCGATCTGAA'),
            'contig2': pyfastaq.sequences.Fasta('contig2', 'AGTGATATCCTGCGATCTATAATTTTTTTCGCGGGATCTTGAACGCGACGATGTTCGATAATTCAATGCAAAGGAGCGACCCGCAAGTACACAGGACTGCAAA')
        }

        hits = [
            ['1', '147', '61', '207', '147', '147', '100.00', '147', '239', '1', '1', 'ref_gene', 'contig1'],
            ['18', '120', '103', '1', '103', '103', '100.00', '147', '103', '1', '-1', 'ref_gene', 'contig2']
        ]
        nucmer_hits = {
            'contig1': [
                pymummer.alignment.Alignment('\t'.join(hits[0])),
            ],
            'contig2': [
                pymummer.alignment.Alignment('\t'.join(hits[1])),
            ]
        }

        assembly_fasta = os.path.join(data_dir, 'cluster_test_nucmer_hits_to_assembled_gene_sequences.assembly.fa')
        tmp_outfile = 'tmp.test_nucmer_hits_to_assembled_gene_sequences.out.fa'
        expected_outfile = os.path.join(data_dir, 'cluster_test_nucmer_hits_to_assembled_gene_sequences.expected.out.fa')
        cluster.Cluster._nucmer_hits_to_assembled_gene_sequences(nucmer_hits, ref_gene, assembly, tmp_outfile)
        self.assertTrue(filecmp.cmp(tmp_outfile, expected_outfile, shallow=False))
        os.unlink(tmp_outfile)


    def test_whole_gene_covered_by_nucmer_hits(self):
        '''test _whole_gene_covered_by_nucmer_hits'''
        cluster_dir = os.path.join(data_dir, 'cluster_test_generic')
        clean_cluster_dir(cluster_dir)
        c = cluster.Cluster(cluster_dir, 'name')
        c.gene = pyfastaq.sequences.Fasta('gene', 'ACGTGTGCAT')
        hit1 = ['1', '10', '1', '10', '10', '10', '100.00', '10', '10', '1', '1', 'gene', 'contig1']
        hit2 = ['1', '5', '1', '5', '5', '5', '100.00', '10', '10', '1', '1', 'gene', 'contig2']
        hit3 = ['6', '10', '6', '10', '5', '5', '100.00', '10', '10', '1', '1', 'gene', 'contig2']
        nucmer_hits = [
            {'contig1': [pymummer.alignment.Alignment('\t'.join(hit1))]},
            {'contig2': [pymummer.alignment.Alignment('\t'.join(hit2))]},
            {'contig2': [pymummer.alignment.Alignment('\t'.join(hit2)), pymummer.alignment.Alignment('\t'.join(hit3))]}
        ]
        expected = [True, False, True]
        for i in range(len(nucmer_hits)):
            c.nucmer_hits = nucmer_hits[i]
            self.assertEqual(expected[i], c._whole_gene_covered_by_nucmer_hits())

        clean_cluster_dir(cluster_dir)


    def test_ref_has_region_assembled_twice(self):
        '''test _ref_has_region_assembled_twice'''
        cluster_dir = os.path.join(data_dir, 'cluster_test_generic')
        clean_cluster_dir(cluster_dir)
        c = cluster.Cluster(cluster_dir, 'name')
        c.gene = pyfastaq.sequences.Fasta('gene', 'ACGTGTGCAT')
        hit1 = ['1', '10', '1', '10', '10', '10', '100.00', '10', '10', '1', '1', 'gene', 'contig1']
        hit2 = ['1', '5', '1', '5', '5', '5', '100.00', '10', '10', '1', '1', 'gene', 'contig2']
        c.nucmer_hits = { 'contig1': [pymummer.alignment.Alignment('\t'.join(hit1))] }
        self.assertTrue(c._gene_coverage_unique())
        c.nucmer_hits['contig2'] = [pymummer.alignment.Alignment('\t'.join(hit2))]
        self.assertFalse(c._gene_coverage_unique())


    def test_ref_covered_by_complete_contig_with_orf(self):
        '''test _ref_covered_by_complete_contig_with_orf'''
        cluster_dir = os.path.join(data_dir, 'cluster_test_generic')
        clean_cluster_dir(cluster_dir)
        c = cluster.Cluster(cluster_dir, 'name')
        gene = pyfastaq.sequences.Fasta('gene', 'GATCGCGAAGCGATGACCCATGAAGCGACCGAACGCTGA')
        gene_no_orf = pyfastaq.sequences.Fasta('gene', 'GATTGAGAAGCGATGACCCATGAAGCGACCGAACGCTGA')
        c.gene = gene
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
            c.final_assembly = assemblies[i]
            c.nucmer_hits = nucmer_hits[i]
            self.assertEqual(c._gene_covered_by_complete_contig_with_orf(), expected[i])
        clean_cluster_dir(cluster_dir)


    def test_ref_covered_by_at_least_one_full_length_contig(self):
        '''test _ref_covered_by_at_least_one_full_length_contig'''
        cluster_dir = os.path.join(data_dir, 'cluster_test_generic')
        clean_cluster_dir(cluster_dir)
        c = cluster.Cluster(cluster_dir, 'name')
        c.gene = pyfastaq.sequences.Fasta('gene', 'GATCGCGAAGCGATGACCCATGAAGCGACCGAACGCTGA')
        hit1 = ['1', '39', '1', '39', '39', '39', '100.00', '39', '39', '1', '1', 'gene', 'contig1']
        hit2 = ['1', '20', '1', '20', '20', '20', '100.00', '39', '39', '1', '1', 'gene', 'contig1']
        nucmer_hits = [
            {'contig1': [pymummer.alignment.Alignment('\t'.join(hit1))]},
            {'contig1': [pymummer.alignment.Alignment('\t'.join(hit2))]},
        ]
        expected = [True, False]
        assert len(expected) == len(nucmer_hits)
        for i in range(len(expected)):
            c.nucmer_hits = nucmer_hits[i]
            self.assertEqual(c._gene_covered_by_at_least_one_full_length_contig(), expected[i])
        clean_cluster_dir(cluster_dir)

