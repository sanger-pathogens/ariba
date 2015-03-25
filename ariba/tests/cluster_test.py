import unittest
import os
import copy
import shutil
import filecmp
import pyfastaq
import pysam
import pymummer
from ariba import cluster, flag

modules_dir = os.path.dirname(os.path.abspath(cluster.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


def clean_cluster_dir(d, exclude=None):
    if not os.path.exists(d):
        return

    '''Cleans up all files made except original ones in a cluster directory'''
    keep = set(['genes.fa', 'reads_1.fq', 'reads_2.fq'])
    if exclude is not None:
        for f in exclude:
            keep.add(f)

    for name in os.listdir(d):
        if name not in keep:
            full_path = os.path.join(d, name)
            if os.path.isdir(full_path):
                shutil.rmtree(full_path)
            else:
                os.unlink(full_path)


def load_gene(filename):
    file_reader = pyfastaq.sequences.file_reader(filename)
    seq = None
    for seq in file_reader:
        pass
    return seq


class TestCluster(unittest.TestCase):
    def test_init_fail_files_missing(self):
        '''test init_fail_files_missing'''
        dirs = [
            'cluster_test_directorynotexist'
            'cluster_test_init_no_genes_fa',
            'cluster_test_init_no_reads_1',
            'cluster_test_init_no_reads_2',
        ]
        dirs = [os.path.join(data_dir, d) for d in dirs]
        for d in dirs:
            clean_cluster_dir(d)
            with self.assertRaises(cluster.Error):
                c = cluster.Cluster(d, 'name')
            clean_cluster_dir(d)


    def test_get_read_counts(self):
        '''test _get_read_counts pass'''
        cluster_dir = os.path.join(data_dir, 'cluster_test_get_read_counts')
        clean_cluster_dir(cluster_dir)
        c = cluster.Cluster(cluster_dir, 'name')
        self.assertEqual(2, c._get_read_counts())
        clean_cluster_dir(cluster_dir)


    def test_get_read_counts_fail(self):
        '''test _get_read_counts fail'''
        cluster_dir = os.path.join(data_dir, 'cluster_test_get_read_counts_fail')
        clean_cluster_dir(cluster_dir)
        c = cluster.Cluster(cluster_dir, 'name')
        with self.assertRaises(cluster.Error):
            c._get_read_counts()
        clean_cluster_dir(cluster_dir)


    def test_get_total_alignment_score(self):
        '''test _get_total_alignment_score'''
        cluster_dir = os.path.join(data_dir, 'cluster_test_get_total_alignment_score')
        clean_cluster_dir(cluster_dir)
        c = cluster.Cluster(cluster_dir, 'name')
        got_score = c._get_total_alignment_score('1')
        expected_score = 1500
        self.assertEqual(got_score, expected_score)
        clean_cluster_dir(cluster_dir)


    def test_get_best_gene_by_alignment_score(self):
        '''test _get_best_gene_by_alignment_score'''
        cluster_dir = os.path.join(data_dir, 'cluster_test_get_best_gene_by_alignment_score')
        clean_cluster_dir(cluster_dir)
        c = cluster.Cluster(cluster_dir, 'name')
        got_name = c._get_best_gene_by_alignment_score()
        self.assertEqual(got_name, '1')
        clean_cluster_dir(cluster_dir)


    def test_choose_best_gene(self):
        '''test _choose_best_gene'''
        cluster_dir = os.path.join(data_dir, 'cluster_test_choose_best_gene')
        clean_cluster_dir(cluster_dir)
        c = cluster.Cluster(cluster_dir, 'name')
        expected_gene = pyfastaq.sequences.Fasta('1', ''.join([
            'AGCGCCTAGCTTTGGCACTTCAGGAGCGCCCGGAAATAATGGCGGGCGATGAAGGTTCTG',
            'TAGGTACGCAAGATCCCTCTTAATCACAGTGGTGTAATCTGCGGGTCAGACCCTGTTAAC',
            'CCGTGGCTTTCACACTCCCTCCTATGGGTAATCAATCCAGAAAGGGGCCGAAATGCAAAA',
            'GTCTTAAGGACTCTGCGAGGCAAAGTACGGGCGAACTAAACCCCCGTGACAGGTCAGACG',
            'TTGTTTCGGCAATCTGTCGCGCTCCCACACCTATAAGCGTACACCGTCTCTTCTGCCAGC',
        ]))
        expected_gene_fa = os.path.join(data_dir, 'cluster_test_choose_best_gene.gene.fa')
        got = c._choose_best_gene()
        self.assertEqual(got, expected_gene)
        self.assertTrue(filecmp.cmp(expected_gene_fa, c.gene_fa, shallow=False))
        clean_cluster_dir(cluster_dir)


    def test_set_assembly_kmer(self):
        '''test _set_assembly_kmer'''
        cluster_dir = os.path.join(data_dir, 'cluster_test_set_assembly_kmer')
        clean_cluster_dir(cluster_dir)
        c = cluster.Cluster(cluster_dir, 'name', assembly_kmer=42)
        self.assertEqual(c.assembly_kmer, 42)
        clean_cluster_dir(cluster_dir)
        c = cluster.Cluster(os.path.join(data_dir, 'cluster_test_set_assembly_kmer'), 'name')
        self.assertEqual(c.assembly_kmer, 5)
        clean_cluster_dir(cluster_dir)


    def test_assemble_with_spades(self):
        '''test _assemble_with_spades'''
        cluster_dir = os.path.join(data_dir, 'cluster_test_assemble_with_spades')
        clean_cluster_dir(cluster_dir)
        c = cluster.Cluster(cluster_dir, 'name')
        shutil.copyfile(os.path.join(data_dir, 'cluster_test_assemble_with_spades.gene.fa'), c.gene_fa)
        c._assemble_with_spades(unittest=True)
        self.assertEqual(c.status_flag.to_number(), 0)
        clean_cluster_dir(cluster_dir)


    def test_assemble_with_spades_fail(self):
        '''test _assemble_with_spades handles spades fail'''
        cluster_dir = os.path.join(data_dir, 'cluster_test_assemble_with_spades')
        clean_cluster_dir(cluster_dir)
        c = cluster.Cluster(cluster_dir, 'name')
        shutil.copyfile(os.path.join(data_dir, 'cluster_test_assemble_with_spades.gene.fa'), c.gene_fa)
        c._assemble_with_spades()
        self.assertEqual(c.status_flag.to_number(), 64)
        clean_cluster_dir(cluster_dir)


    def test_scaffold_with_sspace(self):
        '''test _scaffold_with_sspace'''
        cluster_dir = os.path.join(data_dir, 'cluster_test_scaffold_with_sspace')
        clean_cluster_dir(cluster_dir)
        c = cluster.Cluster(cluster_dir, 'name')
        shutil.copyfile(
            os.path.join(data_dir, 'cluster_test_scaffold_with_sspace.contigs.fa'),
            c.assembly_contigs
        )
        #shutil.copyfile(os.path.join(data_dir, 'cluster_test_scaffold_with_sspace.gene.fa'), c.gene_fa)
        c._scaffold_with_sspace()
        self.assertTrue(os.path.exists(c.scaffolder_scaffolds))
        clean_cluster_dir(cluster_dir)


    def test_gap_fill_with_gapfiller_no_gaps(self):
        '''test _gap_fill_with_gapfiller no gaps'''
        cluster_dir = os.path.join(data_dir, 'cluster_test_gapfill_with_gapfiller')
        clean_cluster_dir(cluster_dir)
        c = cluster.Cluster(cluster_dir, 'name')
        shutil.copyfile(
            os.path.join(data_dir, 'cluster_test_gapfill_with_gapfiller.scaffolds_no_gaps.fa'),
            c.scaffolder_scaffolds
        )
        c.gene = pyfastaq.sequences.Fasta('name_of_gene', 'AAACCCGGGTTT')
        c._gap_fill_with_gapfiller()
        self.assertTrue(os.path.exists(c.gapfilled_scaffolds))
        clean_cluster_dir(cluster_dir)


    def test_gap_fill_with_gapfiller_with_gaps(self):
        '''test _gap_fill_with_gapfiller with gaps'''
        cluster_dir = os.path.join(data_dir, 'cluster_test_gapfill_with_gapfiller')
        clean_cluster_dir(cluster_dir)
        c = cluster.Cluster(cluster_dir, 'name')
        shutil.copyfile(
            os.path.join(data_dir, 'cluster_test_gapfill_with_gapfiller.scaffolds_with_gaps.fa'),
            c.scaffolder_scaffolds
        )
        c.gene = pyfastaq.sequences.Fasta('name_of_gene', 'AAACCCGGGTTT')
        c._gap_fill_with_gapfiller()
        self.assertTrue(os.path.exists(c.gapfilled_scaffolds))
        clean_cluster_dir(cluster_dir)


    def test_rename_scaffolds(self):
        '''test _rename_scaffolds'''
        cluster_dir = os.path.join(data_dir, 'cluster_test_rename_scaffolds')
        clean_cluster_dir(cluster_dir)
        c = cluster.Cluster(cluster_dir, 'name')
        c.gene = pyfastaq.sequences.Fasta('name_of_gene', 'AAACCCGGGTTT')
        infile = os.path.join(data_dir, 'cluster_test_rename_scaffolds.in.fa')
        outfile = os.path.join(data_dir, 'cluster_test_rename_scaffolds.out.fa')
        tmpfile = 'tmp.fa'
        c._rename_scaffolds(infile, tmpfile)
        self.assertTrue(filecmp.cmp(outfile, tmpfile, shallow=False))
        os.unlink(tmpfile)
        clean_cluster_dir(cluster_dir)


    def test_fix_contig_orientation(self):
        '''test _fix_contig_orientation'''
        cluster_dir = os.path.join(data_dir, 'cluster_test_fix_contig_orientation')
        clean_cluster_dir(cluster_dir)
        c = cluster.Cluster(cluster_dir, 'name')
        scaffs_in = os.path.join(data_dir, 'cluster_test_fix_contig_orientation.in.fa')
        scaffs_out = os.path.join(data_dir, 'cluster_test_fix_contig_orientation.out.fa')
        shutil.copyfile(scaffs_in, c.gapfilled_scaffolds)
        shutil.copyfile(os.path.join(data_dir, 'cluster_test_fix_contig_orientation.gene.fa'), c.gene_fa)
        c._fix_contig_orientation()
        self.assertTrue(filecmp.cmp(scaffs_out, c.final_assembly_fa, shallow=False))
        clean_cluster_dir(cluster_dir)


    def test_load_final_contigs(self):
        '''test _load_final_contigs'''
        cluster_dir = os.path.join(data_dir, 'cluster_test_load_final_contigs')
        clean_cluster_dir(cluster_dir)
        c = cluster.Cluster(cluster_dir, 'name')
        contigs_file = os.path.join(data_dir, 'cluster_test_load_final_contigs.contigs.fa')
        shutil.copyfile(contigs_file, c.final_assembly_fa)
        c._load_final_contigs()
        expected = {
            'spam': pyfastaq.sequences.Fasta('spam', 'ACGT'),
            'egg1': pyfastaq.sequences.Fasta('egg1', 'TGCA'),
            'egg2': pyfastaq.sequences.Fasta('egg2', 'AAAA'),
        }
        self.assertEqual(expected, c.final_assembly)
        clean_cluster_dir(cluster_dir)


    def test_parse_assembly_vs_gene_coords(self):
        '''test _parse_assembly_vs_gene_coords'''
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


    def test_parse_assembly_bam(self):
        '''test _parse_assembly_bam'''
        cluster_dir = os.path.join(data_dir, 'cluster_test_parse_assembly_bam')
        clean_cluster_dir(cluster_dir)
        c = cluster.Cluster(cluster_dir, 'name')
        bam = os.path.join(data_dir, 'cluster_test_parse_assembly_bam.bam')
        assembly_fa = os.path.join(data_dir, 'cluster_test_parse_assembly_bam.assembly.fa')
        shutil.copyfile(bam, c.final_assembly_bam)
        shutil.copy(assembly_fa, c.final_assembly_fa)
        c._load_final_contigs()
        c._parse_assembly_bam()
        for e in ['scaff', 'soft_clipped', 'unmapped_mates']:
            self.assertTrue(os.path.exists(c.final_assembly_bam + '.' + e))
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


    def test_nucmer_hits_to_scaff_coords(self):
        '''test _nucmer_hits_to_scaff_coords'''
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


    def test_nucmer_hits_to_gene_cov_per_contig(self):
        '''test _nucmer_hits_to_gene_cov_per_contig'''
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


    def test_gene_coverage_unique(self):
        '''test _gene_coverage_unique'''
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


    def test_gene_covered_by_complete_contig_with_orf(self):
        '''test _gene_covered_by_complete_contig_with_orf'''
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


    def test_gene_covered_by_at_least_one_full_length_contig(self):
        '''test _gene_covered_by_at_least_one_full_length_contig'''
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


    def test_get_mummer_variants(self):
        '''test _get_mummer_variants'''
        cluster_dir = os.path.join(data_dir, 'cluster_test_generic')
        clean_cluster_dir(cluster_dir)
        c = cluster.Cluster(cluster_dir, 'name')
        snp_file = os.path.join(data_dir, 'cluster_test_get_mummer_variants.none.snps')
        shutil.copyfile(snp_file, c.assembly_vs_gene_coords + '.snps')
        c._get_mummer_variants()
        self.assertEqual(c.variants, {})

        clean_cluster_dir(cluster_dir)
        snp_file = os.path.join(data_dir, 'cluster_test_get_mummer_variants.snp.snps')
        shutil.copyfile(snp_file, c.assembly_vs_gene_coords + '.snps')
        v1 = pymummer.variant.Variant(pymummer.snp.Snp('42\tA\tG\t42\t42\t42\t500\t500\t1\t1\tgene\tcontig1'))
        v2 = pymummer.variant.Variant(pymummer.snp.Snp('42\tA\tG\t42\t42\t42\t500\t500\t1\t1\tgene\tcontig2'))
        v3 = pymummer.variant.Variant(pymummer.snp.Snp('40\tT\tC\t40\t42\t42\t500\t500\t1\t1\tgene\tcontig1'))
        v4 = pymummer.variant.Variant(pymummer.snp.Snp('2\tC\tG\t2\t42\t42\t500\t500\t1\t1\tgene\tcontig1'))
        expected = {
            'contig1': [[v4], [v3, v1]],
            'contig2': [[v2]]
        }
        shutil.copyfile(snp_file, c.assembly_vs_gene_coords + '.snps')
        c._get_mummer_variants()
        self.assertEqual(c.variants, expected)
        clean_cluster_dir(cluster_dir)


    def test_filter_mummer_variants(self):
        '''test filter_mummer_variants'''
        cluster_dir = os.path.join(data_dir, 'cluster_test_generic')
        clean_cluster_dir(cluster_dir)
        c = cluster.Cluster(cluster_dir, 'name')
        c.gene = pyfastaq.sequences.Fasta('gene', 'GATCGCGAAGCGATGACCCATGAAGCGACCGAACGCTGA')
        v1 = pymummer.variant.Variant(pymummer.snp.Snp('4\tC\tT\t4\tx\tx\t39\t39\tx\tx\tgene\tcontig'))
        v2 = pymummer.variant.Variant(pymummer.snp.Snp('6\tC\tA\t6\tx\tx\t39\t39\tx\tx\tgene\tcontig'))
        v3 = pymummer.variant.Variant(pymummer.snp.Snp('12\tG\tT\t12\tx\tx\t39\t39\tx\tx\tgene\tcontig'))
        c.variants = {'contig': [[v1, v2], v3]}
        c._filter_mummer_variants()
        expected = {'contig': [[v1, v2]]}
        self.assertEqual(expected, c.variants)
        clean_cluster_dir(cluster_dir)


    def test_get_codon_start(self):
        '''test _get_codon_start'''
        cluster_dir = os.path.join(data_dir, 'cluster_test_generic')
        clean_cluster_dir(cluster_dir)
        c = cluster.Cluster(cluster_dir, 'name')
        tests = [
            (0, 5, 3),
            (0, 0, 0),
            (0, 1, 0),
            (0, 2, 0),
            (1, 3, 1),
            (2, 3, 2),
            (3, 3, 3),
            (3, 6, 6),
            (3, 7, 6),
            (3, 8, 6),
        ]
        for t in tests:
            self.assertEqual(c._get_codon_start(t[0], t[1]), t[2])
        clean_cluster_dir(cluster_dir)


    def test_get_variant_effect(self):
        '''test _get_variant_effect'''
        cluster_dir = os.path.join(data_dir, 'cluster_test_generic')
        clean_cluster_dir(cluster_dir)
        c = cluster.Cluster(cluster_dir, 'name')
        c.gene = pyfastaq.sequences.Fasta('gene', 'GATCGCGAAGCGATGACCCATGAAGCGACCGAACGCTGA')
        v1 = pymummer.variant.Variant(pymummer.snp.Snp('6\tC\tT\t6\tx\tx\t39\t39\tx\tx\tgene\tcontig'))
        v1 = pymummer.variant.Variant(pymummer.snp.Snp('6\tC\tT\t6\tx\tx\t39\t39\tx\tx\tgene\tcontig'))
        v2 = pymummer.variant.Variant(pymummer.snp.Snp('4\tC\tA\t6\tx\tx\t39\t39\tx\tx\tgene\tcontig'))
        v3 = pymummer.variant.Variant(pymummer.snp.Snp('4\tC\tT\t4\tx\tx\t39\t39\tx\tx\tgene\tcontig'))
        v4 = pymummer.variant.Variant(pymummer.snp.Snp('6\tC\tA\t6\tx\tx\t39\t39\tx\tx\tgene\tcontig'))
        v5 = pymummer.variant.Variant(pymummer.snp.Snp('4\tC\t.\t4\tx\tx\t39\t39\tx\tx\tgene\tcontig'))
        v6 = pymummer.variant.Variant(pymummer.snp.Snp('4\t.\tA\t4\tx\tx\t39\t39\tx\tx\tgene\tcontig'))
        v7 = pymummer.variant.Variant(pymummer.snp.Snp('4\t.\tG\t4\tx\tx\t39\t39\tx\tx\tgene\tcontig'))
        v7.qry_base = 'GAT'
        v8 = pymummer.variant.Variant(pymummer.snp.Snp('4\t.\tG\t4\tx\tx\t39\t39\tx\tx\tgene\tcontig'))
        v8.qry_base = 'TGA'
        v9 = pymummer.variant.Variant(pymummer.snp.Snp('4\t.\tG\t4\tx\tx\t39\t39\tx\tx\tgene\tcontig'))
        v9.qry_base = 'ATTCCT'
        v10 = pymummer.variant.Variant(pymummer.snp.Snp('4\tC\t.\t4\tx\tx\t39\t39\tx\tx\tgene\tcontig'))
        v10.ref_base = 'CGC'
        v10.ref_end = 5
        v11 = pymummer.variant.Variant(pymummer.snp.Snp('4\tC\t.\t4\tx\tx\t39\t39\tx\tx\tgene\tcontig'))
        v11.ref_base = 'CGCGAA'
        v11.ref_end = 8

        variants = [
            ([v1], ('SYN', '.')),
            ([v2], ('NONSYN', 'R2S')),
            ([v2, v1], ('NONSYN', 'R2S')),
            ([v3, v4], ('TRUNC', 'R2trunc')),
            ([v5], ('FSHIFT', 'R2fs')),
            ([v6], ('FSHIFT', 'R2fs')),
            ([v7], ('INS', 'R2_E3insD')),
            ([v8], ('TRUNC', 'R2trunc')),
            ([v9], ('INS', 'R2_E3insIP')),
            ([v10], ('DEL', 'R2del')),
            ([v11], ('DEL', 'R2_E3del')),
        ]

        for t in variants:
            self.assertEqual(t[1], c._get_variant_effect(t[0]))

        clean_cluster_dir(cluster_dir)


    def test_make_assembly_vcf(self):
        '''test _make_assembly_vcf'''
        cluster_dir = os.path.join(data_dir, 'cluster_test_generic')
        clean_cluster_dir(cluster_dir)
        c = cluster.Cluster(cluster_dir, 'name')
        c.final_assembly_fa = os.path.join(data_dir, 'cluster_test_make_assembly_vcf.assembly.fa')
        c.final_assembly_bam = os.path.join(data_dir, 'cluster_test_make_assembly_vcf.assembly.bam')
        expected_vcf = os.path.join(data_dir, 'cluster_test_make_assembly_vcf.assembly.vcf')
        expected_depths = os.path.join(data_dir, 'cluster_test_make_assembly_vcf.assembly.read_depths')
        c._make_assembly_vcf()

        def get_vcf_call_lines(fname):
            with open(fname) as f:
                lines = [x for x in f.readlines() if not x.startswith('#')]
            return lines

        expected_lines = get_vcf_call_lines(expected_vcf)
        got_lines = get_vcf_call_lines(c.final_assembly_vcf)
        self.assertEqual(expected_lines, got_lines)
        self.assertTrue(filecmp.cmp(expected_depths, c.final_assembly_read_depths, shallow=False))
        clean_cluster_dir(cluster_dir)


    def test_get_vcf_variant_counts(self):
        '''test _get_vcf_variant_counts'''
        cluster_dir = os.path.join(data_dir, 'cluster_test_generic')
        clean_cluster_dir(cluster_dir)
        c = cluster.Cluster(cluster_dir, 'name')
        hit = ['1', '42', '1', '42', '42', '42', '100.00', '1000', '1000', '1', '1', 'gene', 'scaff1']
        c.nucmer_hits = {
            'scaff1': [pymummer.alignment.Alignment('\t'.join(hit))]
        }

        c.final_assembly_vcf = os.path.join(data_dir, 'cluster_test_get_vcf_variant_counts.vcf')
        c._get_vcf_variant_counts()
        expected = {'scaff1': 1}
        self.assertEqual(expected, c.vcf_variant_counts)
        clean_cluster_dir(cluster_dir)


    def test_make_report_lines(self):
        '''test _make_report_lines'''
        cluster_dir = os.path.join(data_dir, 'cluster_test_generic')
        clean_cluster_dir(cluster_dir)
        c = cluster.Cluster(cluster_dir, 'cluster_name')
        c.gene = pyfastaq.sequences.Fasta('gene', 'GATCGCGAAGCGATGACCCATGAAGCGACCGAACGCTGA')
        v1 = pymummer.variant.Variant(pymummer.snp.Snp('6\tC\tT\t6\tx\tx\t39\t39\tx\tx\tgene\tcontig'))

        nucmer_hit = ['1', '10', '1', '10', '10', '10', '90.00', '1000', '1000', '1', '1', 'gene', 'contig']
        c.nucmer_hits = {'contig': [pymummer.alignment.Alignment('\t'.join(nucmer_hit))]}
        c.variants = {'contig': [[v1]]}
        c.percent_identities = {'contig': 92.42}
        c.status_flag.set_flag(42)
        c._make_report_lines()
        expected = [[
            'gene',
            42,
            2,
            'cluster_name',
            39,
            10,
            92.42,
            'SNP',
            'SYN',
            '.',
            6,
            6,
            'C',
            'contig',
            39,
            6,
            6,
            'T',
        ]]
        self.assertEqual(expected, c.report_lines)
        clean_cluster_dir(cluster_dir)
