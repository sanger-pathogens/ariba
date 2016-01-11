import unittest
import os
import copy
import shutil
import filecmp
import pyfastaq
import pysam
import pymummer
from ariba import assembly

modules_dir = os.path.dirname(os.path.abspath(assembly.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


class TestAssemblye(unittest.TestCase):
    def test_set_assembly_kmer(self):
        '''test _set_assembly_kmer'''
        reads1 = os.path.join(data_dir, 'assembly_test_set_assembly_kmer_reads_1.fq')
        reads2 = os.path.join(data_dir, 'assembly_test_set_assembly_kmer_reads_2.fq')
        got = assembly.Assembly._set_assembly_kmer(0, reads1, reads2)
        self.assertEqual(got, 5)
        got = assembly.Assembly._set_assembly_kmer(42, reads1, reads2)
        self.assertEqual(got, 42)


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


    def test_has_gaps_to_fill(self):
        '''test _has_gaps_to_fill'''
        # FIXME
        pass


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


    def test_parse_bam(self):
        '''test _parse_bam'''
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


