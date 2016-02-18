import unittest
import sys
import os
import shutil
import filecmp
import pyfastaq
from ariba import assembly

modules_dir = os.path.dirname(os.path.abspath(assembly.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


class TestAssembly(unittest.TestCase):
    def test_get_assembly_kmer(self):
        '''test _get_assembly_kmer'''
        reads1 = os.path.join(data_dir, 'assembly_test_set_assembly_kmer_reads_1.fq')
        reads2 = os.path.join(data_dir, 'assembly_test_set_assembly_kmer_reads_2.fq')
        got = assembly.Assembly._get_assembly_kmer(0, reads1, reads2)
        self.assertEqual(got, 5)
        got = assembly.Assembly._get_assembly_kmer(42, reads1, reads2)
        self.assertEqual(got, 42)


    def test_assemble_with_spades(self):
        '''test _assemble_with_spades'''
        reads1 = os.path.join(data_dir, 'assembly_test_assemble_with_spades_reads_1.fq')
        reads2 = os.path.join(data_dir, 'assembly_test_assemble_with_spades_reads_2.fq')
        ref_fasta = os.path.join(data_dir, 'assembly_test_assemble_with_spades_ref.fa')
        tmp_dir = 'tmp.test_assemble_with_spades'
        a = assembly.Assembly(reads1, reads2, ref_fasta, tmp_dir, 'not_needed_for_this_test.fa', 'not_needed_for_this_test.bam', sys.stdout)
        a._assemble_with_spades(unittest=True)
        self.assertTrue(a.assembled_ok)
        shutil.rmtree(tmp_dir)


    def test_assemble_with_spades_fail(self):
        '''test _assemble_with_spades handles spades fail'''
        reads1 = os.path.join(data_dir, 'assembly_test_assemble_with_spades_reads_1.fq')
        reads2 = os.path.join(data_dir, 'assembly_test_assemble_with_spades_reads_2.fq')
        ref_fasta = os.path.join(data_dir, 'assembly_test_assemble_with_spades_ref.fa')
        tmp_dir = 'tmp.test_assemble_with_spades'
        a = assembly.Assembly(reads1, reads2, ref_fasta, tmp_dir, 'not_needed_for_this_test.fa', 'not_needed_for_this_test.bam', sys.stdout)
        a._assemble_with_spades(unittest=False)
        self.assertFalse(a.assembled_ok)
        shutil.rmtree(tmp_dir)


    def test_scaffold_with_sspace(self):
        '''test _scaffold_with_sspace'''
        reads1 = os.path.join(data_dir, 'assembly_test_assemble_with_spades_reads_1.fq')
        reads2 = os.path.join(data_dir, 'assembly_test_assemble_with_spades_reads_2.fq')
        ref_fasta = os.path.join(data_dir, 'assembly_test_assemble_with_spades_ref.fa')
        tmp_dir = 'tmp.test_scaffold_with_sspace'
        a = assembly.Assembly(reads1, reads2, ref_fasta, tmp_dir, 'not_needed_for_this_test.fa', 'not_needed_for_this_test.bam', sys.stdout)
        a.assembly_contigs = os.path.join(data_dir, 'assembly_test_scaffold_with_sspace_contigs.fa')
        a._scaffold_with_sspace()
        self.assertTrue(os.path.exists(a.scaffolder_scaffolds))
        shutil.rmtree(tmp_dir)


    def test_has_gaps_to_fill(self):
        '''test _has_gaps_to_fill'''
        no_gaps = os.path.join(data_dir, 'assembly_test_has_gaps_to_fill.no_gaps.fa')
        has_gaps = os.path.join(data_dir, 'assembly_test_has_gaps_to_fill.has_gaps.fa')
        self.assertTrue(assembly.Assembly._has_gaps_to_fill(has_gaps))
        self.assertFalse(assembly.Assembly._has_gaps_to_fill(no_gaps))


    def test_rename_scaffolds(self):
        '''test _rename_scaffolds'''
        infile = os.path.join(data_dir, 'assembly_test_rename_scaffolds.in.fa')
        outfile = os.path.join(data_dir, 'assembly_test_rename_scaffolds.out.fa')
        tmpfile = 'tmp.fa'
        assembly.Assembly._rename_scaffolds(infile, tmpfile, 'prefix')
        self.assertTrue(filecmp.cmp(outfile, tmpfile, shallow=False))
        os.unlink(tmpfile)


    def test_gap_fill_with_gapfiller_no_gaps(self):
        '''test _gap_fill_with_gapfiller no gaps'''
        reads1 = os.path.join(data_dir, 'assembly_test_gapfill_with_gapfiller_reads_1.fq')
        reads2 = os.path.join(data_dir, 'assembly_test_gapfill_with_gapfiller_reads_2.fq')
        tmp_dir = 'tmp.gap_fill_with_gapfiller_no_gaps'
        a = assembly.Assembly(reads1, reads2, 'ref.fa', tmp_dir, 'not_needed_for_this_test.fa', 'not_needed_for_this_test.bam', sys.stdout)
        a.scaffolder_scaffolds = os.path.join(data_dir, 'assembly_test_gapfill_with_gapfiller.scaffolds_no_gaps.fa')
        a._gap_fill_with_gapfiller()
        self.assertTrue(os.path.exists(a.gapfilled_scaffolds))
        shutil.rmtree(tmp_dir)


    def test_gap_fill_with_gapfiller_with_gaps(self):
        '''test _gap_fill_with_gapfiller with gaps'''
        reads1 = os.path.join(data_dir, 'assembly_test_gapfill_with_gapfiller_reads_1.fq')
        reads2 = os.path.join(data_dir, 'assembly_test_gapfill_with_gapfiller_reads_2.fq')
        tmp_dir = 'tmp.gap_fill_with_gapfiller_with_gaps'
        a = assembly.Assembly(reads1, reads2, 'ref.fa', tmp_dir, 'not_needed_for_this_test.fa', 'not_needed_for_this_test.bam', sys.stdout)
        a.scaffolder_scaffolds = os.path.join(data_dir, 'assembly_test_gapfill_with_gapfiller.scaffolds_with_gaps.fa')
        a._gap_fill_with_gapfiller()
        self.assertTrue(os.path.exists(a.gapfilled_scaffolds))
        shutil.rmtree(tmp_dir)


    def test_fix_contig_orientation(self):
        '''test _fix_contig_orientation'''
        scaffs_in = os.path.join(data_dir, 'assembly_test_fix_contig_orientation.in.fa')
        expected_out = os.path.join(data_dir, 'assembly_test_fix_contig_orientation.out.fa')
        ref_fa = os.path.join(data_dir, 'assembly_test_fix_contig_orientation.ref.fa')
        tmp_out = 'tmp.assembly_test_fix_contig_orientation.out.fa'
        got = assembly.Assembly._fix_contig_orientation(scaffs_in, ref_fa, tmp_out)
        expected = {'match_both_strands'}
        self.assertTrue(filecmp.cmp(expected_out, tmp_out, shallow=False))
        self.assertEqual(expected, got)
        os.unlink(tmp_out)


    def test_parse_bam(self):
        '''test _parse_bam'''
        bam = os.path.join(data_dir, 'assembly_test_parse_assembly_bam.bam')
        assembly_fa = os.path.join(data_dir, 'assembly_test_parse_assembly_bam.assembly.fa')
        assembly_seqs = {}
        pyfastaq.tasks.file_to_dict(assembly_fa, assembly_seqs)
        self.assertTrue(assembly.Assembly._parse_bam(assembly_seqs, bam, 10, 1000))
        os.unlink(bam + '.soft_clipped')
        os.unlink(bam + '.unmapped_mates')
        os.unlink(bam + '.scaff')

