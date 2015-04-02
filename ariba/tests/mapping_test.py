import unittest
import os
import shutil
import pysam
from ariba import mapping

modules_dir = os.path.dirname(os.path.abspath(mapping.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


# different smalt version output slightly different BAMs. Some columns
# should never change, so check just those ones
def get_sam_columns(bamfile):
    sams = []
    sam_reader = pysam.Samfile(bamfile, "rb")
    for sam in sam_reader.fetch(until_eof=True):
        if sam.is_unmapped:
            refname = None
        else:
            refname = sam_reader.getrname(sam.tid)
        sams.append((sam.qname, sam.flag, refname, sam.pos, sam.cigar, sam.seq))
    return sams


class TestMapping(unittest.TestCase):
    def test_bowtie2_in_path(self):
        '''Test that bowtie2 is in the user's path'''
        assert(shutil.which('bowtie2') is not None)


    def test_samtools_in_path(self):
        '''Test that samtools is in the user's path'''
        assert(shutil.which('samtools') is not None)


    def test_run_bowtie2(self):
        '''Test run_bowtie2 unsorted'''
        self.maxDiff = None
        ref = os.path.join(data_dir, 'mapping_test_bowtie2_ref.fa')
        reads1 = os.path.join(data_dir, 'mapping_test_bowtie2_reads_1.fq')
        reads2 = os.path.join(data_dir, 'mapping_test_bowtie2_reads_2.fq')
        out_prefix = 'tmp.out.bowtie2'
        mapping.run_bowtie2(reads1, reads2, ref, out_prefix)
        expected = get_sam_columns(os.path.join(data_dir, 'mapping_test_bowtie2_unsorted.bam'))
        got = get_sam_columns(out_prefix + '.bam')
        self.assertListEqual(expected, got)
        os.unlink(out_prefix + '.bam')


    def test_run_bowtie2_and_sort(self):
        '''Test run_bowtie2 sorted'''
        ref = os.path.join(data_dir, 'mapping_test_bowtie2_ref.fa')
        reads1 = os.path.join(data_dir, 'mapping_test_bowtie2_reads_1.fq')
        reads2 = os.path.join(data_dir, 'mapping_test_bowtie2_reads_2.fq')
        out_prefix = 'tmp.out.bowtie2'
        mapping.run_bowtie2(reads1, reads2, ref, out_prefix, sort=True)
        expected = get_sam_columns(os.path.join(data_dir, 'mapping_test_bowtie2_sorted.bam'))
        got = get_sam_columns(out_prefix + '.bam')
        self.assertListEqual(expected, got)
        os.unlink(out_prefix + '.bam')
        os.unlink(out_prefix + '.bam.bai')
        os.unlink(out_prefix + '.unsorted.bam')


    #def test_run_smalt(self):
    #    '''Test run_smalt unsorted'''
    #    ref = os.path.join(data_dir, 'mapping_test_smalt_ref.fa')
    #    reads1 = os.path.join(data_dir, 'mapping_test_smalt_reads_1.fq')
    #    reads2 = os.path.join(data_dir, 'mapping_test_smalt_reads_2.fq')
    #    out_prefix = 'tmp.out.smalt'
    #    mapping.run_smalt(reads1, reads2, ref, out_prefix)
    #    expected = get_sam_columns(os.path.join(data_dir, 'mapping_test_smalt_unsorted.bam'))
    #    got = get_sam_columns(out_prefix + '.bam')
    #    self.assertListEqual(expected, got)
    #    os.unlink(out_prefix + '.bam')


    #def test_run_smalt_and_sort(self):
    #    '''Test run_smalt sorted'''
    #    ref = os.path.join(data_dir, 'mapping_test_smalt_ref.fa')
    #    reads1 = os.path.join(data_dir, 'mapping_test_smalt_reads_1.fq')
    #    reads2 = os.path.join(data_dir, 'mapping_test_smalt_reads_2.fq')
    #    out_prefix = 'tmp.out.smalt'
    #    mapping.run_smalt(reads1, reads2, ref, out_prefix, sort=True)
    #    expected = get_sam_columns(os.path.join(data_dir, 'mapping_test_smalt_sorted.bam'))
    #    got = get_sam_columns(out_prefix + '.bam')
    #    self.assertListEqual(expected, got)
    #    os.unlink(out_prefix + '.bam')
    #    os.unlink(out_prefix + '.bam.bai')
    #    os.unlink(out_prefix + '.unsorted.bam')


    def test_get_total_alignment_score(self):
        '''Test get_total_alignment_score'''
        bam = os.path.join(data_dir, 'mapping_test_get_total_alignment_score.bam')
        expected = 219
        got = mapping.get_total_alignment_score(bam)
        self.assertEqual(got, expected)

