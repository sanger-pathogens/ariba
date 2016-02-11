import unittest
import os
import shutil
import pysam
from ariba import mapping, external_progs

modules_dir = os.path.dirname(os.path.abspath(mapping.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')
extern_progs = external_progs.ExternalProgs()


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
    def test_run_bowtie2(self):
        '''Test run_bowtie2 unsorted'''
        self.maxDiff = None
        ref = os.path.join(data_dir, 'mapping_test_bowtie2_ref.fa')
        reads1 = os.path.join(data_dir, 'mapping_test_bowtie2_reads_1.fq')
        reads2 = os.path.join(data_dir, 'mapping_test_bowtie2_reads_2.fq')
        out_prefix = 'tmp.out.bowtie2'
        mapping.run_bowtie2(
            reads1,
            reads2,
            ref,
            out_prefix,
            samtools=extern_progs.exe('samtools'),
            bowtie2=extern_progs.exe('bowtie2'),
        )
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
        mapping.run_bowtie2(
            reads1,
            reads2,
            ref,
            out_prefix,
            sort=True,
            samtools=extern_progs.exe('samtools'),
            bowtie2=extern_progs.exe('bowtie2'),
        )
        expected = get_sam_columns(os.path.join(data_dir, 'mapping_test_bowtie2_sorted.bam'))
        got = get_sam_columns(out_prefix + '.bam')
        self.assertListEqual(expected, got)
        os.unlink(out_prefix + '.bam')
        os.unlink(out_prefix + '.bam.bai')


    def test_get_total_alignment_score(self):
        '''Test get_total_alignment_score'''
        bam = os.path.join(data_dir, 'mapping_test_get_total_alignment_score.bam')
        expected = 219
        got = mapping.get_total_alignment_score(bam)
        self.assertEqual(got, expected)

