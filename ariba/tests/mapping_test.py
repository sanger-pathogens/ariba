import unittest
import os
import pysam
import pyfastaq
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
    def test_bowtie2_index(self):
        '''test bowtie2_index'''
        tmp_ref = 'tmp.test_bowtie2_index.ref.fa'
        with open(tmp_ref, 'w') as f:
            print('>ref', file=f)
            print('ATCATACTACTCATACTGACTCATCATCATCATGACGTATG', file=f)

        tmp_out = 'tmp.test_bowtie2_index.ref.out'
        mapping.bowtie2_index(tmp_ref, tmp_out, bowtie2=extern_progs.exe('bowtie2'))

        expected_files = [tmp_out + '.' + x + '.bt2' for x in ['1', '2', '3', '4', 'rev.1', 'rev.2']]
        for filename in expected_files:
            self.assertTrue(os.path.exists(filename))
            os.unlink(filename)

        os.unlink(tmp_ref)


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
            bowtie2=extern_progs.exe('bowtie2'),
            bowtie2_version=extern_progs.version('bowtie2'),
        )
        expected = get_sam_columns(os.path.join(data_dir, 'mapping_test_bowtie2_unsorted.bam'))
        got = get_sam_columns(out_prefix + '.bam')
        self.assertListEqual(expected, got)
        os.unlink(out_prefix + '.bam')


    def test_run_bowtie2_remove_both_unmapped(self):
        '''Test run_bowtie2 unsorted remove both unmapped'''
        self.maxDiff = None
        ref = os.path.join(data_dir, 'mapping_test_bowtie2_ref.fa')
        reads1 = os.path.join(data_dir, 'mapping_test_bowtie2_remove_both_unmapped_reads_1.fq')
        reads2 = os.path.join(data_dir, 'mapping_test_bowtie2_remove_both_unmapped_reads_2.fq')
        out_prefix = 'tmp.out.bowtie2_remove_both_unmapped'
        mapping.run_bowtie2(
            reads1,
            reads2,
            ref,
            out_prefix,
            bowtie2=extern_progs.exe('bowtie2'),
            bowtie2_version=extern_progs.version('bowtie2'),
            remove_both_unmapped=True,
        )
        expected = get_sam_columns(os.path.join(data_dir, 'mapping_test_bowtie2_remove_both_unmapped_reads.bam'))
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
            bowtie2=extern_progs.exe('bowtie2'),
            bowtie2_version=extern_progs.version('bowtie2'),
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


    def test_sam_to_fastq(self):
        '''test sam_to_fastq'''
        expected = [
            pyfastaq.sequences.Fastq('read1/1', 'GTATGAGTAGATATAAAGTCCGGAACTGTGATCGGGGGCGATTTATTTACTGGCCGTCCC', 'GHIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII'),
            pyfastaq.sequences.Fastq('read1/2', 'TCCCATACGTTGCAATCTGCAGACGCCACTCTTCCACGTCGGACGAACGCAACGTCAGGA', 'IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIHGEDCBA')
        ]


        sam_reader = pysam.Samfile(os.path.join(data_dir, 'mapping_test_sam_to_fastq.bam'), "rb")
        i = 0
        for s in sam_reader.fetch(until_eof=True):
            self.assertEqual(expected[i], mapping.sam_to_fastq(s))
            i += 1


    def test_sam_pair_to_insert(self):
        '''test sam_pair_to_insert'''
        expected = [
            None, # both unmapped
            None, # read 1 unmapped
            None, # read 2 unmpapped
            None, # mapped to different seqs
            None, # same seqs, wrond orientation
            660
        ]

        sam1 = None
        i = 0
        sam_reader = pysam.Samfile(os.path.join(data_dir, 'mapping_test_sam_pair_to_insert.bam'), 'rb')
        for s in sam_reader.fetch(until_eof=True):
            if sam1 is None:
                sam1 = s
                continue

            self.assertEqual(mapping.sam_pair_to_insert(s, sam1), expected[i])
            sam1 = None
            i += 1


