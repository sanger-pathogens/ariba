import unittest
import os
import filecmp
import pyfastaq
from ariba import assembly, common
from ariba import external_progs

modules_dir = os.path.dirname(os.path.abspath(assembly.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')
extern_progs = external_progs.ExternalProgs(using_spades=True)

class TestAssembly(unittest.TestCase):
    def test_run_fermilite(self):
        '''test _run_fermilite'''
        reads = os.path.join(data_dir, 'assembly_run_fermilite.reads.fq')
        tmp_fa = 'tmp.test_run_fermilite.fa'
        tmp_log = 'tmp.test_run_fermilite.log'
        expected_fa = os.path.join(data_dir, 'assembly_run_fermilite.expected.fa')
        expected_log = os.path.join(data_dir, 'assembly_run_fermilite.expected.log')
        got = assembly.Assembly._run_fermilite(reads, tmp_fa, tmp_log, 'contig')
        self.assertEqual(0, got)
        self.assertTrue(filecmp.cmp(expected_fa, tmp_fa, shallow=False))
        self.assertTrue(filecmp.cmp(expected_log, tmp_log, shallow=False))
        os.unlink(tmp_fa)
        os.unlink(tmp_log)


    def test_run_fermilite_fails(self):
        '''test _run_fermilite when it fails'''
        reads = os.path.join(data_dir, 'assembly_run_fermilite_fail.reads.fq')
        tmp_fa = 'tmp.test_run_fermilite_fails.fa'
        tmp_log = 'tmp.test_run_fermilite_fails.log'
        expected_log = os.path.join(data_dir, 'assembly_run_fermilite_fails.expected.log')
        got = assembly.Assembly._run_fermilite(reads, tmp_fa, tmp_log, 'contig')
        self.assertEqual(1, got)
        self.assertFalse(os.path.exists(tmp_fa))
        self.assertTrue(filecmp.cmp(expected_log, tmp_log, shallow=False))
        os.unlink(tmp_log)


    def test_assemble_with_fermilite(self):
        '''test _assemble_with_fermilite'''
        reads1 = os.path.join(data_dir, 'assembly_assemble_with_fermilite.reads_1.fq')
        reads2 = os.path.join(data_dir, 'assembly_assemble_with_fermilite.reads_2.fq')
        expected_log = os.path.join(data_dir, 'assembly_assemble_with_fermilite.expected.log')
        expected_fa = os.path.join(data_dir, 'assembly_assemble_with_fermilite.expected.fa')
        tmp_dir = 'tmp.test_assemble_with_fermilite'
        tmp_log = 'tmp.test_assemble_with_fermilite.log'
        tmp_log_fh = open(tmp_log, 'w')
        print('First line', file=tmp_log_fh)
        a = assembly.Assembly(reads1, reads2, 'not needed', 'not needed', tmp_dir, 'not_needed_for_this_test.fa', 'not_needed_for_this_test.bam', tmp_log_fh, 'not needed')
        a._assemble_with_fermilite()
        self.assertTrue(a.assembled_ok)
        tmp_log_fh.close()
        self.assertTrue(filecmp.cmp(expected_log, tmp_log, shallow=False))
        self.assertTrue(filecmp.cmp(expected_fa, os.path.join(tmp_dir, 'debug_all_contigs.fa'), shallow=False))
        common.rmtree(tmp_dir)
        os.unlink(tmp_log)


    def test_assemble_with_fermilite_fails(self):
        '''test _assemble_with_fermilite fails'''
        reads1 = os.path.join(data_dir, 'assembly_assemble_with_fermilite_fails.reads_1.fq')
        reads2 = os.path.join(data_dir, 'assembly_assemble_with_fermilite_fails.reads_2.fq')
        expected_log = os.path.join(data_dir, 'assembly_assemble_with_fermilite_fails.expected.log')
        tmp_dir = 'tmp.test_assemble_with_fermilite_fails'
        tmp_log = 'tmp.test_assemble_with_fermilite_fails.log'
        tmp_log_fh = open(tmp_log, 'w')
        print('First line', file=tmp_log_fh)
        a = assembly.Assembly(reads1, reads2, 'not needed', 'not needed', tmp_dir, 'not_needed_for_this_test.fa', 'not_needed_for_this_test.bam', tmp_log_fh, 'not needed')
        a._assemble_with_fermilite()
        self.assertFalse(a.assembled_ok)
        tmp_log_fh.close()
        self.assertTrue(filecmp.cmp(expected_log, tmp_log, shallow=False))
        self.assertFalse(os.path.exists(os.path.join(tmp_dir, 'contigs.fa')))
        common.rmtree(tmp_dir)
        os.unlink(tmp_log)


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

    def test_check_spades_log_file(self):
        '''test _check_spades_log_file'''
        good_file = os.path.join(data_dir, 'assembly_test_check_spades_log_file.log.good')
        bad_file = os.path.join(data_dir, 'assembly_test_check_spades_log_file.log.bad')
        self.assertTrue(assembly.Assembly._check_spades_log_file(good_file))
        with self.assertRaises(assembly.Error):
            self.assertTrue(assembly.Assembly._check_spades_log_file(bad_file))

    @unittest.skipUnless(extern_progs.exe('spades'), "Spades assembler is optional and is not configured")
    def test_assemble_with_spades(self):
        '''test _assemble_with_spades'''
        reads1 = os.path.join(data_dir, 'assembly_test_assemble_with_spades_reads_1.fq')
        reads2 = os.path.join(data_dir, 'assembly_test_assemble_with_spades_reads_2.fq')
        tmp_dir = 'tmp.test_assemble_with_spades'
        tmp_log = 'tmp.test_assemble_with_spades.log'
        with open(tmp_log, 'w') as tmp_log_fh:
            print('First line', file=tmp_log_fh)
            common.rmtree(tmp_dir)
            #using spades_options=" --only-assembler" because error correction cannot determine quality offset on this
            #artificial dataset
            a = assembly.Assembly(reads1, reads2, 'not needed', 'not needed', tmp_dir, 'not_needed_for_this_test.fa',
                                  'not_needed_for_this_test.bam', tmp_log_fh, 'not needed',
                                  assembler="spades", spades_options=" --only-assembler")
            a._assemble_with_spades()
        self.assertTrue(a.assembled_ok)
        common.rmtree(tmp_dir)
        os.unlink(tmp_log)

    @unittest.skipUnless(extern_progs.exe('spades'), "Spades assembler is optional and is not configured")
    def test_assemble_with_spades_fail(self):
        '''test _assemble_with_spades handles spades fail'''
        reads1 = os.path.join(data_dir, 'assembly_test_assemble_with_spades_fails_reads_1.fq')
        reads2 = os.path.join(data_dir, 'assembly_test_assemble_with_spades_fails_reads_2.fq')
        tmp_dir = 'tmp.test_assemble_with_spades_fail'
        tmp_log = 'tmp.test_assemble_with_spades_fail.log'
        with open(tmp_log, 'w') as tmp_log_fh:
            print('First line', file=tmp_log_fh)
            common.rmtree(tmp_dir)
            a = assembly.Assembly(reads1, reads2, 'not needed', 'not needed', tmp_dir, 'not_needed_for_this_test.fa',
                                  'not_needed_for_this_test.bam', tmp_log_fh, 'not needed',
                                  assembler="spades", spades_options=" --only-assembler")
            a._assemble_with_spades()
        self.assertFalse(a.assembled_ok)
        common.rmtree(tmp_dir)
        os.unlink(tmp_log)
