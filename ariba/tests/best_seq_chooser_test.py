import unittest
import sys
import os
import pyfastaq
from ariba import best_seq_chooser, external_progs

modules_dir = os.path.dirname(os.path.abspath(best_seq_chooser.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')
extern_progs = external_progs.ExternalProgs()


class TestBestSeqChooser(unittest.TestCase):
    def test_total_alignment_score(self):
        '''test _total_alignment_score'''
        reads1 = os.path.join(data_dir, 'best_seq_chooser_total_alignment_score_reads_1.fq')
        reads2 = os.path.join(data_dir, 'best_seq_chooser_total_alignment_score_reads_2.fq')
        ref = os.path.join(data_dir, 'best_seq_chooser_total_alignment_score_ref_seqs.fa')
        chooser = best_seq_chooser.BestSeqChooser(
            reads1,
            reads2,
            ref,
            sys.stdout,
            samtools_exe=extern_progs.exe('samtools'),
            bowtie2_exe=extern_progs.exe('bowtie2'),
        )
        self.assertEqual(3000, chooser._total_alignment_score('1'))


    def test_get_best_seq_by_alignment_score(self):
        '''test _get_best_seq_by_alignment_score'''
        reads1 = os.path.join(data_dir, 'best_seq_chooser_get_best_seq_by_alignment_score_reads_1.fq')
        reads2 = os.path.join(data_dir, 'best_seq_chooser_get_best_seq_by_alignment_score_reads_2.fq')
        ref = os.path.join(data_dir, 'best_seq_chooser_get_best_seq_by_alignment_score_ref.fa')
        chooser = best_seq_chooser.BestSeqChooser(
            reads1,
            reads2,
            ref,
            sys.stdout,
            samtools_exe=extern_progs.exe('samtools'),
            bowtie2_exe=extern_progs.exe('bowtie2'),
        )
        self.assertEqual('1', chooser._get_best_seq_by_alignment_score())


    def test_best_seq(self):
        '''test best_seq'''
        reads1 = os.path.join(data_dir, 'best_seq_chooser_best_seq_reads_1.fq')
        reads2 = os.path.join(data_dir, 'best_seq_chooser_best_seq_reads_2.fq')
        ref = os.path.join(data_dir, 'best_seq_chooser_best_seq_ref.fa')
        expected_seq = pyfastaq.sequences.Fasta('1', ''.join([
            'AGCGCCTAGCTTTGGCACTTCAGGAGCGCCCGGAAATAATGGCGGGCGATGAAGGTTCTG',
            'TAGGTACGCAAGATCCCTCTTAATCACAGTGGTGTAATCTGCGGGTCAGACCCTGTTAAC',
            'CCGTGGCTTTCACACTCCCTCCTATGGGTAATCAATCCAGAAAGGGGCCGAAATGCAAAA',
            'GTCTTAAGGACTCTGCGAGGCAAAGTACGGGCGAACTAAACCCCCGTGACAGGTCAGACG',
            'TTGTTTCGGCAATCTGTCGCGCTCCCACACCTATAAGCGTACACCGTCTCTTCTGCCAGC',
        ]))

        tmp_file = 'tmp.best_seq.fa'
        chooser = best_seq_chooser.BestSeqChooser(
            reads1,
            reads2,
            ref,
            sys.stdout,
            samtools_exe=extern_progs.exe('samtools'),
            bowtie2_exe=extern_progs.exe('bowtie2'),
        )
        got_seq = chooser.best_seq(tmp_file)
        self.assertEqual(expected_seq, got_seq)
        os.unlink(tmp_file)
