import unittest
import os
import pyfastaq
import pysam
from ariba import bam_parse, link

modules_dir = os.path.dirname(os.path.abspath(bam_parse.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


def file_to_lines(filename):
    f = pyfastaq.utils.open_file_read(filename)
    lines = [l.rstrip() for l in f]
    pyfastaq.utils.close(f)
    return lines


class TestBamParse(unittest.TestCase):
    def test_sam_to_soft_clipped(self):
        '''test _sam_to_soft_clipped'''
        bam = os.path.join(data_dir, 'bam_parse_test_sam_to_soft_clipped.bam')
        bp = bam_parse.Parser(bam, {'x':pyfastaq.sequences.Fasta('x', 'ACGT')})
        sam_reader = pysam.Samfile(bam, "rb")
        i = 0
        expected = [
            None,
            (False, False),
            (True, False),
            (False, True),
            (True, True)
        ]

        for sam in sam_reader.fetch(until_eof=True):
            if i == 0:
                with self.assertRaises(bam_parse.Error):
                    bp._sam_to_soft_clipped(sam)
            else:
                self.assertEqual(expected[i], bp._sam_to_soft_clipped(sam))

            i += 1

    def test_update_soft_clipped_from_sam(self):
        '''test _update_soft_clipped_from_sam'''
        ref_seqs = {}
        ref_fasta = os.path.join(data_dir, 'bam_parse_test_update_soft_clipped_from_sam.ref.fa')
        bam = os.path.join(data_dir, 'bam_parse_test_update_soft_clipped_from_sam.bam')
        pyfastaq.tasks.file_to_dict(ref_fasta, ref_seqs)
        bp = bam_parse.Parser(bam, ref_seqs)
        expected = [
            {},
            {'ref1': {61: [1, 0]}},
            {'ref1': {61: [2, 0], 118: [0, 1]}},
        ]
        i = 0
        sam_reader = pysam.Samfile(bam, "rb")
        for sam in sam_reader.fetch(until_eof=True):
            bp._update_soft_clipped_from_sam(sam)
            self.assertEqual(expected[i], bp.soft_clipped)
            i += 1


    def test_update_unmapped_mates_from_sam(self):
        '''test _update_unmapped_mates_from_sam'''
        ref_seqs = {}
        ref_fasta = os.path.join(data_dir, 'bam_parse_test_update_unmapped_mates_from_sam.ref.fa')
        bam = os.path.join(data_dir, 'bam_parse_test_update_unmapped_mates_from_sam.bam')
        pyfastaq.tasks.file_to_dict(ref_fasta, ref_seqs)
        bp = bam_parse.Parser(bam, ref_seqs)
        expected = [
            {},
            {'ref': {240: 1}},
            {'ref': {240: 1, 360: 1}},
            {'ref': {240: 1, 360: 1}},
            {'ref': {240: 1, 360: 1}},
            {'ref': {240: 1, 360: 1}},
            {'ref': {240: 1, 360: 1}},
            {'ref': {240: 1, 360: 1}},
        ]
        i = 0
        sam_reader = pysam.Samfile(bam, "rb")
        for sam in sam_reader.fetch(until_eof=True):
            bp._update_unmapped_mates_from_sam(sam)
            self.assertEqual(expected[i], bp.unmapped_mates)
            i += 1


    def test_parse(self):
        '''test parse'''
        ref_seqs = {}
        ref_fasta = os.path.join(data_dir, 'bam_parse_test_parse.ref.fa')
        bam = os.path.join(data_dir, 'bam_parse_test_parse.bam')
        pyfastaq.tasks.file_to_dict(ref_fasta, ref_seqs)
        bp = bam_parse.Parser(bam, ref_seqs)
        bp.parse()
        expected_soft_clipped = {'ref1': {2: [1, 0], 58: [0, 1]}}
        expected_unmapped_mates = {'ref1': {2: 1}}
        l = link.Link(None, None, None, s='\t'.join(['ref1', '500', 'R', '240', 'ref2', '500', 'L', '299']))
        expected_link_keys = [('ref1', 'ref2')]
        self.assertEqual(expected_soft_clipped, bp.soft_clipped)
        self.assertEqual(expected_unmapped_mates, bp.unmapped_mates)
        self.assertEqual(list(bp.scaff_graph.links.keys()), expected_link_keys)
        self.assertEqual(len(bp.scaff_graph.links[expected_link_keys[0]]), 1)
        self.assertEqual(bp.scaff_graph.links[expected_link_keys[0]][0], l)


    def test_write_soft_clipped_to_file(self):
        '''test _write_soft_clipped_to_file'''
        tmp_file = 'tmp.soft_clipped'
        bam = os.path.join(data_dir, 'dummy.bam')
        bp = bam_parse.Parser(bam, {'x':pyfastaq.sequences.Fasta('x', 'ACGT')})
        soft_clipped = {'ref1': {42: [142, 242], 43: [44, 45]}}
        bp.soft_clipped = soft_clipped
        bp._write_soft_clipped_to_file(tmp_file)
        expected_lines = [
            'ref1\t43\t142\t242',
            'ref1\t44\t44\t45'
        ]
        got_lines = file_to_lines(tmp_file)
        self.assertEqual(got_lines, expected_lines)
        os.unlink(tmp_file)


    def test_write_unmapped_mates_to_file(self):
        '''test _write_unmapped_mates_to_file'''
        tmp_file = 'tmp.unmapped_mates'
        bam = os.path.join(data_dir, 'dummy.bam')
        bp = bam_parse.Parser(bam, {'x':pyfastaq.sequences.Fasta('x', 'ACGT')})
        bp.unmapped_mates = {'ref': {41: 1, 42: 43}}
        bp._write_unmapped_mates_to_file(tmp_file)
        expected_lines = ['ref\t42\t1', 'ref\t43\t43']
        got_lines = file_to_lines(tmp_file)
        self.assertEqual(got_lines, expected_lines)
        os.unlink(tmp_file)
