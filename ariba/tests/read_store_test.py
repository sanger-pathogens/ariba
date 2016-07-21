import unittest
import os
import shutil
import filecmp
import pyfastaq
from ariba import read_store

modules_dir = os.path.dirname(os.path.abspath(read_store.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


def file_to_list(infile):
    f = pyfastaq.utils.open_file_read(infile)
    lines = [x for x in f.readlines()]
    pyfastaq.utils.close(f)
    return lines


class TestReadStore(unittest.TestCase):
    def test_sort_file(self):
        '''test _sort_file'''
        infile = os.path.join(data_dir, 'read_store_test_sort_file.in')
        expected = os.path.join(data_dir, 'read_store_test_sort_file.out')
        tmpfile = 'tmp.read_store_test_sort_file.out'
        read_store.ReadStore._sort_file(infile, tmpfile)
        self.assertTrue(filecmp.cmp(expected, tmpfile, shallow=False))
        os.unlink(tmpfile)


    def test_compress_and_index_file(self):
        '''Test _compress_and_index_file'''
        infile = os.path.join(data_dir, 'read_store_test_compress_and_index_file.in')
        tmpfile = 'tmp.test_compress_and_index_file.in'
        tmpfile_gz = 'tmp.test_compress_and_index_file.in.gz'
        shutil.copyfile(infile, tmpfile)
        read_store.ReadStore._compress_and_index_file(tmpfile)
        self.assertTrue(os.path.exists(tmpfile_gz))
        expected_lines = file_to_list(infile)
        got_lines = file_to_list(tmpfile_gz)
        self.assertEqual(expected_lines, got_lines)
        self.assertTrue(os.path.exists(tmpfile_gz + '.tbi'))
        os.unlink(tmpfile)
        os.unlink(tmpfile_gz)
        os.unlink(tmpfile_gz + '.tbi')


    def test_get_reads_fq_pair(self):
        '''Test get_reads fastq pair'''
        infile = os.path.join(data_dir, 'read_store_test_get_reads.in')
        expected1 = os.path.join(data_dir, 'read_store_test_get_reads.expected.reads_1.fq')
        expected2 = os.path.join(data_dir, 'read_store_test_get_reads.expected.reads_2.fq')
        outprefix = 'tmp.read_store_test_get_reads'
        reads1 = outprefix + '.reads_1.fq'
        reads2 = outprefix + '.reads_2.fq'
        rstore = read_store.ReadStore(infile, outprefix)
        got_reads, got_bases = rstore.get_reads('cluster2', reads1, out2=reads2)
        self.assertEqual(6, got_reads)
        self.assertEqual(24, got_bases)
        self.assertTrue(filecmp.cmp(expected1, reads1))
        self.assertTrue(filecmp.cmp(expected2, reads2))
        os.unlink(outprefix + '.gz')
        os.unlink(outprefix + '.gz.tbi')
        os.unlink(reads1)
        os.unlink(reads2)


    def test_get_reads_fq_interleave(self):
        '''Test get_reads fastq interleaved'''
        infile = os.path.join(data_dir, 'read_store_test_get_reads.in')
        expected = os.path.join(data_dir, 'read_store_test_get_reads.expected.reads.fq')
        outprefix = 'tmp.read_store_test_get_reads'
        reads = outprefix + '.reads_1.fq'
        rstore = read_store.ReadStore(infile, outprefix)
        got_reads, got_bases = rstore.get_reads('cluster2', reads)
        self.assertEqual(6, got_reads)
        self.assertEqual(24, got_bases)
        self.assertTrue(filecmp.cmp(expected, reads))
        os.unlink(outprefix + '.gz')
        os.unlink(outprefix + '.gz.tbi')
        os.unlink(reads)


    def test_get_reads_fa_pair(self):
        '''Test get_reads fasta pair'''
        infile = os.path.join(data_dir, 'read_store_test_get_reads.in')
        expected1 = os.path.join(data_dir, 'read_store_test_get_reads.expected.reads_1.fa')
        expected2 = os.path.join(data_dir, 'read_store_test_get_reads.expected.reads_2.fa')
        outprefix = 'tmp.read_store_test_get_reads'
        reads1 = outprefix + '.reads_1.fa'
        reads2 = outprefix + '.reads_2.fa'
        rstore = read_store.ReadStore(infile, outprefix)
        got_reads, got_bases = rstore.get_reads('cluster2', reads1, out2=reads2, fasta=True)
        self.assertEqual(6, got_reads)
        self.assertEqual(24, got_bases)
        self.assertTrue(filecmp.cmp(expected1, reads1))
        self.assertTrue(filecmp.cmp(expected2, reads2))
        os.unlink(outprefix + '.gz')
        os.unlink(outprefix + '.gz.tbi')
        os.unlink(reads1)
        os.unlink(reads2)


    def test_get_reads_subset(self):
        '''Test get_reads subset'''
        infile = os.path.join(data_dir, 'read_store_test_get_reads.in')
        expected1 = os.path.join(data_dir, 'read_store_test_get_reads.expected.reads_subset.1.fq')
        expected2 = os.path.join(data_dir, 'read_store_test_get_reads.expected.reads_subset.2.fq')
        wanted_ids = {1, 11}
        outprefix = 'tmp.read_store_test_get_reads'
        reads1 = outprefix + '.reads_1.fq'
        reads2 = outprefix + '.reads_2.fq'
        rstore = read_store.ReadStore(infile, outprefix)
        got_reads, got_bases = rstore.get_reads('cluster2', reads1, out2=reads2, wanted_ids=wanted_ids)
        self.assertEqual(4, got_reads)
        self.assertEqual(16, got_bases)
        self.assertTrue(filecmp.cmp(expected1, reads1))
        self.assertTrue(filecmp.cmp(expected2, reads2))
        os.unlink(outprefix + '.gz')
        os.unlink(outprefix + '.gz.tbi')
        os.unlink(reads1)
        os.unlink(reads2)


    def test_clean(self):
        '''Test clean'''
        infile = os.path.join(data_dir, 'read_store_test_clean.in')
        outprefix = 'tmp.read_store_test_clean'
        self.assertFalse(os.path.exists(outprefix))
        self.assertFalse(os.path.exists(outprefix + '.gz'))
        self.assertFalse(os.path.exists(outprefix + '.gz.tbi'))
        rstore = read_store.ReadStore(infile, outprefix)
        self.assertFalse(os.path.exists(outprefix))
        self.assertTrue(os.path.exists(outprefix + '.gz'))
        self.assertTrue(os.path.exists(outprefix + '.gz.tbi'))
        rstore.clean()
        self.assertFalse(os.path.exists(outprefix))
        self.assertFalse(os.path.exists(outprefix + '.gz'))
        self.assertFalse(os.path.exists(outprefix + '.gz.tbi'))
