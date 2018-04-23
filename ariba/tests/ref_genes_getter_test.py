import unittest
import os
import filecmp
from ariba import ref_genes_getter

modules_dir = os.path.dirname(os.path.abspath(ref_genes_getter.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


class TestRefGenesGetter(unittest.TestCase):
    def test_fix_virulencefinder_fasta_file(self):
        '''test _fix_virulencefinder_fasta_file'''
        infile = os.path.join(data_dir, 'ref_genes_getter.fix_virulencefinder_fasta_file.in.fa')
        tmp_file = 'tmp.test.ref_genes_getter.fix_virulencefinder_fasta_file.out.fa'
        expected_file = os.path.join(data_dir, 'ref_genes_getter.fix_virulencefinder_fasta_file.out.fa')
        ref_genes_getter.RefGenesGetter._fix_virulencefinder_fasta_file(infile, tmp_file)
        self.assertTrue(filecmp.cmp(expected_file, tmp_file, shallow=False))
        os.unlink(tmp_file)
