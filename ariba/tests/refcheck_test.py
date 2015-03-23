import unittest
import os
import filecmp
import pyfastaq
from ariba import refcheck

modules_dir = os.path.dirname(os.path.abspath(refcheck.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


class TestRefcheck(unittest.TestCase):
    def test_check_pass(self):
        '''test check file OK'''
        infile = os.path.join(data_dir, 'refcheck_test_check_ok.fa')
        c = refcheck.Checker(infile)
        self.assertEqual(c.check(), (True, None, None))


    def test_check_file_fail_not_gene(self):
        '''test check file fail not a gene'''
        infile = os.path.join(data_dir, 'refcheck_test_check_not_gene.fa')
        c = refcheck.Checker(infile)
        seq = pyfastaq.sequences.Fasta('gene1', 'TTGTGATGA')
        self.assertEqual(c.check(), (False, 'Not a gene', seq))


    def test_check_file_fail_too_short(self):
        '''test check file fail short gene'''
        infile = os.path.join(data_dir, 'refcheck_test_check_too_short.fa')
        c = refcheck.Checker(infile, min_length=10)
        seq = pyfastaq.sequences.Fasta('gene1', 'TTGTGGTGA')
        self.assertEqual(c.check(), (False, 'Too short', seq))


    def test_check_fix(self):
        '''test fix'''
        infile = os.path.join(data_dir, 'refcheck_test_fix_in.fa')
        tmp_prefix = 'tmp.refcheck_test_fix.out'
        c = refcheck.Checker(infile, min_length=10)
        c.fix(tmp_prefix)
        for x in ['fa', 'log', 'rename']:
            expected = os.path.join(data_dir, 'refcheck_test_fix_out.' + x)
            got = tmp_prefix + '.' + x
            self.assertTrue(filecmp.cmp(expected, got, shallow=False))
            os.unlink(got)
