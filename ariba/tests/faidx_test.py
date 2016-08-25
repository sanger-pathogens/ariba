import unittest
import filecmp
import os
from ariba import faidx, external_progs

modules_dir = os.path.dirname(os.path.abspath(faidx.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')
extern_progs = external_progs.ExternalProgs()


class TestFaidx(unittest.TestCase):
    def test_write_fa_subset(self):
        '''test write_fa_subset'''
        infile = os.path.join(data_dir, 'faidx_test_write_fa_subset.in.fa')
        expected = os.path.join(data_dir, 'faidx_test_write_fa_subset.out.fa')
        tmpfile = 'tmp.test_write_fa_subset.out.fa'
        faidx.write_fa_subset(['seq1', 'seq3', 'seq4'], infile, tmpfile)
        self.assertTrue(filecmp.cmp(expected, tmpfile, shallow=False))
        os.unlink(tmpfile)
