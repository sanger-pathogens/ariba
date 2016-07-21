import unittest
import os
from ariba import mash

modules_dir = os.path.dirname(os.path.abspath(mash.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


class TestMash(unittest.TestCase):
    def test_run(self):
        '''test run'''
        ref_in = os.path.join(data_dir, 'mash_test_run.in.ref.fa')
        qry_in = os.path.join(data_dir, 'mash_test_run.in.qry.fa')
        masher = mash.Masher(ref_in, qry_in)
        tmp_out = 'tmp.mash_test_run.out'
        got = masher.run(tmp_out)
        self.assertEqual('ref2', got)
        os.unlink(tmp_out)
        os.unlink(ref_in + '.msh')
        os.unlink(qry_in + '.msh')

