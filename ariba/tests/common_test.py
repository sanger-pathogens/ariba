import unittest
import os
import filecmp
from ariba import common

modules_dir = os.path.dirname(os.path.abspath(common.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


class TestCommon(unittest.TestCase):
    def test_cat_files(self):
        '''test cat_files'''
        infiles = [
            os.path.join(data_dir, 'test_common_cat_files.in.1'),
            os.path.join(data_dir, 'test_common_cat_files.in.2'),
            os.path.join(data_dir, 'test_common_cat_files.in.3'),
        ]
        tmp_out = 'tmp.test.common_cat_files.out'
        expected = os.path.join(data_dir, 'test_common_cat_files.out')
        common.cat_files(infiles, tmp_out)
        self.assertTrue(filecmp.cmp(expected, tmp_out, shallow=False))
        os.unlink(tmp_out)


    def test_rmtree(self):
        '''test rmtree'''
        tmp_dir = 'tmp.rmtree'
        os.mkdir(tmp_dir)
        with open (os.path.join(tmp_dir, 'foo'), 'w') as f:
            pass

        self.assertTrue(os.path.exists(tmp_dir))
        common.rmtree(tmp_dir)
        self.assertFalse(os.path.exists(tmp_dir))


