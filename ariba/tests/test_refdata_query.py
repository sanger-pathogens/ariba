import unittest
import os
import copy
import shutil
import filecmp
import pyfastaq
from ariba import refdata_query

modules_dir = os.path.dirname(os.path.abspath(refdata_query.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


class TestRefdataQuery(unittest.TestCase):
    def setUp(self):
        prepareref_dir = os.path.join(data_dir, 'refdata_query_prepareref')
        self.rquery = refdata_query.RefdataQuery(prepareref_dir)


    def test_query_with_unknown_query(self):
        with self.assertRaises(refdata_query.Error):
            self.rquery.query('notaquery', 'spam')


    def test_cluster2seqs(self):
        '''test _cluster2seqs'''
        expected = ['Sequences belonging to cluster 0:', '1', '2']
        got = self.rquery._cluster2seqs('0')
        self.assertEqual(expected, got)

