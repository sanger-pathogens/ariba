#!/usr/bin/env python3
import unittest
import os
from ariba.ref_genes_getter import RefGenesGetter

class TestNcbiGetter(unittest.TestCase):
    def setUp(self):
        self.ncbi_db = RefGenesGetter('ncbi')._get_from_ncbi('ncbi.test', 'test')
        # self.ncbi_db = RefGenesGetter.run('ncbi')
        
    def test_ncbi(self):
        '''
        Test that more than 4000 records have been found on NCBI AMR DB.
        '''
        self.assertTrue(len(self.ncbi_db) > 4000)

if __name__ == '__main__':
    unittest.main()