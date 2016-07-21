import unittest
import os
from ariba import refdata_query

modules_dir = os.path.dirname(os.path.abspath(refdata_query.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


class TestRefdataQuery(unittest.TestCase):
    def setUp(self):
        self.prepareref_dir = os.path.join(data_dir, 'refdata_query_prepareref')
        self.rquery = refdata_query.RefdataQuery(self.prepareref_dir)


    def test_query_with_unknown_query(self):
        '''test query with unknown query'''
        with self.assertRaises(refdata_query.Error):
            self.rquery.query('notaquery', 'spam')


    def test_load_clusters(self):
        '''test _load_clusters'''
        infile = os.path.join(self.prepareref_dir, '02.cdhit.clusters.pickle')
        got = refdata_query.RefdataQuery._load_clusters(infile)
        expected = {
            '0': {'noncoding1', 'noncoding2'},
            '1': {'noncoding3'},
            '2': {'noncoding4.varonly'},
            '3': {'gene4'},
            '4': {'gene1_foo_', 'gene2', 'gene3'},
            '5': {'gene5.varonly', 'gene6.varonly'},
        }
        self.assertEqual(expected, got)


    def test_seq2cluster(self):
        '''test _seq2cluster'''
        clusters = {
            '0': {'seq1', 'seq2'},
            '1': {'seq3', 'seq4'},
        }
        self.assertEqual(None, refdata_query.RefdataQuery._seq2cluster(clusters, 'not_there'))
        self.assertEqual('0', refdata_query.RefdataQuery._seq2cluster(clusters, 'seq1'))
        self.assertEqual('0', refdata_query.RefdataQuery._seq2cluster(clusters, 'seq2'))
        self.assertEqual('1', refdata_query.RefdataQuery._seq2cluster(clusters, 'seq3'))
        self.assertEqual('1', refdata_query.RefdataQuery._seq2cluster(clusters, 'seq4'))


    def test_cluster2seqs(self):
        '''test _cluster2seqs'''
        expected = ['Sequences belonging to cluster 0:', 'noncoding1', 'noncoding2']
        got = self.rquery._cluster2seqs('0')
        self.assertEqual(expected, got)

        expected = ['Cluster name "fortytwo" not found']
        got = self.rquery._cluster2seqs('fortytwo')
        self.assertEqual(expected, got)


    def test_seqinfo(self):
        '''test _seqinfo'''
        expected = ['Sequence "fortytwo" not found']
        got = self.rquery._seqinfo('fortytwo')
        self.assertEqual(expected, got)

        expected = [
            'Name\tgene5.varonly',
            'Is gene\t1',
            'Variant only\t1',
            'Cluster\t5',
            'Description\tGeneric description of gene5.varonly',
            'Variant\tH7I\t.\tGeneric description of gene5.varonly H7I',
            'Variant\tH7J\t.\tGeneric description of gene5.varonly H7J',
            'Sequence\tATGCTGCAATCACTCAACCATCTGACCCTCGCGGTCAGCGACCTGCAAAAAAGCGTTACCTTCTGGCACGAGCTGCTGGGGCTGACGCTGCACGCCCGCTGGAATACCGGGGCCTATCTTACCTGCGGCGATCTGTGGGTCTGCCTGTCCTATGACGAGGCGCGCGGTTACGTGCCGCCGCAGGAGAGCGACTATACCCATTACGCGTTTACCGTTGCGGCGGAAGATTTTGAGCCGTTCTCGCACAAGCTGGAGCAGGCGGGCGTTACCGTCTGGAAGCAAAACAAAAGTGAGGGGGCATCGTTCTATTTTCTCGACCCGGACGGGCACAAGCTGGAGCTGCACGTGGGCAGCCTCGCCGCGCGGCTGGCGGCGTGCCGGGAGAAACCCTATGCCGGAATGGTCTTCACCTCAGACGAGGCTTGA',
        ]

        got = self.rquery._seqinfo('gene5.varonly')
        self.assertEqual(expected, got)

        expected = [
            'Name\tnoncoding1',
            'Is gene\t0',
            'Variant only\t0',
            'Cluster\t0',
            'Description\tGeneric description of noncoding1',
            'Sequence\tTCGACGTCGTCAGTACGTCACGTACGTACGTACGTACGTAGTCAGTCAGTCAGTCGTTAATACCTACTGACTGACTGATCGACGTACGTCTGACTGAGTCAGTCAGTCAGTCAGTCTAAACATCTACTGATGTACGAGCAGCTACGTACCGTACGTACGTACGTACTCATATCGTCGTTAACACATACTCTGATCGATCACGTAGTCGTCGTACGTACGTAGTCATATAT',
        ]
        got = self.rquery._seqinfo('noncoding1')
        self.assertEqual(expected, got)

