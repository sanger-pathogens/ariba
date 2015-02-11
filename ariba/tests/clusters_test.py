import unittest
import shutil
import os
import pysam
import pyfastaq
import filecmp
from ariba import clusters

modules_dir = os.path.dirname(os.path.abspath(clusters.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


class TestClusters(unittest.TestCase):
    def setUp(self):
        self.cluster_dir = 'tmp.Cluster'
        reads1 = os.path.join(data_dir, 'clusters_test_dummy_reads_1.fq')
        reads2 = os.path.join(data_dir, 'clusters_test_dummy_reads_2.fq')
        db = os.path.join(data_dir, 'clusters_test_dummy_db.fa')
        self.clusters = clusters.Clusters(db, reads1, reads2, self.cluster_dir)


    def tearDown(self):
        shutil.rmtree(self.cluster_dir)


    def test_sam_to_fastq(self):
        '''test _sam_to_fastq'''
        expected = [
            pyfastaq.sequences.Fastq('read1/1', 'GTATGAGTAGATATAAAGTCCGGAACTGTGATCGGGGGCGATTTATTTACTGGCCGTCCC', 'GHIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII'),
            pyfastaq.sequences.Fastq('read1/2', 'TCCCATACGTTGCAATCTGCAGACGCCACTCTTCCACGTCGGACGAACGCAACGTCAGGA', 'IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIHGEDCBA')
        ]


        sam_reader = pysam.Samfile(os.path.join(data_dir, 'clusters_test_sam_to_fastq.bam'), "rb")
        i = 0
        for s in sam_reader.fetch(until_eof=True):
            self.assertEqual(expected[i], self.clusters._sam_to_fastq(s))
            i += 1


    def test_sam_pair_to_insert(self):
        '''test _sam_pair_to_insert'''
        expected = [
            None, # both unmapped
            None, # read 1 unmapped
            None, # read 2 unmpapped
            None, # mapped to different seqs
            None, # same seqs, wrond orientation
            660
        ]

        sam1 = None
        i = 0
        sam_reader = pysam.Samfile(os.path.join(data_dir, 'clusters_test_sam_pair_to_insert.bam'), 'rb')
        for s in sam_reader.fetch(until_eof=True):
            if sam1 is None:
                sam1 = s
                continue

            self.assertEqual(self.clusters._sam_pair_to_insert(s, sam1), expected[i])
            sam1 = None
            i += 1


    def test_bam_to_clusters_reads(self):
        '''test _bam_to_clusters_reads'''
        clusters_dir = 'tmp.Cluster.test_bam_to_clusters_reads'
        reads1 = os.path.join(data_dir, 'clusters_test_bam_to_clusters_reads.reads_1.fq')    
        reads2 = os.path.join(data_dir, 'clusters_test_bam_to_clusters_reads.reads_2.fq')    
        ref = os.path.join(data_dir, 'clusters_test_bam_to_clusters_reads.db.fa')
        c = clusters.Clusters(ref, reads1, reads2, clusters_dir)
        shutil.copyfile(os.path.join(data_dir, 'clusters_test_bam_to_clusters_reads.bam'), c.bam)
        c._bam_to_clusters_reads()
        expected = [
            os.path.join(data_dir, 'clusters_test_bam_to_clusters.out.ref1.reads_1.fq'),
            os.path.join(data_dir, 'clusters_test_bam_to_clusters.out.ref1.reads_2.fq'),
            os.path.join(data_dir, 'clusters_test_bam_to_clusters.out.ref2.reads_1.fq'),
            os.path.join(data_dir, 'clusters_test_bam_to_clusters.out.ref2.reads_2.fq'),
        ]

        got = [
            os.path.join(clusters_dir, 'ref1/reads_1.fq'),
            os.path.join(clusters_dir, 'ref1/reads_2.fq'),
            os.path.join(clusters_dir, 'ref2/reads_1.fq'),
            os.path.join(clusters_dir, 'ref2/reads_2.fq'),
        ]


        for i in range(len(got)):
            self.assertTrue(os.path.exists(got[i]))
            self.assertTrue(filecmp.cmp(expected[i], got[i], shallow=False))

        self.assertEqual({780:1}, c.insert_hist.bins)

        shutil.rmtree(clusters_dir)


    def test_set_insert_size_data(self):
        '''test _set_insert_size_data'''
        self.clusters.insert_hist.bins = {
            1: 1,
            2: 1,
            3: 3,
            4: 3,
            5: 5,
            6: 3,
            7: 2,
            8: 2,
            9: 1,
            10: 1,
        }
        self.clusters.insert_hist.bin_width=1
       
        self.clusters._set_insert_size_data()
        self.assertEqual(self.clusters.insert_size, 5.5)
        self.assertEqual(self.clusters.insert_sspace_sd, 0.91)


    def test_write_gene_fa(self):
        '''Test _write_gene_fa'''
        self.clusters.db_fasta = os.path.join(data_dir, 'clusters_test_write_gene_fa.db.fa')
        expected = os.path.join(data_dir, 'clusters_test_write_gene_fa.out.fa')
        tmp_file = 'tmp.test_write_gene_fa.fa'
        self.clusters._write_gene_fa('gene2', tmp_file)
        self.assertTrue(filecmp.cmp(expected, tmp_file, shallow=False))
        os.unlink(tmp_file)


    def test_write_report(self):
        class FakeCluster:
            def __init__(self, lines):
                self.report_lines = lines

        self.clusters.clusters = {
            'gene1': FakeCluster(['gene1 line1']),
            'gene2': FakeCluster(['gene2 line2'])
        }

        self.clusters._write_report()
        expected = os.path.join(data_dir, 'clusters_test_write_report.tsv')
        self.assertTrue(filecmp.cmp(expected, self.clusters.report_file, shallow=False))

