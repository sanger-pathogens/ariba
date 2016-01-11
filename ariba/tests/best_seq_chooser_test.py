import unittest
from ariba import best_seq_chooser

modules_dir = os.path.dirname(os.path.abspath(best_seq_chooser.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')
cluster.unittest = True


class TestBestSeqChooser(unittest.TestCase):
    def test_total_alignment_score(self):
        '''test _total_alignment_score'''
        cluster_dir = os.path.join(data_dir, 'cluster_test_get_total_alignment_score')
        clean_cluster_dir(cluster_dir)
        c = cluster.Cluster(cluster_dir, 'name')
        got_score = c._get_total_alignment_score('1')
        expected_score = 3000
        self.assertEqual(got_score, expected_score)
        clean_cluster_dir(cluster_dir)


    def test_get_best_seq_by_alignment_score(self):
        '''test _get_best_seq_by_alignment_score'''
        cluster_dir = os.path.join(data_dir, 'cluster_test_get_best_seq_by_alignment_score')
        clean_cluster_dir(cluster_dir)
        c = cluster.Cluster(cluster_dir, 'name')
        got_name = c._get_best_seq_by_alignment_score()
        self.assertEqual(got_name, '1')
        clean_cluster_dir(cluster_dir)


    def test_best_seq(self):
        '''test best_seq'''
        cluster_dir = os.path.join(data_dir, 'cluster_test_choose_best_seq')
        clean_cluster_dir(cluster_dir)
        c = cluster.Cluster(cluster_dir, 'name')
        expected_seq = pyfastaq.sequences.Fasta('1', ''.join([
            'AGCGCCTAGCTTTGGCACTTCAGGAGCGCCCGGAAATAATGGCGGGCGATGAAGGTTCTG',
            'TAGGTACGCAAGATCCCTCTTAATCACAGTGGTGTAATCTGCGGGTCAGACCCTGTTAAC',
            'CCGTGGCTTTCACACTCCCTCCTATGGGTAATCAATCCAGAAAGGGGCCGAAATGCAAAA',
            'GTCTTAAGGACTCTGCGAGGCAAAGTACGGGCGAACTAAACCCCCGTGACAGGTCAGACG',
            'TTGTTTCGGCAATCTGTCGCGCTCCCACACCTATAAGCGTACACCGTCTCTTCTGCCAGC',
        ]))
        expected_seq_fa = os.path.join(data_dir, 'cluster_test_choose_best_seq.seq.fa')
        got = c._choose_best_seq()
        self.assertEqual(got, expected_seq)
        self.assertTrue(filecmp.cmp(expected_seq_fa, c.seq_fa, shallow=False))
        clean_cluster_dir(cluster_dir)


