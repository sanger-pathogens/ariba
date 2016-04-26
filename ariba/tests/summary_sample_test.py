import unittest
import os
from ariba import summary_cluster, summary_sample

modules_dir = os.path.dirname(os.path.abspath(summary_sample.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


class TestSummarySample(unittest.TestCase):
    def test_load_file(self):
        '''Test _load_file'''
        infile = os.path.join(data_dir, 'summary_sample_test_load_file.in.tsv')
        with open(infile) as f:
             lines = [x.rstrip() for x in f]

        dicts = [summary_cluster.SummaryCluster.line2dict(x) for x in lines[1:]]
        cluster1 = summary_cluster.SummaryCluster()
        cluster1.add_data_dict(dicts[0])
        cluster1.add_data_dict(dicts[1])
        cluster1.add_data_dict(dicts[2])
        cluster2 = summary_cluster.SummaryCluster()
        cluster2.add_data_dict(dicts[3])
        cluster2.add_data_dict(dicts[4])
        cluster3 = summary_cluster.SummaryCluster()
        cluster3.add_data_dict(dicts[5])

        expected = {
            'cluster.n': cluster1,
            'cluster.p': cluster2,
            'cluster.v': cluster3
        }

        got = summary_sample.SummarySample._load_file(infile, 90)
        self.assertEqual(expected, got)
