import unittest
import os
from ariba import histogram

modules_dir = os.path.dirname(os.path.abspath(histogram.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')

class TestCluster(unittest.TestCase):
    def test_to_bin(self):
        '''Test _to_bin'''
        h = histogram.Histogram(3)
        tests = [
            (0, 0),
            (1, 0),
            (2, 0),
            (3, 3),
            (4, 3),
            (5, 3),
            (6, 6),
        ]

        for t in tests:
            self.assertEqual(h._to_bin(t[0]), t[1])


    def test_len(self):
        '''Test __len__'''
        h = histogram.Histogram(1)
        self.assertEqual(0, len(h))
        h.add(1)
        self.assertEqual(1, len(h))
        h.add(1)
        self.assertEqual(2, len(h))
        h.add(2)
        self.assertEqual(3, len(h))


    def test_add(self):
        '''Test add'''
        h = histogram.Histogram(3)
        h.add(4)
        self.assertEqual({3:1}, h.bins)
        h.add(4)
        self.assertEqual({3:2}, h.bins)
        h.add(42)
        self.assertEqual({3:2, 42:1}, h.bins)
        h.add(42, count=11)
        self.assertEqual({3:2, 42:12}, h.bins)


    def test_stats(self):
        '''Test stats'''
        h = histogram.Histogram(1)
        for i in range(10):
            h.add(i+1)
        h.bins[3] = 3
        h.bins[4] = 3
        h.bins[5] = 5
        h.bins[6] = 3
        h.bins[7] = 2
        h.bins[8] = 2
        self.assertEqual((2.5, 5.5, 10.5, 0.91), h.stats())
