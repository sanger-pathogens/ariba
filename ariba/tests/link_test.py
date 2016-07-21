import unittest
import os
import copy
import pyfastaq
import pysam
from ariba import link

modules_dir = os.path.dirname(os.path.abspath(link.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


class TestLink(unittest.TestCase):
    def test_init_no_links(self):
        '''test link __init__ no links made'''
        sam_reader = pysam.Samfile(os.path.join(data_dir, 'link_test_init.reads.no_link.bam'), "rb")
        ref_lengths = {}
        pyfastaq.tasks.lengths_from_fai(os.path.join(data_dir, 'link_test_init.ref.fa.fai'), ref_lengths)
        for sam in sam_reader.fetch(until_eof=True):
            with self.assertRaises(link.Error):
                link.Link(sam, sam_reader, ref_lengths)


    def test_init_with_link(self):
        '''test link __init__ link made'''
        sam_reader = pysam.Samfile(os.path.join(data_dir, 'link_test_init.reads.make_link.bam'), "rb")
        ref_lengths = {}
        pyfastaq.tasks.lengths_from_fai(os.path.join(data_dir, 'link_test_init.ref.fa.fai'), ref_lengths)
        links_from_bam = []
        for sam in sam_reader.fetch(until_eof=True):
            links_from_bam.append(link.Link(sam, sam_reader, ref_lengths))

        expected = [
            '\t'.join(['1', '506', 'R', '300', '2', '500', 'L', '.']),
            '\t'.join(['1', '506', 'R', '.', '2', '500', 'L', '359']),
        ]

        self.assertEqual(len(expected), len(links_from_bam))

        for i in range(len(expected)):
            self.assertEqual(expected[i], str(links_from_bam[i]))


    def test_lt(self):
        '''test __lt__'''
        l = [
            (['1', '10', 'R', '5', '2', '20', 'L', '3'], ['1', '10', 'R', '5', '2', '20', 'L', '3'], False),
            (['1', '10', 'R', '5', '2', '20', 'L', '3'], ['1', '10', 'R', '5', '2', '20', 'L', '4'], True),
            (['1', '10', 'R', '5', '2', '20', 'L', '3'], ['1', '10', 'R', '6', '2', '20', 'L', '3'], True),
            (['1', '10', 'R', '5', '2', '20', 'L', '3'], ['2', '10', 'R', '6', '2', '20', 'L', '3'], True),
            (['2', '10', 'R', '5', '2', '20', 'L', '3'], ['1', '10', 'R', '6', '2', '20', 'L', '3'], False),
        ]

        for t in l:
            link1 = link.Link(None, None, None, '\t'.join(t[0]))
            link2 = link.Link(None, None, None, '\t'.join(t[1]))
            self.assertEqual(link1 < link2, t[2])


    def test_swap(self):
        '''test _swap'''
        original_str = '\t'.join(['1', '506', 'R', '300', '2', '500', 'L', '359'])
        original_link = link.Link(None, None, None, s=original_str)
        l = copy.copy(original_link)
        expected_swap = link.Link(None, None, None, s='\t'.join(['2', '500', 'L', '359', '1', '506', 'R', '300']))
        l._swap()
        self.assertEqual(l, expected_swap)
        l._swap()
        self.assertEqual(l, original_link)


    def test_sort(self):
        '''test sort'''
        links = [
            link.Link(None, None, None, s='\t'.join(['ref1', '500', 'L', '359', 'ref2', '506', 'R', '300'])),
            link.Link(None, None, None, s='\t'.join(['ref2', '500', 'L', '359', 'ref1', '506', 'R', '300']))
        ]

        expected = [
            link.Link(None, None, None, s='\t'.join(['ref1', '500', 'L', '359', 'ref2', '506', 'R', '300'])),
            link.Link(None, None, None, s='\t'.join(['ref1', '506', 'R', '300', 'ref2', '500', 'L', '359'])),
        ]

        assert len(links) == len(expected)
        for i in range(len(links)):
            links[i].sort()
            self.assertEqual(links[i], expected[i])


    def test_merge(self):
        '''test merge'''
        link1 = link.Link(None, None, None, '\t'.join(['1', '506', 'R', '300', '2', '500', 'L', '.']))
        link2 = link.Link(None, None, None, '\t'.join(['1', '506', 'R', '.', '2', '500', 'L', '359']))
        merged = link.Link(None, None, None, '\t'.join(['1', '506', 'R', '300', '2', '500', 'L', '359']))
        with self.assertRaises(link.Error):
            link1.merge(link1)
            link2.merge(link2)

        link1.merge(link2)
        self.assertEqual(link1, merged)


    def test_distance_to_contig_end(self):
        '''test _distance_to_contig_end'''
        links = [
            link.Link(None, None, None, '\t'.join(['1', '506', 'R', '300', '2', '500', 'L', '.'])),
            link.Link(None, None, None, '\t'.join(['1', '506', 'R', '.', '2', '500', 'L', '359']))
        ]

        self.assertEqual(links[0]._distance_to_contig_end(1), 205)
        self.assertEqual(links[1]._distance_to_contig_end(2), 359)
        with self.assertRaises(link.Error):
            links[0]._distance_to_contig_end(2)

        with self.assertRaises(link.Error):
            links[1]._distance_to_contig_end(1)


    def test_distance_to_contig_ends(self):
        '''test _distance_to_contig_ends'''
        link1 = link.Link(None, None, None, '\t'.join(['1', '506', 'R', '300', '2', '500', 'L', '.']))
        link2 = link.Link(None, None, None, '\t'.join(['1', '506', 'R', '.', '2', '500', 'L', '359']))
        link3 = link.Link(None, None, None, '\t'.join(['1', '506', 'R', '300', '2', '500', 'L', '359']))

        with self.assertRaises(link.Error):
            link1._distance_to_contig_ends()

        with self.assertRaises(link.Error):
            link2._distance_to_contig_ends()

        self.assertEqual(link3._distance_to_contig_ends(), (205, 359))


    def test_insert_size(self):
        '''test insert_size'''
        link1 = link.Link(None, None, None, '\t'.join(['1', '506', 'R', '300', '2', '500', 'L', '.']))
        link2 = link.Link(None, None, None, '\t'.join(['1', '506', 'R', '.', '2', '500', 'L', '359']))
        link3 = link.Link(None, None, None, '\t'.join(['1', '506', 'R', '300', '2', '500', 'L', '359']))

        with self.assertRaises(link.Error):
            link1.insert_size()

        with self.assertRaises(link.Error):
            link2.insert_size()

        self.assertEqual(link3.insert_size(), 564)

