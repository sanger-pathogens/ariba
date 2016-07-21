import unittest
import os
import filecmp
import pyfastaq
import pysam
from ariba import link, scaffold_graph

modules_dir = os.path.dirname(os.path.abspath(scaffold_graph.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


class TestScaffoldGraph(unittest.TestCase):
    def test_update_from_sam(self):
        '''test update_from_sam'''
        sam_reader = pysam.Samfile(os.path.join(data_dir, 'graph_test_update_from_sam.bam'), "rb")
        ref_lengths = {}
        pyfastaq.tasks.lengths_from_fai(os.path.join(data_dir, 'graph_test_update_from_sam.ref.fa.fai'), ref_lengths)
        graph = scaffold_graph.Graph(ref_lengths)

        for sam in sam_reader.fetch(until_eof=True):
            graph.update_from_sam(sam, sam_reader)

        self.assertEqual(len(graph.partial_links), 0)
        self.assertEqual(len(graph.links), 1)
        key = ('ref1', 'ref2')
        self.assertTrue(key in graph.links)
        expected_links = [link.Link(None, None, None, s='\t'.join(['ref1', '506', 'R', '300', 'ref2', '500', 'L', '359']))]
        self.assertListEqual(graph.links[key], expected_links)


    def test_make_graph(self):
        '''test _make_graph'''
        ref_lengths = {'ref1':100, 'ref2':200, 'ref3':300}
        g = scaffold_graph.Graph(ref_lengths)

        g.partial_links = {42:42}
        with self.assertRaises(scaffold_graph.Error):
            g._make_graph(1000)
        g.partial_links = {}

        g.links[('ref1', 'ref2')] = [
            link.Link(None, None, None, '\t'.join(['ref1', '100', 'R', '50', 'ref2', '200', 'L', '10']))
        ]

        g._make_graph(10)
        self.assertEqual(len(g.contig_links), 0)

        g._make_graph(1000)
        self.assertEqual(len(g.contig_links), 1)
        expected_contig_links = {
            ('ref1', 'ref2'): {'RL': 1}
        }
        self.assertDictEqual(g.contig_links, expected_contig_links)

        g.links[('ref1', 'ref2')].append(link.Link(None, None, None, '\t'.join(['ref1', '100', 'R', '50', 'ref2', '200', 'L', '10'])))
        g.links[('ref1', 'ref2')].append(link.Link(None, None, None, '\t'.join(['ref1', '100', 'L', '50', 'ref2', '200', 'R', '10'])))
        g._make_graph(1000)
        expected_contig_links = {
            ('ref1', 'ref2'): {'RL': 2, 'LR': 1}
        }
        self.assertDictEqual(g.contig_links, expected_contig_links)


    def test_remove_low_cov_links(self):
        '''test _remove_low_cov_links'''
        ref_lengths = {'ref1':100, 'ref2':200, 'ref3':300}
        g = scaffold_graph.Graph(ref_lengths)
        g.contig_links = {
            ('ref1', 'ref2'): {'LR': 42, 'RL':41},
            ('ref2', 'ref3'): {'LR': 43, 'RL':42},
            ('ref3', 'ref4'): {'LR': 1, 'RL':2}
        }

        g._remove_low_cov_links(42)
        expected = {
            ('ref1', 'ref2'): {'LR': 42},
            ('ref2', 'ref3'): {'LR': 43, 'RL':42}
        }

        self.assertDictEqual(expected, g.contig_links)


    def test_contig_graph_is_consistent(self):
        '''test _contig_graph_is_consistent'''
        ref_lengths = {'ref1':100, 'ref2':200, 'ref3':300}
        g = scaffold_graph.Graph(ref_lengths)

        good_links = [
            {('ref1', 'ref2'): {'RL': 42}},
            {('ref1', 'ref2'): {'RL': 42}, ('ref2', 'ref3'): {'RL': 42}},
        ]

        bad_links = [
            {('ref1', 'ref2'): {'LR': 42, 'RL':42}},
            {('ref1', 'ref2'): {'LR': 42}, ('ref1', 'ref3'): {'LR': 42}}
        ]


        for d in good_links:
            g.contig_links = d
            self.assertTrue(g._contig_graph_is_consistent())

        for d in bad_links:
            g.contig_links = d
            self.assertFalse(g._contig_graph_is_consistent())


    def test_write_all_links_to_file(self):
        '''test write_all_links_to_file'''
        ref_lengths = {'ref1':100, 'ref2':200, 'ref3':300}
        g = scaffold_graph.Graph(ref_lengths)
        g.links[('ref1', 'ref2')] = [
            link.Link(None, None, None, '\t'.join(['ref1', '100', 'R', '50', 'ref2', '200', 'L', '10']))
        ]


        g.links[('ref1', 'ref2')].append(link.Link(None, None, None, '\t'.join(['ref1', '100', 'R', '50', 'ref2', '200', 'L', '10'])))
        g.links[('ref1', 'ref2')].append(link.Link(None, None, None, '\t'.join(['ref1', '100', 'L', '50', 'ref2', '200', 'R', '10'])))


        g.links[('ref3', 'ref4')] = [
            link.Link(None, None, None, '\t'.join(['ref3', '100', 'R', '42', 'ref4', '200', 'L', '42']))
        ]

        tmp_file = 'tmp.contig_links'
        g.write_all_links_to_file(tmp_file)
        self.assertTrue(filecmp.cmp(os.path.join(data_dir, 'graph_test_write_all_links_to_file.out'), tmp_file, shallow=False))
        os.unlink(tmp_file)


