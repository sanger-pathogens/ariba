import unittest
import sys
import os
import filecmp
import pymummer
from ariba import ref_seq_chooser

modules_dir = os.path.dirname(os.path.abspath(ref_seq_chooser.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


class TestRefSeqChooser(unittest.TestCase):
    def test_make_matching_contig_pieces_fasta(self):
        '''test _make_matching_contig_pieces_fasta'''
        contigs_fasta =  os.path.join(data_dir, 'ref_seq_chooser_matching_contig_pieces.ctg.fa')
        coords_file = os.path.join(data_dir, 'ref_seq_chooser_matching_contig_pieces.coords')
        coords_list = [hit for hit in pymummer.coords_file.reader(coords_file)]
        tmp_file = 'tmp.test.ref_seq_chooser_matching_contig_pieces.fa'
        expected_file = os.path.join(data_dir, 'ref_seq_chooser_matching_contig_pieces.expect.fa')
        ref_seq_chooser.RefSeqChooser._make_matching_contig_pieces_fasta(contigs_fasta, coords_list, tmp_file)
        self.assertTrue(filecmp.cmp(expected_file, tmp_file, shallow=False))
        os.unlink(tmp_file)


    def _sequence_is_in_fasta_file(self):
        '''Test _sequence_is_in_fasta_file'''
        fasta = os.path.join(data_dir, 'ref_seq_chooser_sequence_is_in_fasta_file.fa')
        self.assertFalse(ref_seq_chooser.RefSeqChooser._sequence_is_in_fasta_file('not there', fasta))
        self.assertTrue(ref_seq_chooser.RefSeqChooser._sequence_is_in_fasta_file('contig42', fasta))


    def test_run_no_match_in_cluster(self):
        '''Test full run when nearest match is not in the cluster'''
        all_ref_fasta = os.path.join(data_dir, 'ref_seq_chooser_full_run_not_in_cluster.allrefs.fa')
        cluster_fasta = os.path.join(data_dir, 'ref_seq_chooser_full_run_not_in_cluster.clusterrefs.fa')
        contig_fasta = os.path.join(data_dir, 'ref_seq_chooser_full_run_not_in_cluster.contigs.fa')
        tmp_out = 'tmp.ref_seq_chooser_full_run_not_in_cluster.fa'
        refchooser = ref_seq_chooser.RefSeqChooser(cluster_fasta, all_ref_fasta, contig_fasta, tmp_out, sys.stdout)
        refchooser.run()
        self.assertEqual(None, refchooser.closest_ref_from_all_refs)
        self.assertFalse(refchooser.closest_ref_is_in_cluster)
        self.assertFalse(os.path.exists(tmp_out))


    def test_run_no_nucmer_match(self):
        '''Test full run when there is nearest match in cluster, but no nucmer matches'''
        all_ref_fasta = os.path.join(data_dir, 'ref_seq_chooser_full_run_no_nucmer_match.allrefs.fa')
        cluster_fasta = os.path.join(data_dir, 'ref_seq_chooser_full_run_no_nucmer_match.clusterrefs.fa')
        contig_fasta = os.path.join(data_dir, 'ref_seq_chooser_full_run_no_nucmer_match.contigs.fa')
        tmp_out = 'tmp.ref_seq_chooser_full_run_no_nucmer_match.fa'
        refchooser = ref_seq_chooser.RefSeqChooser(cluster_fasta, all_ref_fasta, contig_fasta, tmp_out, sys.stdout)
        refchooser.run()
        self.assertEqual(None, refchooser.closest_ref_from_all_refs)
        self.assertFalse(refchooser.closest_ref_is_in_cluster)
        self.assertFalse(os.path.exists(tmp_out))


    def test_run_best_match_not_in_cluster(self):
        '''Test full run where there is a match in cluster, but better match to seq not in cluster'''
        all_ref_fasta = os.path.join(data_dir, 'ref_seq_chooser_full_run_best_match_not_in_cluster.allrefs.fa')
        cluster_fasta = os.path.join(data_dir, 'ref_seq_chooser_full_run_best_match_not_in_cluster.clusterrefs.fa')
        contig_fasta = os.path.join(data_dir, 'ref_seq_chooser_full_run_best_match_not_in_cluster.contigs.fa')
        tmp_out = 'tmp.ref_seq_chooser_full_run_best_match_not_in_cluster.fa'
        refchooser = ref_seq_chooser.RefSeqChooser(cluster_fasta, all_ref_fasta, contig_fasta, tmp_out, sys.stdout)
        refchooser.run()
        self.assertEqual('ref2', refchooser.closest_ref_from_all_refs)
        self.assertFalse(refchooser.closest_ref_is_in_cluster)
        self.assertFalse(os.path.exists(tmp_out))


    def test_run_best_match_is_in_cluster(self):
        '''Test full run where the best match is in the cluster'''
        all_ref_fasta = os.path.join(data_dir, 'ref_seq_chooser_full_run_best_match_is_in_cluster.allrefs.fa')
        cluster_fasta = os.path.join(data_dir, 'ref_seq_chooser_full_run_best_match_is_in_cluster.clusterrefs.fa')
        contig_fasta = os.path.join(data_dir, 'ref_seq_chooser_full_run_best_match_is_in_cluster.contigs.fa')
        tmp_out = 'tmp.ref_seq_chooser_full_run_best_match_is_in_cluster.fa'
        refchooser = ref_seq_chooser.RefSeqChooser(cluster_fasta, all_ref_fasta, contig_fasta, tmp_out, sys.stdout)
        refchooser.run()
        self.assertEqual('ref1', refchooser.closest_ref_from_all_refs)
        self.assertTrue(refchooser.closest_ref_is_in_cluster)
        self.assertTrue(os.path.exists(tmp_out))
        os.unlink(tmp_out)


    def test_run_contained_ref_seq(self):
        '''Test full run where ref seq completely contains another seq outside cluster'''
        all_ref_fasta = os.path.join(data_dir, 'ref_seq_chooser_full_run_contained_ref_seq.all_refs.fa')
        cluster_fasta = os.path.join(data_dir, 'ref_seq_chooser_full_run_contained_ref_seq.cluster_refs.fa')
        contig_fasta = os.path.join(data_dir, 'ref_seq_chooser_full_run_contained_ref_seq.contigs.fa')
        tmp_out = 'tmp.ref_seq_chooser_full_run_contained_ref_seq.fa'
        refchooser = ref_seq_chooser.RefSeqChooser(cluster_fasta, all_ref_fasta, contig_fasta, tmp_out, sys.stdout)
        refchooser.run()
        self.assertEqual('ref2', refchooser.closest_ref_from_all_refs)
        self.assertTrue(refchooser.closest_ref_is_in_cluster)
        self.assertTrue(os.path.exists(tmp_out))
        os.unlink(tmp_out)


    def test_run_flanking_different(self):
        '''Test full run where amount of flanking seq varies'''
        all_ref_fasta = os.path.join(data_dir, 'ref_seq_chooser_test_flanking.all_refs.fa')
        cluster_fasta = os.path.join(data_dir, 'ref_seq_chooser_test_flanking.cluster_refs.fa')
        contig_fasta = os.path.join(data_dir, 'ref_seq_chooser_test_flanking.contigs.fa')
        expected_fa = os.path.join(data_dir, 'ref_seq_chooser_test_flanking.expected_contigs.fa')
        tmp_out = 'tmp.ref_seq_chooser_test_flanking.fa'
        refchooser = ref_seq_chooser.RefSeqChooser(cluster_fasta, all_ref_fasta, contig_fasta, tmp_out, sys.stdout)
        refchooser.run()
        self.assertEqual('ref1', refchooser.closest_ref_from_all_refs)
        self.assertTrue(refchooser.closest_ref_is_in_cluster)
        self.assertTrue(filecmp.cmp(expected_fa, tmp_out, shallow=False))
        os.unlink(tmp_out)
