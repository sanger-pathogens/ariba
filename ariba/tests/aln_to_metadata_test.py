import unittest
import os
import filecmp
import pyfastaq
from ariba import aln_to_metadata, sequence_variant

modules_dir = os.path.dirname(os.path.abspath(aln_to_metadata.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


class TestAlnToMetadata(unittest.TestCase):
    def test_load_aln_file(self):
        '''test _load_aln_file'''
        aln_file = os.path.join(data_dir, 'aln_to_metadata_load_aln_file.in.fa')
        expected = {
            'seq1': pyfastaq.sequences.Fasta('seq1', 'ABC-DE'),
            'seq2': pyfastaq.sequences.Fasta('seq2', 'ABCQDE'),
        }
        got = aln_to_metadata.AlnToMetadata._load_aln_file(aln_file)
        self.assertEqual(expected, got)


    def test_load_vars_file_good_file(self):
        '''test _load_vars_file good input file'''
        infile = os.path.join(data_dir, 'aln_to_metadata_load_vars_file_good.tsv')
        variant1 = sequence_variant.Variant('p', 'A42B', 'id1')
        variant2 = sequence_variant.Variant('p', 'C43D', 'id2')
        variant3 = sequence_variant.Variant('p', 'E100F', 'id3')
        expected = {
            'seq1': [(variant1, 'description 1')],
            'seq2': [(variant2, 'description 2'), (variant3, 'description 3')]
        }
        got = aln_to_metadata.AlnToMetadata._load_vars_file(infile, True)
        self.assertEqual(expected, got)


    def test_load_vars_bad_files(self):
        '''test _load_vars_file bad input files'''
        infiles = [
            os.path.join(data_dir, 'aln_to_metadata_load_vars_file_bad.1.tsv'),
            os.path.join(data_dir, 'aln_to_metadata_load_vars_file_bad.2.tsv')
        ]

        for infile in infiles:
            with self.assertRaises(aln_to_metadata.Error):
                aln_to_metadata.AlnToMetadata._load_vars_file(infile, True)


    def test_make_unpadded_seqs(self):
        '''test _make_unpadded_seqs'''
        padded = {
            'seq1': pyfastaq.sequences.Fasta('seq1', 'acg---t'),
            'seq2': pyfastaq.sequences.Fasta('seq2', '---a-cgt-'),
        }
        expected = {
            'seq1': pyfastaq.sequences.Fasta('seq1', 'acgt'),
            'seq2': pyfastaq.sequences.Fasta('seq2', 'acgt'),
        }
        got = aln_to_metadata.AlnToMetadata._make_unpadded_seqs(padded)
        self.assertEqual(expected, got)


    def test_check_seq_lengths_same(self):
        '''test _check_seq_lengths_same'''
        seqs = {
            'seq1': pyfastaq.sequences.Fasta('seq1', 'acgt'),
            'seq2': pyfastaq.sequences.Fasta('seq2', 'acgt'),
        }

        self.assertTrue(aln_to_metadata.AlnToMetadata._check_seq_lengths_same(seqs))
        seqs['seq1'].seq = 'a'
        with self.assertRaises(aln_to_metadata.Error):
            aln_to_metadata.AlnToMetadata._check_seq_lengths_same(seqs)


    def test_insertion_coords(self):
        '''test _insertion_coords'''
        tests = [
            ('acgt', []),
            ('-a', [pyfastaq.intervals.Interval(0, 0)]),
            ('a---cgt--', [pyfastaq.intervals.Interval(1, 3), pyfastaq.intervals.Interval(7, 8)]),
        ]

        for seq, expected in tests:
            fa = pyfastaq.sequences.Fasta('x', seq)
            got = aln_to_metadata.AlnToMetadata._insertion_coords(fa)
            self.assertEqual(expected, got)


    def test_make_unpadded_insertion_coords(self):
        '''test _make_unpadded_insertion_coords'''
        seqs = {
            'seq1': pyfastaq.sequences.Fasta('seq1', 'acgt'),
            'seq2': pyfastaq.sequences.Fasta('seq2', 'ac-gt'),
            'seq3': pyfastaq.sequences.Fasta('seq3', '--acg-t'),
        }

        expected = {
            'seq1': [],
            'seq2': [pyfastaq.intervals.Interval(2, 2)],
            'seq3': [pyfastaq.intervals.Interval(0, 1), pyfastaq.intervals.Interval(5, 5)],

        }
        got = aln_to_metadata.AlnToMetadata._make_unpadded_insertion_coords(seqs)
        self.assertEqual(expected, got)


    def test_check_insertion_coords(self):
        '''test _check_insertion_coords'''
        seq = pyfastaq.sequences.Fasta('name', 'AAA---GGG------TTT---')
        self.assertTrue(aln_to_metadata.AlnToMetadata._check_insertion_coords(seq))

        bad_seqs = [
            pyfastaq.sequences.Fasta('name', 'AAA--GGG'),  # bad length
            pyfastaq.sequences.Fasta('name', 'A---AA'),  # bad start position
            pyfastaq.sequences.Fasta('name', 'AA---AA'), # bad start position
        ]

        for seq in bad_seqs:
            with self.assertRaises(aln_to_metadata.Error):
                aln_to_metadata.AlnToMetadata._check_insertion_coords(seq)


    def test_check_coding_seq(self):
        '''test _check_coding_seq'''
        seq = pyfastaq.sequences.Fasta('name', 'ATGCTTTAG')
        self.assertTrue(aln_to_metadata.AlnToMetadata._check_coding_seq(seq))

        bad_seqs = [
            pyfastaq.sequences.Fasta('name', 'TTGCTTAG'), # length not a mutliple of 3
            pyfastaq.sequences.Fasta('name', 'TTTCTTTAG'), # no start codon
            pyfastaq.sequences.Fasta('name', 'ATGTAGCTTTAG'), # stop codon in middle
            pyfastaq.sequences.Fasta('name', 'TTGCTTTTT'), # no stop at end
        ]

        for seq in bad_seqs:
            with self.assertRaises(aln_to_metadata.Error):
                aln_to_metadata.AlnToMetadata._check_coding_seq(seq)


    def test_check_sequences_non_coding(self):
        '''test _check_sequences with noncoding seqs'''
        padded_sequences = {
            'seq1': pyfastaq.sequences.Fasta('seq1', 'AC-T')
        }

        unpadded_sequences = aln_to_metadata.AlnToMetadata._make_unpadded_seqs(padded_sequences)
        self.assertTrue(aln_to_metadata.AlnToMetadata._check_sequences(padded_sequences, unpadded_sequences, False))
        padded_sequences['seq2'] = pyfastaq.sequences.Fasta('seq2', 'AC-')
        unpadded_sequences = aln_to_metadata.AlnToMetadata._make_unpadded_seqs(padded_sequences)
        with self.assertRaises(aln_to_metadata.Error):
            aln_to_metadata.AlnToMetadata._check_sequences(padded_sequences, unpadded_sequences, False)


    def test_check_sequences_coding(self):
        '''test _check_sequences with coding seqs'''
        padded_sequences = {
            'seq1': pyfastaq.sequences.Fasta('seq1', 'ATGCTTTAG'),
            'seq2': pyfastaq.sequences.Fasta('seq2', 'ATG---TAG')
        }

        unpadded_sequences = aln_to_metadata.AlnToMetadata._make_unpadded_seqs(padded_sequences)

        self.assertTrue(aln_to_metadata.AlnToMetadata._check_sequences(padded_sequences, unpadded_sequences, True))

        bad_seqs = [
            'ATGCTTAG', # length not a mutliple of 3
            'TTTCTTTAG', # no start codon
            'ATGTAGCTTTAG', # stop codon in middle
            'ATGTTTTTT', # no stop at end
            'ATGC---TTTAG', # bad insertion
            'ATGCT---TTAG', # bad insertion
            'ATG-CTTTAG', # bad insertion
            'ATG--CTTTAG', # bad insertion
            'ATG----CTTTAG', # bad insertion
        ]

        for seq in bad_seqs:
            padded_sequences['seq2'] = pyfastaq.sequences.Fasta('seq2', seq)
            unpadded_sequences = aln_to_metadata.AlnToMetadata._make_unpadded_seqs(padded_sequences)
            with self.assertRaises(aln_to_metadata.Error):
                aln_to_metadata.AlnToMetadata._check_sequences(padded_sequences, unpadded_sequences, True)


    def test_check_variants_match_sequences(self):
        '''test _check_variants_match_sequences'''
        seqs = {
            'seq1': pyfastaq.sequences.Fasta('seq1', 'ATGCTTTAG'),
            'seq2': pyfastaq.sequences.Fasta('seq2', 'ATGCTTCTTTAG'),
            'seq3': pyfastaq.sequences.Fasta('seq3', 'ATG---TAG')
        }

        variants = {'seq1': [(sequence_variant.Variant('p', 'L2M', 'id1'), 'description1')]}
        self.assertTrue(aln_to_metadata.AlnToMetadata._check_variants_match_sequences(seqs, variants, True))
        variants = {'seq1': [(sequence_variant.Variant('p', 'M2L', 'id1'), 'description1')]}
        self.assertTrue(aln_to_metadata.AlnToMetadata._check_variants_match_sequences(seqs, variants, True))

        variants = {'seq1': [(sequence_variant.Variant('p', 'A2M', 'id1'), 'description1')]}
        with self.assertRaises(aln_to_metadata.Error):
            self.assertTrue(aln_to_metadata.AlnToMetadata._check_variants_match_sequences(seqs, variants, True))

        variants = {'seq4': [(sequence_variant.Variant('p', 'A2M', 'id1'), 'description1')]}
        with self.assertRaises(aln_to_metadata.Error):
            self.assertTrue(aln_to_metadata.AlnToMetadata._check_variants_match_sequences(seqs, variants, True))


    def test_variant_ids_are_unique(self):
        '''test variant_ids_are_unique'''
        variants = {
            'seq1': [(sequence_variant.Variant('p', 'L2M', 'id1'), 'description1')],
            'seq2': [(sequence_variant.Variant('p', 'L2M', 'id2'), 'description2')]
        }

        self.assertTrue(aln_to_metadata.AlnToMetadata._variant_ids_are_unique(variants))
        variants['seq2'].append((sequence_variant.Variant('p', 'I3K', 'id1'), 'description3'))
        with self.assertRaises(aln_to_metadata.Error):
            self.assertTrue(aln_to_metadata.AlnToMetadata._variant_ids_are_unique(variants))


    def test_unpadded_to_padded_nt_position(self):
        '''test _unpadded_to_padded_nt_position'''
        ivl = pyfastaq.intervals.Interval

        tests = [
            (0, [], 0),
            (1, [], 1),
            (2, [], 2),
            (0, [ivl(3, 5)], 0),
            (1, [ivl(3, 5)], 1),
            (2, [ivl(3, 5)], 2),
            (3, [ivl(3, 5)], 6),
            (4, [ivl(3, 5)], 7),
            (5, [ivl(3, 5)], 8),
            (0, [ivl(3, 5), ivl(9,14)], 0),
            (1, [ivl(3, 5), ivl(9,14)], 1),
            (2, [ivl(3, 5), ivl(9,14)], 2),
            (3, [ivl(3, 5), ivl(9,14)], 6),
            (4, [ivl(3, 5), ivl(9,14)], 7),
            (5, [ivl(3, 5), ivl(9,14)], 8),
            (6, [ivl(3, 5), ivl(9,14)], 15),
            (7, [ivl(3, 5), ivl(9,14)], 16),
            (8, [ivl(3, 5), ivl(9,14)], 17),
        ]

        for position, insertions, expected in tests:
            got = aln_to_metadata.AlnToMetadata._unpadded_to_padded_nt_position(position, insertions)
            self.assertEqual(expected, got)


    def test_padded_to_unpadded_nt_position(self):
        '''test _padded_to_unpadded_nt_position'''
        ivl = pyfastaq.intervals.Interval

        tests = [
            (0, [], 0),
            (1, [], 1),
            (2, [], 2),
            (0, [ivl(3, 5)], 0),
            (1, [ivl(3, 5)], 1),
            (2, [ivl(3, 5)], 2),
            (3, [ivl(3, 5)], None),
            (4, [ivl(3, 5)], None),
            (5, [ivl(3, 5)], None),
            (6, [ivl(3, 5)], 3),
            (7, [ivl(3, 5)], 4),
            (8, [ivl(3, 5)], 5),
            (0, [ivl(3, 5), ivl(7,10)], 0),
            (1, [ivl(3, 5), ivl(7,10)], 1),
            (2, [ivl(3, 5), ivl(7,10)], 2),
            (3, [ivl(3, 5), ivl(7,10)], None),
            (4, [ivl(3, 5), ivl(7,10)], None),
            (5, [ivl(3, 5), ivl(7,10)], None),
            (6, [ivl(3, 5), ivl(7,10)], 3),
            (7, [ivl(3, 5), ivl(7,10)], None),
            (8, [ivl(3, 5), ivl(7,10)], None),
            (9, [ivl(3, 5), ivl(7,10)], None),
            (10, [ivl(3, 5), ivl(7,10)], None),
            (11, [ivl(3, 5), ivl(7,10)], 4),
            (12, [ivl(3, 5), ivl(7,10)], 5),
        ]

        for position, insertions, expected in tests:
            got = aln_to_metadata.AlnToMetadata._padded_to_unpadded_nt_position(position, insertions)
            self.assertEqual(expected, got)


    def test_variants_to_tsv_lines_coding(self):
        '''test _variants_to_tsv_lines coding sequences'''
        padded_seqs = {
            'seq1': pyfastaq.sequences.Fasta('seq1', 'ATG---GCTAATTAG'), # M-AN*
            'seq2': pyfastaq.sequences.Fasta('seq2', 'ATG---GCTAATTAG'), # MFAN*
            'seq3': pyfastaq.sequences.Fasta('seq3', 'ATGTTT---AATTAG'), # MF-N*
            'seq4': pyfastaq.sequences.Fasta('seq4', 'ATGTTTTGTAATTAG'), # MFCN*
            'seq5': pyfastaq.sequences.Fasta('seq5', 'ATGTTTGATAATTAG'), # MFDN*
        }

        unpadded_seqs = aln_to_metadata.AlnToMetadata._make_unpadded_seqs(padded_seqs)
        insertions = aln_to_metadata.AlnToMetadata._make_unpadded_insertion_coords(padded_seqs)

        variant1 = sequence_variant.Variant('p', 'A2D', 'id1')
        variant2 = sequence_variant.Variant('p', 'F2E', 'id2')
        variants = {
            'seq1': [(variant1, 'description 1')],
            'seq5': [(variant2, 'description 2')],
        }

        expected = [
            'seq1\t1\t0\tA2D\tid1\tdescription 1',
            'seq2\t1\t0\tA2D\tid1\tdescription 1',
            'seq4\t1\t0\tC3D\tid1\tdescription 1',
            'seq5\t1\t0\tA3D\tid1\tdescription 1',
            'seq5\t1\t0\tF2E\tid2\tdescription 2',
            'seq3\t1\t0\tF2E\tid2\tdescription 2',
            'seq4\t1\t0\tF2E\tid2\tdescription 2',
        ]

        got = aln_to_metadata.AlnToMetadata._variants_to_tsv_lines(variants, unpadded_seqs, padded_seqs, insertions, True, False)
        self.assertEqual(expected, got)


    def test_variants_to_tsv_lines_noncoding(self):
        '''test _variants_to_tsv_lines noncoding sequences'''
        padded_seqs = {
            'seq1': pyfastaq.sequences.Fasta('seq1', 'ATG---GCTAATTAG'),
            'seq2': pyfastaq.sequences.Fasta('seq2', 'ATG---GCTAATTAG'),
            'seq3': pyfastaq.sequences.Fasta('seq3', 'ATGTAT---AATTAG'),
            'seq4': pyfastaq.sequences.Fasta('seq4', 'ATGTGTTGTAATTAG'),
            'seq5': pyfastaq.sequences.Fasta('seq5', 'ATGTTTGATAATTAG'),
        }

        unpadded_seqs = aln_to_metadata.AlnToMetadata._make_unpadded_seqs(padded_seqs)
        insertions = aln_to_metadata.AlnToMetadata._make_unpadded_insertion_coords(padded_seqs)

        variant1 = sequence_variant.Variant('n', 'C5T', 'id1')
        variant2 = sequence_variant.Variant('n', 'A5T', 'id2')
        variants = {
            'seq1': [(variant1, 'description 1')],
            'seq5': [(variant2, 'description 2')],
        }

        expected = [
            'seq1\t0\t1\tC5T\tid1\tdescription 1',
            'seq2\t0\t1\tC5T\tid1\tdescription 1',
            'seq4\t0\t1\tG8T\tid1\tdescription 1',
            'seq5\t0\t1\tA8T\tid1\tdescription 1',
            'seq5\t0\t1\tA5T\tid2\tdescription 2',
            'seq3\t0\t1\tA5T\tid2\tdescription 2',
            'seq4\t0\t1\tG5T\tid2\tdescription 2',
        ]

        got = aln_to_metadata.AlnToMetadata._variants_to_tsv_lines(variants, unpadded_seqs, padded_seqs, insertions, False, True)
        self.assertEqual(expected, got)


    def test_make_cluster_file(self):
        '''test _make_cluster_file'''
        seqs = {
            'seq1': pyfastaq.sequences.Fasta('seq1', 'a'),
            'seq2': pyfastaq.sequences.Fasta('seq2', 'c'),
            'seq3': pyfastaq.sequences.Fasta('seq3', 'g'),
        }
        tmpfile = 'tmp.aln_to_meta_test_make_cluster_file.out'
        expected_file = os.path.join(data_dir, 'aln_to_metadata_make_cluster_file.out')
        aln_to_metadata.AlnToMetadata._make_cluster_file(seqs, tmpfile)
        self.assertTrue(filecmp.cmp(expected_file, tmpfile, shallow=False))
        os.unlink(tmpfile)


    def test_run_coding(self):
        '''test run coding sequences'''
        fa_in = os.path.join(data_dir, 'aln_to_metadata_run_coding.in.fa')
        fa_expected = os.path.join(data_dir, 'aln_to_metadata_run_coding.out.fa')
        tsv_in = os.path.join(data_dir, 'aln_to_metadata_run_coding.in.tsv')
        tsv_expected = os.path.join(data_dir, 'aln_to_metadata_run_coding.out.tsv')
        cluster_expected = os.path.join(data_dir, 'aln_to_metadata_run_coding.out.cluster')
        a_to_m = aln_to_metadata.AlnToMetadata(fa_in, tsv_in, True, False)
        outprefix = 'tmp.test.aln_to_metadata.run_coding'
        a_to_m.run(outprefix)
        self.assertTrue(filecmp.cmp(tsv_expected, outprefix + '.tsv', shallow=False))
        self.assertTrue(filecmp.cmp(fa_expected, outprefix + '.fa', shallow=False))
        self.assertTrue(filecmp.cmp(cluster_expected, outprefix + '.cluster', shallow=False))
        os.unlink(outprefix + '.tsv')
        os.unlink(outprefix + '.fa')
        os.unlink(outprefix + '.cluster')


    def test_run_noncoding(self):
        '''test run noncoding sequences'''
        fa_in = os.path.join(data_dir, 'aln_to_metadata_run_noncoding.in.fa')
        fa_expected = os.path.join(data_dir, 'aln_to_metadata_run_noncoding.out.fa')
        tsv_in = os.path.join(data_dir, 'aln_to_metadata_run_noncoding.in.tsv')
        tsv_expected = os.path.join(data_dir, 'aln_to_metadata_run_noncoding.out.tsv')
        cluster_expected = os.path.join(data_dir, 'aln_to_metadata_run_noncoding.out.cluster')
        a_to_m = aln_to_metadata.AlnToMetadata(fa_in, tsv_in, False, True)
        outprefix = 'tmp.test.aln_to_metadata.run_noncoding'
        a_to_m.run(outprefix)
        self.assertTrue(filecmp.cmp(tsv_expected, outprefix + '.tsv', shallow=False))
        self.assertTrue(filecmp.cmp(fa_expected, outprefix + '.fa', shallow=False))
        self.assertTrue(filecmp.cmp(cluster_expected, outprefix + '.cluster', shallow=False))
        os.unlink(outprefix + '.tsv')
        os.unlink(outprefix + '.fa')
        os.unlink(outprefix + '.cluster')

