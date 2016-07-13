import unittest
import filecmp
import os
import pyfastaq
from ariba import reference_data, sequence_metadata

modules_dir = os.path.dirname(os.path.abspath(reference_data.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


class TestReferenceData(unittest.TestCase):
    def test_init_fails(self):
        '''Test __init__ fails when it should'''
        empty_fasta = os.path.join(data_dir, 'reference_data_init_fails.empty.fa')
        empty_tsv = os.path.join(data_dir, 'reference_data_init_fails.empty.tsv')
        fasta = os.path.join(data_dir, 'reference_data_init_fails.in.fa')

        with self.assertRaises(reference_data.Error):
            ref_data = reference_data.ReferenceData([empty_fasta], [empty_tsv])
            ref_data = reference_data.ReferenceData([fasta], [empty_tsv])


    def test_init_ok(self):
        '''Test init with good input'''
        fasta_in = os.path.join(data_dir, 'reference_data_init_ok.in.fa')
        tsv_in = os.path.join(data_dir, 'reference_data_init_ok.in.tsv')
        meta1 = sequence_metadata.SequenceMetadata('gene1\t1\t0\tR2S\t.\tconfers killer rabbit resistance')
        meta2 = sequence_metadata.SequenceMetadata("gene2\t1\t0\tI42L\t.\tremoves tardigrade's space-living capability")

        expected_metadata = {
            'gene1': {
                'seq_type': 'p',
                'variant_only': False,
                'n': {},
                'p': {1: {meta1}},
                '.': set(),
            },
            'gene2': {
                'seq_type': 'p',
                'variant_only': False,
                'n': {},
                'p': {41: {meta2}},
                '.': set(),
            }
        }
        ref_data = reference_data.ReferenceData([fasta_in], [tsv_in])
        self.assertEqual(expected_metadata, ref_data.metadata)

        expected_seqs_dict = {
            'gene1': pyfastaq.sequences.Fasta('gene1', 'CATCGTCGTCTATCGTCGTCCTAG'),
            'gene2': pyfastaq.sequences.Fasta('gene2', 'AAAAACCCCGGGGTTTT')
        }

        self.assertEqual(expected_seqs_dict, ref_data.sequences)


    def test_load_metadata_tsv(self):
        '''Test _load_metadata_tsv'''
        meta1 = sequence_metadata.SequenceMetadata('gene1\t0\t0\tA42G\t.\tfree text')
        meta2 = sequence_metadata.SequenceMetadata('gene1\t0\t0\tG13T\t.\tconfers killer rabbit resistance')
        meta3 = sequence_metadata.SequenceMetadata("gene2\t1\t1\tI42L\t.\tremoves tardigrade's space-living capability")
        expected = {
            'gene1': {
                'seq_type': 'n',
                'variant_only': False,
                'n': {12: {meta2}, 41: {meta1}},
                'p': {},
                '.': set(),
            },
            'gene2': {
                'seq_type': 'p',
                'variant_only': True,
                'n': {},
                'p': {41: {meta3}},
                '.': set(),
            }
        }

        got = {}
        tsv_file = os.path.join(data_dir, 'reference_data_load_metadata_tsv.tsv')
        reference_data.ReferenceData._load_metadata_tsv(tsv_file, got)
        self.assertEqual(expected, got)


    def test_load_all_metadata_tsvs(self):
       '''Test _load_all_metadata_tsvs'''
       input_files = [os.path.join(data_dir, 'reference_data_load_all_metadata_tsvs.' + x + '.tsv') for x in ['1', '2']]
       meta1 = sequence_metadata.SequenceMetadata('gene1\t0\t0\tA42G\t.\tfree text')
       meta2 = sequence_metadata.SequenceMetadata('gene1\t0\t0\tG13T\t.\tconfers killer rabbit resistance')
       meta3 = sequence_metadata.SequenceMetadata("gene2\t1\t0\tI42L\t.\tremoves tardigrade's space-living capability")
       expected = {
           'gene1': {
                'seq_type': 'n',
                'variant_only': False,
               'n': {12: {meta2}, 41: {meta1}},
               'p': {},
               '.': set(),
           },
           'gene2': {
                'seq_type': 'p',
                'variant_only': False,
               'n': {},
               'p': {41: {meta3}},
               '.': set(),
           }
       }

       got = reference_data.ReferenceData._load_all_metadata_tsvs(input_files)
       self.assertEqual(expected, got)


    def test_load_fasta_file(self):
        '''Test _load_fasta_file'''
        got = {}
        expected = {'seq1': pyfastaq.sequences.Fasta('seq1', 'ACGT')}
        filename = os.path.join(data_dir, 'reference_data_load_fasta_file.fa')
        reference_data.ReferenceData._load_fasta_file(filename, got)
        self.assertEqual(expected, got)


    def test_load_all_fasta_files(self):
        '''Test _load_all_fasta_files'''
        filenames = [os.path.join(data_dir, 'reference_data_load_all_fasta_files.in.' + x) for x in ['1', '2']]
        expected = {
            'seq1': pyfastaq.sequences.Fasta('seq1', 'ACGT'),
            'seq2': pyfastaq.sequences.Fasta('seq2', 'TTTT')
        }
        got = reference_data.ReferenceData._load_all_fasta_files(filenames)
        self.assertEqual(expected, got)


    def test_load_input_files_and_check_seq_names_ok(self):
        '''Test _load_input_files_and_check_seq_names with good input'''
        fasta_files = [os.path.join(data_dir, 'reference_data_load_input_files_and_check_seq_names.good.in.fa.' + x) for x in ['1', '2']]
        metadata_files = [os.path.join(data_dir, 'reference_data_load_input_files_and_check_seq_names.good.in.csv.' + x) for x in ['1', '2']]
        expected_seqs = {
             'seq1': pyfastaq.sequences.Fasta('seq1', 'ACGT'),
             'seq2': pyfastaq.sequences.Fasta('seq2', 'TTTT')
        }
        meta1 = sequence_metadata.SequenceMetadata('seq1\t0\t0\tA1G\t.\tfree text')
        meta2 = sequence_metadata.SequenceMetadata("seq2\t0\t0\t.\t.\tspam eggs")
        expected_meta = {
            'seq1': {
               'seq_type': 'n',
               'variant_only': False,
               'n': {0: {meta1}},
               'p': {},
               '.': set(),
            },
            'seq2': {
               'seq_type': 'n',
               'variant_only': False,
               'n': {},
               'p': {},
               '.': {meta2},
            }
        }
        got_seqs, got_meta = reference_data.ReferenceData._load_input_files_and_check_seq_names(fasta_files, metadata_files)
        self.assertEqual(expected_seqs, got_seqs)
        self.assertEqual(expected_meta, got_meta)


    def test_load_input_files_and_check_seq_names_bad(self):
        '''Test _load_input_files_and_check_seq_names with bad input'''
        fasta_files = [os.path.join(data_dir, 'reference_data_load_input_files_and_check_seq_names.bad.in.fa.' + x) for x in ['1', '2']]
        metadata_files = [os.path.join(data_dir, 'reference_data_load_input_files_and_check_seq_names.bad.in.csv.' + x) for x in ['1', '2']]
        with self.assertRaises(reference_data.Error):
            reference_data.ReferenceData._load_input_files_and_check_seq_names(fasta_files, metadata_files)


    def test_find_gene_in_seqs(self):
        '''Test _find_gene_in_seqs'''
        seqs_dict = {
            'dict1': {'name1': 'seq1', 'name2': 'seq2'},
            'dict2': {'name3': 'seq3'}
        }
        self.assertEqual(None, reference_data.ReferenceData._find_gene_in_seqs('name42', seqs_dict))
        self.assertEqual('dict1', reference_data.ReferenceData._find_gene_in_seqs('name1', seqs_dict))
        self.assertEqual('dict1', reference_data.ReferenceData._find_gene_in_seqs('name2', seqs_dict))
        self.assertEqual('dict2', reference_data.ReferenceData._find_gene_in_seqs('name3', seqs_dict))


    def test_write_metadata_tsv(self):
        '''Test _write_metadata_tsv'''
        metadata_tsv_in = os.path.join(data_dir, 'reference_data_write_metadata_tsv.tsv')
        metadata_tsv_expected = os.path.join(data_dir, 'reference_data_write_metadata_tsv.expected.tsv')
        tmp_tsv = 'tmp.test_write_metadata_tsv.out.tsv'
        metadata = reference_data.ReferenceData._load_all_metadata_tsvs([metadata_tsv_in])
        reference_data.ReferenceData._write_metadata_tsv(metadata, tmp_tsv)
        self.assertTrue(filecmp.cmp(metadata_tsv_expected, tmp_tsv, shallow=False))
        os.unlink(tmp_tsv)


    def test_write_sequences_to_files(self):
        '''Test _write_sequences_to_files'''
        sequences = {
            'seq1': pyfastaq.sequences.Fasta('seq1', 'ACGT'),
            'seq2': pyfastaq.sequences.Fasta('seq2', 'ACGTA'),
            'seq3': pyfastaq.sequences.Fasta('seq3', 'ACGTAC'),
            'seq4': pyfastaq.sequences.Fasta('seq4', 'ACGTAAA'),
            'seq5': pyfastaq.sequences.Fasta('seq5', 'ACGTCCC'),
        }
        metadata = {
            'seq1': {'seq_type': 'n', 'variant_only': False},
            'seq2': {'seq_type': 'n', 'variant_only': True},
            'seq3': {'seq_type': 'p', 'variant_only': False},
            'seq4': {'seq_type': 'p', 'variant_only': True},
            'seq5': {'seq_type': 'n', 'variant_only': False},
        }
        tmp_prefix = 'tmp.test_write_sequences_to_files'
        reference_data.ReferenceData._write_sequences_to_files(sequences, metadata, tmp_prefix)
        expected_prefix = os.path.join(data_dir, 'reference_data_write_sequences_to_files')
        for suffix in ['gene.fa', 'gene.varonly.fa', 'noncoding.fa', 'noncoding.varonly.fa', 'all.fa']:
            expected = expected_prefix + '.' + suffix
            got = tmp_prefix + '.' + suffix
            self.assertTrue(filecmp.cmp(expected, got, shallow=False))
            os.unlink(got)


    def test_filter_bad_variant_data(self):
        '''Test _filter_bad_variant_data'''
        fasta_in = os.path.join(data_dir, 'reference_data_filter_bad_data.in.fa')
        metadata_tsv = os.path.join(data_dir, 'reference_data_filter_bad_data_metadata.in.tsv')
        sequences, metadata = reference_data.ReferenceData._load_input_files_and_check_seq_names([fasta_in], [metadata_tsv])
        tmp_prefix = 'tmp.test_filter_bad_variant_data'
        reference_data.ReferenceData._filter_bad_variant_data(sequences, metadata, tmp_prefix, set())
        expected_prefix = os.path.join(data_dir, 'reference_data_filter_bad_data.expected')

        for suffix in ['gene.fa', 'gene.varonly.fa', 'noncoding.fa', 'noncoding.varonly.fa', 'all.fa', 'log', 'metadata.tsv']:
            expected = expected_prefix + '.' + suffix
            got = tmp_prefix + '.' + suffix
            self.assertTrue(filecmp.cmp(expected, got, shallow=False))
            os.unlink(got)


    def test_try_to_get_gene_seq(self):
        '''Test _try_to_get_gene_seq'''
        tests = [
            (pyfastaq.sequences.Fasta('x', 'ACGTG'), None, 'Remove: too short. Length: 5'),
            (pyfastaq.sequences.Fasta('x', 'A' * 100), None, 'Remove: too long. Length: 100'),
            (pyfastaq.sequences.Fasta('x', 'GAGGAGCCG'), None, 'Does not look like a gene (tried both strands and all reading frames) GAGGAGCCG'),
            (pyfastaq.sequences.Fasta('x', 'ATGTAACCT'), None, 'Does not look like a gene (tried both strands and all reading frames) ATGTAACCT'),
            (pyfastaq.sequences.Fasta('x', 'ATGCCTTAA'), pyfastaq.sequences.Fasta('x', 'ATGCCTTAA'), 'Made x into gene. strand=+, frame=0')
        ]

        for seq, got_seq, message in tests:
            self.assertEqual((got_seq, message), reference_data.ReferenceData._try_to_get_gene_seq(seq, 6, 99))


    def test_remove_bad_genes(self):
        '''Test _remove_bad_genes'''
        test_seq_dict = {}
        fasta_file = os.path.join(data_dir, 'reference_data_remove_bad_genes.in.fa')
        metadata_file = os.path.join(data_dir, 'reference_data_remove_bad_genes.in.tsv')
        metadata = reference_data.ReferenceData._load_all_metadata_tsvs([metadata_file])
        pyfastaq.tasks.file_to_dict(fasta_file, test_seq_dict)
        tmp_log = 'tmp.test_remove_bad_genes.log'
        expected_removed = {'g1', 'g2', 'g3', 'g4'}
        got_removed = reference_data.ReferenceData._remove_bad_genes(test_seq_dict, metadata, tmp_log, min_gene_length=6, max_gene_length=99)
        self.assertEqual(expected_removed, got_removed)
        expected_dict = {
            'g5': pyfastaq.sequences.Fasta('g5', 'ATGCCTTAA'),
            'noncoding1': pyfastaq.sequences.Fasta('noncoding1', 'AAAAAAAAAAAAAAAAAAAAAAA')
        }
        self.assertEqual(expected_dict, test_seq_dict)
        expected_log = os.path.join(data_dir, 'reference_data_test_remove_bad_genes.log')
        self.assertTrue(filecmp.cmp(expected_log, tmp_log, shallow=False))
        os.unlink(tmp_log)


    def test_new_seq_name(self):
        '''Test _new_seq_name'''
        tests = [
            ('name', 'name'),
            ('name ', 'name'),
            ('name xyz', 'name'),
            ('name_a', 'name_a'),
            ('name.a', 'name.a'),
            ('name-a', 'name-a'),
            ('name spam eggs foo', 'name'),
            ('name!', 'name_'),
            ('name:foo', 'name_foo'),
            ('name:!@foo', 'name___foo'),
        ]

        for name, expected in tests:
            self.assertEqual(expected, reference_data.ReferenceData._new_seq_name(name))


    def test_seq_names_to_rename_dict(self):
        '''Test _seq_names_to_rename_dict'''
        names = {'foo', 'foo abc', 'foo xyz', 'bar!', 'bar:', 'spam abc', 'eggs'}
        got = reference_data.ReferenceData._seq_names_to_rename_dict(names)
        expected = {
            'foo abc': 'foo_1',
            'foo xyz': 'foo_2',
            'bar!': 'bar_',
            'bar:': 'bar__1',
            'spam abc': 'spam'
        }
        self.assertEqual(expected, got)


    def test_rename_names_in_seq_dict(self):
        '''Test _rename_names_in_seq_dict'''
        original_seqs = {
            'pa abc': pyfastaq.sequences.Fasta('pa abc', 'AAAA'),
            'pa 1': pyfastaq.sequences.Fasta('pa 1', 'CCC'),
            'vo:': pyfastaq.sequences.Fasta('vo:', 'GGG'),
            'nonc': pyfastaq.sequences.Fasta('nonc', 'TTT'),
        }
        rename_dict = {
            'pa abc': 'pa',
            'pa 1': 'pa_1',
            'vo:': 'vo_',
        }
        expected = {
            'pa': pyfastaq.sequences.Fasta('pa', 'AAAA'),
            'pa_1': pyfastaq.sequences.Fasta('pa_1', 'CCC'),
            'vo_': pyfastaq.sequences.Fasta('vo_', 'GGG'),
            'nonc': pyfastaq.sequences.Fasta('nonc', 'TTT'),
        }

        got = reference_data.ReferenceData._rename_names_in_seq_dict(original_seqs, rename_dict)
        self.assertEqual(expected, got)


    def test_rename_metadata_set(self):
        '''Test _rename_metadata_set'''
        metaset = {
            sequence_metadata.SequenceMetadata('foo 1\t1\t0\t.\t.\tdescription'),
            sequence_metadata.SequenceMetadata('foo 1\t1\t0\tI42L\t.\tspam eggs')
        }

        expected = {
            sequence_metadata.SequenceMetadata('new_name\t1\t0\t.\t.\tdescription'),
            sequence_metadata.SequenceMetadata('new_name\t1\t0\tI42L\t.\tspam eggs')
        }
        got = reference_data.ReferenceData._rename_metadata_set(metaset, 'new_name')
        self.assertEqual(expected, got)


    def test_rename_names_in_metadata(self):
        '''Test _rename_names_in_metadata'''
        meta1 = sequence_metadata.SequenceMetadata('gene1\t0\t0\tA42G\t.\tfree text')
        meta2 = sequence_metadata.SequenceMetadata('gene1\t0\t0\tA42T\t.\tfree text2')
        meta3 = sequence_metadata.SequenceMetadata('gene1\t0\t0\t.\t.\tfree text3')
        meta4 = sequence_metadata.SequenceMetadata('gene1\t0\t0\tG13T\t.\tconfers killer rabbit resistance')
        meta5 = sequence_metadata.SequenceMetadata("gene2\t1\t0\tI42L\t.\tremoves tardigrade's space-living capability")
        meta1rename = sequence_metadata.SequenceMetadata('new_gene1\t0\t0\tA42G\t.\tfree text')
        meta2rename = sequence_metadata.SequenceMetadata('new_gene1\t0\t0\tA42T\t.\tfree text2')
        meta3rename = sequence_metadata.SequenceMetadata('new_gene1\t0\t0\t.\t.\tfree text3')
        meta4rename = sequence_metadata.SequenceMetadata('new_gene1\t0\t0\tG13T\t.\tconfers killer rabbit resistance')

        metadata = {
            'gene1': {
                'n': {12: {meta4}, 41: {meta1, meta2}},
                'p': {},
                '.': {meta3},
            },
            'gene2': {
                'n': {},
                'p': {41: {meta5}},
                '.': set(),
            }
        }

        expected = {
            'new_gene1': {
                'n': {12: {meta4rename}, 41: {meta1rename, meta2rename}},
                'p': {},
                '.': {meta3rename},
            },
            'gene2': {
                'n': {},
                'p': {41: {meta5}},
                '.': set(),
            }
        }

        rename_dict = {'gene1': 'new_gene1'}
        got = reference_data.ReferenceData._rename_names_in_metadata(metadata, rename_dict)
        self.assertEqual(expected, got)


    def test_rename_sequences(self):
        '''Test rename_sequences'''
        fasta_in = os.path.join(data_dir, 'reference_data_rename_sequences.fa')
        tsv_in = os.path.join(data_dir, 'reference_data_rename_sequences_metadata.tsv')
        refdata = reference_data.ReferenceData([fasta_in], [tsv_in])
        tmp_out = 'tmp.test_rename_sequences.out'
        refdata.rename_sequences(tmp_out)
        expected_file = os.path.join(data_dir, 'reference_data_test_rename_sequences.out')
        self.assertTrue(filecmp.cmp(expected_file, tmp_out, shallow=False))
        os.unlink(tmp_out)

        meta1 = sequence_metadata.SequenceMetadata('noncoding1\t0\t0\t.\t.\toriginal name "noncoding1"')
        meta2 = sequence_metadata.SequenceMetadata('noncoding1_1\t0\t0\t.\t.\toriginal name "noncoding1 blah"')
        meta3 = sequence_metadata.SequenceMetadata('pres_abs1_2\t0\t0\t.\t.\toriginal name "pres_abs1 foo bar spam eggs"')
        meta4 = sequence_metadata.SequenceMetadata('pres_abs1_1\t0\t0\t.\t.\toriginal name "pres_abs1 blah"')
        meta5 = sequence_metadata.SequenceMetadata('pres_abs1\t0\t0\t.\t.\toriginal name "pres\'abs1"')
        meta6 = sequence_metadata.SequenceMetadata('pres_abs2\t0\t0\t.\t.\toriginal name "pres_abs2"')
        meta7 = sequence_metadata.SequenceMetadata('pres_abs3\t0\t0\t.\t.\toriginal name "pres!abs3"')
        meta8 = sequence_metadata.SequenceMetadata('var_only1_2\t0\t0\t.\t.\toriginal name "var_only1 hello"')
        meta9 = sequence_metadata.SequenceMetadata('var_only1\t0\t0\t.\t.\toriginal name "var:only1 boo"')
        meta10 = sequence_metadata.SequenceMetadata('var_only1_1\t0\t0\t.\t.\toriginal name "var_only1"')
        meta11 = sequence_metadata.SequenceMetadata('var_only2\t0\t0\t.\t.\toriginal name "var_only2"')

        expected_meta = {
            'noncoding1': {'seq_type': 'n', 'variant_only': False, 'n': {}, 'p': {}, '.': {meta1}},
            'noncoding1_1': {'seq_type': 'n', 'variant_only': False, 'n': {}, 'p': {}, '.': {meta2}},
            'pres_abs1_2': {'seq_type': 'n', 'variant_only': False, 'n': {}, 'p': {}, '.': {meta3}},
            'pres_abs1_1': {'seq_type': 'n', 'variant_only': False, 'n': {}, 'p': {}, '.': {meta4}},
            'pres_abs1': {'seq_type': 'n', 'variant_only': False, 'n': {}, 'p': {}, '.': {meta5}},
            'pres_abs2': {'seq_type': 'n', 'variant_only': False, 'n': {}, 'p': {}, '.': {meta6}},
            'pres_abs3': {'seq_type': 'n', 'variant_only': False, 'n': {}, 'p': {}, '.': {meta7}},
            'var_only1_2': {'seq_type': 'n', 'variant_only': False, 'n': {}, 'p': {}, '.': {meta8}},
            'var_only1': {'seq_type': 'n', 'variant_only': False, 'n': {}, 'p': {}, '.': {meta9}},
            'var_only1_1': {'seq_type': 'n', 'variant_only': False, 'n': {}, 'p': {}, '.': {meta10}},
            'var_only2': {'seq_type': 'n', 'variant_only': False, 'n': {}, 'p': {}, '.': {meta11}},
        }

        self.assertEqual(expected_meta, refdata.metadata)

        expected_seqs_dict = {
            'noncoding1': pyfastaq.sequences.Fasta('noncoding1', 'AAAA'),
            'noncoding1_1': pyfastaq.sequences.Fasta('noncoding1_1', 'CCCC'),
            'pres_abs1_2': pyfastaq.sequences.Fasta('pres_abs1_2', 'ACGT'),
            'pres_abs1_1': pyfastaq.sequences.Fasta('pres_abs1_1', 'AAAA'),
            'pres_abs1': pyfastaq.sequences.Fasta('pres_abs1', 'CCCC'),
            'pres_abs2': pyfastaq.sequences.Fasta('pres_abs2', 'TTTT'),
            'pres_abs3': pyfastaq.sequences.Fasta('pres_abs3', 'GGGG'),
            'var_only1_2': pyfastaq.sequences.Fasta('var_only1_2', 'AAAA'),
            'var_only1': pyfastaq.sequences.Fasta('var_only1', 'CCCC'),
            'var_only1_1': pyfastaq.sequences.Fasta('var_only1_1', 'GGGG'),
            'var_only2': pyfastaq.sequences.Fasta('var_only2', 'TTTT'),
        }

        self.assertEqual(expected_seqs_dict, refdata.sequences)


    def test_sequence_type(self):
        '''Test sequence_type'''
        presence_absence_fa = os.path.join(data_dir, 'reference_data_sequence_type.presence_absence.fa')
        variants_only_fa = os.path.join(data_dir, 'reference_data_sequence_type.variants_only.fa')
        noncoding_fa = os.path.join(data_dir, 'reference_data_sequence_type.noncoding.fa')
        refdata = reference_data.ReferenceData(
            presence_absence_fa=presence_absence_fa,
            variants_only_fa=variants_only_fa,
            non_coding_fa=noncoding_fa
        )

        tests = [
            ('pa', 'presence_absence'),
            ('var_only', 'variants_only'),
            ('noncoding', 'non_coding'),
            ('not_there', None)
        ]

        for name, expected in tests:
            self.assertEqual(expected, refdata.sequence_type(name))


    def test_sequence(self):
        '''Test sequence'''
        fasta_in = os.path.join(data_dir, 'reference_data_sequence.in.fa')
        tsv_in = os.path.join(data_dir, 'reference_data_sequence.in.tsv')
        expected = pyfastaq.sequences.Fasta('seq1', 'ATGTTTTAA')
        refdata = reference_data.ReferenceData([fasta_in], [tsv_in])
        self.assertEqual(expected, refdata.sequence('seq1'))


    def test_all_non_wild_type_variants(self):
        '''Test all_non_wild_type_variants'''
        tsv_file = os.path.join(data_dir, 'reference_data_test_all_non_wild_type_variants.tsv')
        fasta_in = os.path.join(data_dir, 'reference_data_test_all_non_wild_type_variants.ref.fa')
        refdata = reference_data.ReferenceData([fasta_in], [tsv_file])

        v1 = sequence_metadata.SequenceMetadata('var_only_gene\t1\t1\tP3Q\t.\tref has wild type P')
        v2 = sequence_metadata.SequenceMetadata('var_only_gene\t1\t1\tG4I\t.\tref has wild type F')
        v3 = sequence_metadata.SequenceMetadata('var_only_gene\t1\t1\tI5V\t.\tref has variant V instead of I')
        v4 = sequence_metadata.SequenceMetadata('var_only_gene\t1\t1\tF6I\t.\tref has wild type F')
        p1 = sequence_metadata.SequenceMetadata('presence_absence_gene\t1\t0\tN2I\t.\tref has wild type N')
        p2 = sequence_metadata.SequenceMetadata('presence_absence_gene\t1\t0\tA4G\t.\tref has variant G instead of A')
        n1 = sequence_metadata.SequenceMetadata('non_coding\t0\t0\tA2C\t.\tref has wild type A')
        n2 = sequence_metadata.SequenceMetadata('non_coding\t0\t0\tC4T\t.\tref has variant T instead of C')

        var_only_expected = {
             'n': {},
             'p': {2: {v1}, 3: {v2}, 4: {v3}, 5: {v4}}
        }

        pres_abs_expected = {
            'n': {},
            'p': {1: {p1}, 3: {p2}},
        }

        non_coding_expected = {
            'n': {1: {n1}, 3: {n2}},
            'p': {}
        }

        self.assertEqual(var_only_expected, refdata.all_non_wild_type_variants('var_only_gene'))
        self.assertEqual(pres_abs_expected, refdata.all_non_wild_type_variants('presence_absence_gene'))
        self.assertEqual(non_coding_expected, refdata.all_non_wild_type_variants('non_coding'))
        self.assertEqual({'n': {}, 'p': {}}, refdata.all_non_wild_type_variants('not_a_known_sequence'))


    def test_write_cluster_allocation_file(self):
        '''Test write_cluster_allocation_file'''
        clusters = {
            '0': {'cluster0.1', 'cluster0.2'},
            '1': {'cluster1.1', 'cluster1.2'},
            '11': {'cluster11.1', 'cluster11.2'},
            '2': {'cluster2.1'}
        }

        tmpfile = 'tmp.test_write_cluster_allocation_file.out'
        reference_data.ReferenceData.write_cluster_allocation_file(clusters, tmpfile)
        expected_file = os.path.join(data_dir, 'reference_data_test_write_cluster_allocation_file.expected')
        self.assertTrue(filecmp.cmp(expected_file, tmpfile, shallow=False))
        os.unlink(tmpfile)


    def test_cluster_with_cdhit(self):
        '''Test cluster_with_cd_hit'''
        fasta_in = os.path.join(data_dir, 'reference_data_test_cluster_with_cdhit.in.fa')
        tsv_in = os.path.join(data_dir, 'reference_data_test_cluster_with_cdhit.in.tsv')
        refdata = reference_data.ReferenceData([fasta_in], [tsv_in])
        outprefix = 'tmp.test_cluster_with_cdhit'

        expected_clusters = {
            '0': {'noncoding1'},
            '1': {'presence_absence1', 'presence_absence2'},
            '2': {'presence_absence3', 'presence_absence4'},
        }

        got_clusters = refdata.cluster_with_cdhit(outprefix)
        self.assertEqual(expected_clusters, got_clusters)

        expected_clusters_file = os.path.join(data_dir, 'reference_data_test_cluster_with_cdhit.expected.clusters.tsv')
        got_clusters_file = outprefix + '.clusters.tsv'
        self.assertTrue(filecmp.cmp(expected_clusters_file, got_clusters_file, shallow=False))

        os.unlink(got_clusters_file)
        os.unlink(outprefix + '.all.fa')
        os.unlink(outprefix + '.gene.fa')
        os.unlink(outprefix + '.gene.varonly.fa')
        os.unlink(outprefix + '.noncoding.fa')
        os.unlink(outprefix + '.noncoding.varonly.fa')


    def test_cluster_with_cdhit_clusters_in_file(self):
        '''Test cluster_with_cd_hit clusters from file'''
        fasta_in = os.path.join(data_dir, 'reference_data_test_cluster_with_cdhit_clusters_in_file.in.fa')
        meta_tsv_in = os.path.join(data_dir, 'reference_data_test_cluster_with_cdhit_clusters_in_file.in.meta.tsv')
        cluster_tsv_in = os.path.join(data_dir, 'reference_data_test_cluster_with_cdhit_clusters_in_file.in.clusters.tsv')
        refdata = reference_data.ReferenceData([fasta_in], [meta_tsv_in])
        outprefix = 'tmp.test_cluster_with_cdhit_clusters_in_file'

        expected_clusters = {
            '0': {'presence_absence1', 'presence_absence3', 'presence_absence4'},
            '1': {'presence_absence2'},
            '2': {'noncoding1'},
            '3': {'noncoding2'},
        }

        got_clusters = refdata.cluster_with_cdhit(outprefix, clusters_file=cluster_tsv_in)
        print(got_clusters)
        self.assertEqual(expected_clusters, got_clusters)

        expected_clusters_file = os.path.join(data_dir, 'reference_data_test_cluster_with_cdhit_clusters_in_file.expected.clusters.tsv')
        got_clusters_file = outprefix + '.clusters.tsv'
        self.assertTrue(filecmp.cmp(expected_clusters_file, got_clusters_file, shallow=False))

        os.unlink(got_clusters_file)
        os.unlink(outprefix + '.all.fa')
        os.unlink(outprefix + '.gene.fa')
        os.unlink(outprefix + '.gene.varonly.fa')
        os.unlink(outprefix + '.noncoding.fa')
        os.unlink(outprefix + '.noncoding.varonly.fa')


    def test_cluster_with_cdhit_nocluster(self):
        '''Test cluster_with_cd_hit do not run cdhit'''
        fasta_in = os.path.join(data_dir, 'reference_data_test_cluster_with_cdhit_nocluster.in.fa')
        tsv_in = os.path.join(data_dir, 'reference_data_test_cluster_with_cdhit_nocluster.in.tsv')
        refdata = reference_data.ReferenceData([fasta_in], [tsv_in])
        outprefix = 'tmp.test_cluster_with_cdhit_nocluster'

        expected_clusters = {
            '0': {'noncoding1'},
            '1': {'noncoding2'},
            '2': {'presence_absence1'},
            '3': {'presence_absence2'},
            '4': {'presence_absence3'},
            '5': {'presence_absence4'},
        }

        got_clusters = refdata.cluster_with_cdhit(outprefix, nocluster=True)
        self.assertEqual(expected_clusters, got_clusters)

        expected_clusters_file = os.path.join(data_dir, 'reference_data_test_cluster_with_cdhit_nocluster.expected.clusters.tsv')
        got_clusters_file = outprefix + '.clusters.tsv'
        self.assertTrue(filecmp.cmp(expected_clusters_file, got_clusters_file, shallow=False))

        os.unlink(got_clusters_file)
        os.unlink(outprefix + '.all.fa')
        os.unlink(outprefix + '.gene.fa')
        os.unlink(outprefix + '.gene.varonly.fa')
        os.unlink(outprefix + '.noncoding.fa')
        os.unlink(outprefix + '.noncoding.varonly.fa')



    def test_write_seqs_to_fasta(self):
        '''Test write_seqs_to_fasta'''
        refdata = reference_data.ReferenceData(presence_absence_fa=os.path.join(data_dir, 'reference_data_test_write_seqs_to_fasta.in.fa'))
        expected_outfile = os.path.join(data_dir, 'reference_data_test_write_seqs_to_fasta.expected.fa')
        tmpfile = 'tmp.test.reference_data.write_seqs_to_fasta.out.fa'
        refdata.write_seqs_to_fasta(tmpfile, {'seq1', 'seq4', 'seq5'})
        self.assertTrue(filecmp.cmp(expected_outfile, tmpfile, shallow=False))
        os.unlink(tmpfile)

