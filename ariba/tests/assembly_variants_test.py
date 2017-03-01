import unittest
import os
import pymummer
import pyfastaq
from ariba import assembly_variants, reference_data, sequence_metadata

modules_dir = os.path.dirname(os.path.abspath(assembly_variants.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


class TestAssemblyVariants(unittest.TestCase):
    def test_get_codon_start(self):
        '''test _get_codon_start'''
        tests = [
            (0, 5, 3),
            (0, 0, 0),
            (0, 1, 0),
            (0, 2, 0),
            (1, 3, 1),
            (2, 3, 2),
            (3, 3, 3),
            (3, 6, 6),
            (3, 7, 6),
            (3, 8, 6),
        ]
        for start, position, expected in tests:
            self.assertEqual(expected, assembly_variants.AssemblyVariants._get_codon_start(start, position))


    def test_get_mummer_variants_no_variants(self):
        '''test _get_mummer_variants when no variants'''
        snp_file = os.path.join(data_dir, 'assembly_variants_test_get_mummer_variants.none.snps')
        got = assembly_variants.AssemblyVariants._get_mummer_variants(snp_file)
        self.assertEqual({}, got)


    def test_get_mummer_variants_has_variants(self):
        '''test _get_mummer_variants when there are variants'''
        snp_file = os.path.join(data_dir, 'assembly_variants_test_get_mummer_variants.snp.snps')
        v1 = pymummer.variant.Variant(pymummer.snp.Snp('42\tA\tG\t42\t42\t42\t500\t500\t1\t1\tgene\tcontig1'))
        v2 = pymummer.variant.Variant(pymummer.snp.Snp('42\tA\tG\t42\t42\t42\t500\t500\t1\t1\tgene\tcontig2'))
        v3 = pymummer.variant.Variant(pymummer.snp.Snp('40\tT\tC\t40\t42\t42\t500\t500\t1\t1\tgene\tcontig1'))
        v4 = pymummer.variant.Variant(pymummer.snp.Snp('2\tC\tG\t2\t42\t42\t500\t500\t1\t1\tgene\tcontig1'))
        expected = {
            'contig1': [[v4], [v3, v1]],
            'contig2': [[v2]]
        }
        got = assembly_variants.AssemblyVariants._get_mummer_variants(snp_file)
        self.assertEqual(expected, got)


    def test_get_variant_effect(self):
        '''test _get_variant_effect'''
        ref_seq = pyfastaq.sequences.Fasta('gene', 'GATCGCGAAGCGATGACCCATGAAGCGACCGAACGCTGA')
        v1 = pymummer.variant.Variant(pymummer.snp.Snp('6\tC\tT\t6\tx\tx\t39\t39\tx\tx\tgene\tcontig'))
        v1 = pymummer.variant.Variant(pymummer.snp.Snp('6\tC\tT\t6\tx\tx\t39\t39\tx\tx\tgene\tcontig'))
        v2 = pymummer.variant.Variant(pymummer.snp.Snp('4\tC\tA\t6\tx\tx\t39\t39\tx\tx\tgene\tcontig'))
        v3 = pymummer.variant.Variant(pymummer.snp.Snp('4\tC\tT\t4\tx\tx\t39\t39\tx\tx\tgene\tcontig'))
        v4 = pymummer.variant.Variant(pymummer.snp.Snp('6\tC\tA\t6\tx\tx\t39\t39\tx\tx\tgene\tcontig'))
        v5 = pymummer.variant.Variant(pymummer.snp.Snp('4\tC\t.\t4\tx\tx\t39\t39\tx\tx\tgene\tcontig'))
        v6 = pymummer.variant.Variant(pymummer.snp.Snp('4\t.\tA\t4\tx\tx\t39\t39\tx\tx\tgene\tcontig'))
        v7 = pymummer.variant.Variant(pymummer.snp.Snp('4\t.\tG\t4\tx\tx\t39\t39\tx\tx\tgene\tcontig'))
        v7.qry_base = 'GAT'
        v8 = pymummer.variant.Variant(pymummer.snp.Snp('4\t.\tG\t4\tx\tx\t39\t39\tx\tx\tgene\tcontig'))
        v8.qry_base = 'TGA'
        v9 = pymummer.variant.Variant(pymummer.snp.Snp('4\t.\tG\t4\tx\tx\t39\t39\tx\tx\tgene\tcontig'))
        v9.qry_base = 'ATTCCT'
        v10 = pymummer.variant.Variant(pymummer.snp.Snp('4\tC\t.\t4\tx\tx\t39\t39\tx\tx\tgene\tcontig'))
        v10.ref_base = 'CGC'
        v10.ref_end = 5
        v11 = pymummer.variant.Variant(pymummer.snp.Snp('4\tC\t.\t4\tx\tx\t39\t39\tx\tx\tgene\tcontig'))
        v11.ref_base = 'CGCGAA'
        v11.ref_end = 8

        variants = [
            ([v1], ('SYN', '.', 1)),
            ([v2], ('NONSYN', 'R2S', 1)),
            ([v2, v1], ('NONSYN', 'R2S', 1)),
            ([v3, v4], ('TRUNC', 'R2trunc', 1)),
            ([v5], ('FSHIFT', 'R2fs', 1)),
            ([v6], ('FSHIFT', 'R2fs', 1)),
            ([v7], ('INS', 'R2_E3insD', 1)),
            ([v8], ('TRUNC', 'R2trunc', 1)),
            ([v9], ('INS', 'R2_E3insIP', 1)),
            ([v10], ('DEL', 'R2del', 1)),
            ([v11], ('DEL', 'R2_E3del', 1)),
        ]

        for variant_list, expected in variants:
            self.assertEqual(expected, assembly_variants.AssemblyVariants._get_variant_effect(variant_list, ref_seq))


    def test_filter_mummer_variants(self):
        '''test filter_mummer_variants'''
        ref_seq = pyfastaq.sequences.Fasta('gene', 'GATCGCGAAGCGATGACCCATGAAGCGACCGAACGCTGA')
        v1 = pymummer.variant.Variant(pymummer.snp.Snp('4\tC\tT\t4\tx\tx\t39\t39\tx\tx\tgene\tcontig'))
        v2 = pymummer.variant.Variant(pymummer.snp.Snp('6\tC\tA\t6\tx\tx\t39\t39\tx\tx\tgene\tcontig'))
        v3 = pymummer.variant.Variant(pymummer.snp.Snp('12\tG\tT\t12\tx\tx\t39\t39\tx\tx\tgene\tcontig'))
        mummer_variants = {'contig': [[v1, v2], v3]}
        assembly_variants.AssemblyVariants._filter_mummer_variants(mummer_variants, ref_seq)
        expected = {'contig': [[v1, v2]]}
        self.assertEqual(expected, mummer_variants)


    def test_one_var_one_ctg_noncdg(self):
        '''test _get_one_variant_for_one_contig_non_coding'''
        fasta_in = os.path.join(data_dir, 'assembly_variants_one_var_one_ctg_noncdg.fa')
        tsv_in = os.path.join(data_dir, 'assembly_variants_one_var_one_ctg_noncdg.tsv')
        refdata = reference_data.ReferenceData([fasta_in], [tsv_in])
        ref_sequence_name = 'non_coding'
        refdata_var_dict = refdata.metadata[ref_sequence_name]

        v0 = pymummer.variant.Variant(pymummer.snp.Snp('2\tT\tA\t2\tx\tx\t42\t42\tx\tx\tnon_coding\tcontig'))

        # ref has A at position 3, which is variant type. This gives contig the wild type C. Shouldn't report
        v1 = pymummer.variant.Variant(pymummer.snp.Snp('3\tA\tC\t3\tx\tx\t42\t42\tx\tx\tnon_coding\tcontig'))

        # ref has T at position 5, which is wild type. This gives contig variant type A. Should report
        v2 = pymummer.variant.Variant(pymummer.snp.Snp('5\tT\tA\t5\tx\tx\t42\t42\tx\tx\tnon_coding\tcontig'))

        meta0 = sequence_metadata.SequenceMetadata('non_coding\t0\t0\tC3A\tid1\tref has variant type A')
        meta2 = sequence_metadata.SequenceMetadata('non_coding\t0\t0\tT5A\tid1\tref has wild type T')

        mummer_variants = [v0, v1, v2]

        expected_tuples = [
            (1, 'n', 'T2A', 'SNP', [v0], set(), set()),   #0
            None,                                     #1
            (4, 'n', 'T5A', 'SNP', [v2], {meta2}, set()), #2
        ]

        expected_used_variants = [
            set(),     #0
            {meta0},   #1
            {meta2},   #2
        ]

        assert len(mummer_variants) == len(expected_tuples) == len(expected_used_variants)


        for i in range(len(mummer_variants)):
            got_tuple, got_used_variants = assembly_variants.AssemblyVariants._get_one_variant_for_one_contig_non_coding(refdata_var_dict, mummer_variants[i])
            self.assertEqual(expected_tuples[i], got_tuple)
            self.assertEqual(expected_used_variants[i], got_used_variants)


    def test_one_var_one_ctg_cdg(self):
        '''test _get_one_variant_for_one_contig_coding'''
        fasta_in = os.path.join(data_dir, 'assembly_variants_one_var_one_ctg_cdg.fa')
        tsv_in = os.path.join(data_dir, 'assembly_variants_one_var_one_ctg_cdg.tsv')
        refdata = reference_data.ReferenceData([fasta_in], [tsv_in])
        ref_sequence_name = 'presence_absence'
        ref_sequence = refdata.sequence(ref_sequence_name)
        refdata_var_dict = refdata.metadata[ref_sequence_name]

        v0 = pymummer.variant.Variant(pymummer.snp.Snp('6\tT\tA\t6\tx\tx\t42\t42\tx\tx\tpresence_absence\tcontig'))
        v1 = pymummer.variant.Variant(pymummer.snp.Snp('9\tA\tT\t9\tx\tx\t42\t42\tx\tx\tpresence_absence\tcontig'))
        v2 = pymummer.variant.Variant(pymummer.snp.Snp('18\tG\tT\t18\tx\tx\t42\t42\tx\tx\tpresence_absence\tcontig'))
        v3 = pymummer.variant.Variant(pymummer.snp.Snp('21\tC\tT\t21\tx\tx\t42\t42\tx\tx\tpresence_absence\tcontig'))
        v4 = pymummer.variant.Variant(pymummer.snp.Snp('7\tA\tT\t7\tx\tx\t42\t42\tx\tx\tpresence_absence\tcontig'))
        v5 = pymummer.variant.Variant(pymummer.snp.Snp('12\tA\tC\t11\tx\tx\t42\t42\tx\tx\tpresence_absence\tcontig'))

        v6 = pymummer.variant.Variant(pymummer.snp.Snp('4\tG\t.\t4\tx\tx\t42\t42\tx\tx\tpresence_absence\tcontig'))
        self.assertTrue(v6.update_indel(pymummer.snp.Snp('5\tA\t.\t4\tx\tx\t42\t42\tx\tx\tpresence_absence\tcontig')))

        v7 = pymummer.variant.Variant(pymummer.snp.Snp('4\t.\tA\t4\tx\tx\t42\t42\tx\tx\tpresence_absence\tcontig'))
        self.assertTrue(v7.update_indel(pymummer.snp.Snp('4\t.\tA\t5\tx\tx\t42\t42\tx\tx\tpresence_absence\tcontig')))

        v8 = pymummer.variant.Variant(pymummer.snp.Snp('4\tG\t.\t4\tx\tx\t42\t42\tx\tx\tpresence_absence\tcontig'))
        self.assertTrue(v8.update_indel(pymummer.snp.Snp('5\tA\t.\t4\tx\tx\t42\t42\tx\tx\tpresence_absence\tcontig')))
        self.assertTrue(v8.update_indel(pymummer.snp.Snp('6\tT\t.\t4\tx\tx\t42\t42\tx\tx\tpresence_absence\tcontig')))

        v9 = pymummer.variant.Variant(pymummer.snp.Snp('4\tG\t.\t4\tx\tx\t42\t42\tx\tx\tpresence_absence\tcontig'))
        self.assertTrue(v9.update_indel(pymummer.snp.Snp('5\tA\t.\t4\tx\tx\t42\t42\tx\tx\tpresence_absence\tcontig')))
        self.assertTrue(v9.update_indel(pymummer.snp.Snp('6\tT\t.\t4\tx\tx\t42\t42\tx\tx\tpresence_absence\tcontig')))
        self.assertTrue(v9.update_indel(pymummer.snp.Snp('7\tA\t.\t4\tx\tx\t42\t42\tx\tx\tpresence_absence\tcontig')))
        self.assertTrue(v9.update_indel(pymummer.snp.Snp('8\tG\t.\t4\tx\tx\t42\t42\tx\tx\tpresence_absence\tcontig')))
        self.assertTrue(v9.update_indel(pymummer.snp.Snp('9\tA\t.\t4\tx\tx\t42\t42\tx\tx\tpresence_absence\tcontig')))

        v10 = pymummer.variant.Variant(pymummer.snp.Snp('4\t.\tA\t4\tx\tx\t42\t42\tx\tx\tpresence_absence\tcontig'))
        self.assertTrue(v10.update_indel(pymummer.snp.Snp('4\t.\tT\t5\tx\tx\t42\t42\tx\tx\tpresence_absence\tcontig')))
        self.assertTrue(v10.update_indel(pymummer.snp.Snp('4\t.\tT\t6\tx\tx\t42\t42\tx\tx\tpresence_absence\tcontig')))

        mummer_variants = [[v0], [v1], [v2], [v3], [v4], [v5], [v6], [v7], [v8], [v9], [v10]]

        meta0 = sequence_metadata.SequenceMetadata('presence_absence\t1\t0\tD2E\tid1\tref has wild type D (GAT=D, GAA=E)')
        meta4 = sequence_metadata.SequenceMetadata('presence_absence\t1\t0\tS3R\tid1\tref has variant type R (AGA=R, AGT=S)')

        expected_tuples = [
            (1, 'p', 'D2E', 'NONSYN', [v0], {meta0}, set()),    #0
            None,                                               #1
            (5, 'p', 'M6I', 'NONSYN', [v2], set(), set()),      #2
            (6, 'p', '.', 'SYN', [v3], set(), set()),           #3
            (2, 'p', 'R3trunc', 'TRUNC', [v4], set(), {meta4}), #4
            None,                                               #5
            (1, 'p', 'D2fs', 'FSHIFT', [v6], set(), {meta0}),   #6
            (1, 'p', 'D2fs', 'FSHIFT', [v7], set(), {meta0}),   #7
            (1, 'p', 'D2del', 'DEL', [v8], set(), {meta0}),     #8
            (1, 'p', 'D2_R3del', 'DEL', [v9], set(), {meta0}),  #9
            (1, 'p', 'D2_R3insI', 'INS', [v10], set(), {meta0}) #10
        ]

        expected_used_variants = [
            refdata_var_dict['p'][1], #0
            refdata_var_dict['p'][2], #1
            set(),                    #2
            set(),                    #3
            refdata_var_dict['p'][2], #4
            refdata_var_dict['p'][3], #5
            refdata_var_dict['p'][1], #6
            refdata_var_dict['p'][1], #7
            refdata_var_dict['p'][1], #8
            refdata_var_dict['p'][1], #9
            refdata_var_dict['p'][1], #10
        ]

        assert len(mummer_variants) == len(expected_tuples) == len(expected_used_variants)

        for i in range(len(mummer_variants)):
            got_tuple, got_used_variants = assembly_variants.AssemblyVariants._get_one_variant_for_one_contig_coding(ref_sequence, refdata_var_dict, mummer_variants[i])
            self.assertEqual(expected_tuples[i], got_tuple)
            self.assertEqual(expected_used_variants[i], got_used_variants)


    def test_get_remaining_known_ref_variants_amino_acids(self):
        '''test _get_remaining_known_ref_variants with amino acids'''
        ref_var1 = sequence_metadata.SequenceMetadata('gene1\t1\t0\tD2E\tid1\tfoo bar')
        ref_var2 = sequence_metadata.SequenceMetadata('gene1\t1\t0\tD3E\tid1\tfoo bar baz')
        ref_var3 = sequence_metadata.SequenceMetadata('gene1\t1\t0\tD3I\tid1\tfoo bar baz spam')
        ref_var4 = sequence_metadata.SequenceMetadata('gene1\t1\t0\tD10E\tid1\tfoo bar baz spam egg')
        ref_var5 = sequence_metadata.SequenceMetadata('gene1\t1\t0\tD14E\tid1\tfoo bar baz spam egg chips')
        ref_var6 = sequence_metadata.SequenceMetadata('gene1\t1\t0\tD15E\tid1\tfoo bar baz spam egg chips')
        ref_var7 = sequence_metadata.SequenceMetadata('gene1\t1\t0\tD40E\tid1\tfoo bar baz spam egg chips')

        known_ref_variants = {
            1: {ref_var1},
            2: {ref_var2, ref_var3},
            9: {ref_var4},
            13: {ref_var5},
            14: {ref_var6},
            39: {ref_var7}
        }

        used_ref_variants = {ref_var3, ref_var5}

        nucmer_coords = [
            pyfastaq.intervals.Interval(6, 25),
            pyfastaq.intervals.Interval(30, 100)
        ]

        expected = [(None, 'p', None, None, None, {x}, set()) for x in [ref_var2, ref_var6]]
        got = assembly_variants.AssemblyVariants._get_remaining_known_ref_variants(known_ref_variants, used_ref_variants, nucmer_coords)
        self.assertEqual(expected, got)


    def test_get_remaining_known_ref_variants_nucleotides(self):
        '''test _get_remaining_known_ref_variants with nucleotides'''
        ref_var1 = sequence_metadata.SequenceMetadata('gene1\t0\t0\tA2C\tid1\tfoo bar')
        ref_var2 = sequence_metadata.SequenceMetadata('gene1\t0\t0\tA3C\tid1\tfoo bar baz')
        ref_var3 = sequence_metadata.SequenceMetadata('gene1\t0\t0\tA3T\tid1\tfoo bar baz spam')
        ref_var4 = sequence_metadata.SequenceMetadata('gene1\t0\t0\tA10C\tid1\tfoo bar baz spam egg')
        ref_var5 = sequence_metadata.SequenceMetadata('gene1\t0\t0\tA14C\tid1\tfoo bar baz spam egg chips')
        ref_var6 = sequence_metadata.SequenceMetadata('gene1\t0\t0\tA15C\tid1\tfoo bar baz spam egg chips')
        ref_var7 = sequence_metadata.SequenceMetadata('gene1\t0\t0\tA40C\tid1\tfoo bar baz spam egg chips')

        known_ref_variants = {
            1: {ref_var1},
            2: {ref_var2, ref_var3},
            9: {ref_var4},
            13: {ref_var5},
            14: {ref_var6},
            39: {ref_var7}
        }

        used_ref_variants = {ref_var3, ref_var5}

        nucmer_coords = [
            pyfastaq.intervals.Interval(2, 13),
            pyfastaq.intervals.Interval(30, 100)
        ]

        expected = [(None, 'n', None, None, None, {x}, set()) for x in [ref_var2, ref_var4, ref_var7]]
        got = assembly_variants.AssemblyVariants._get_remaining_known_ref_variants(known_ref_variants, used_ref_variants, nucmer_coords)
        self.assertEqual(expected, got)


    def test_get_variants_presence_absence(self):
        '''test get_variants presence absence genes'''
        meta1 = sequence_metadata.SequenceMetadata('presence_absence\t1\t0\tD2E\tid1\tref has wild type D, contig has var (GAT=D, GAA=E)')
        meta2 = sequence_metadata.SequenceMetadata('presence_absence\t1\t0\tS3R\tid1\tref has variant type R, contig has wild (AGA=R, AGT=S)')
        meta3 = sequence_metadata.SequenceMetadata('presence_absence\t1\t0\tD4E\tid1\tref has variant type E, contig has var (GAA=E, GAC=D)')
        meta4 = sequence_metadata.SequenceMetadata('presence_absence\t1\t0\tA5D\tid1\tref has wild type A, contig has var (GCG=A, GAC=D)')
        meta5 = sequence_metadata.SequenceMetadata('presence_absence\t1\t0\tR13S\tid1\tref and qry have wild type')

        metadata_tsv = 'tmp.test_get_variants_presence_absence.metadata.tsv'
        with open(metadata_tsv, 'w') as f:
            print(meta1, file=f)
            print(meta2, file=f)
            print(meta3, file=f)
            print(meta4, file=f)
            print(meta5, file=f)

        fasta_in = os.path.join(data_dir, 'assembly_variants_test_get_variants_presence_absence.fa')
        refdata = reference_data.ReferenceData([fasta_in], [metadata_tsv])

        os.unlink(metadata_tsv)

        nucmer_snp_file = os.path.join(data_dir, 'assembly_variants_test_get_variants_presence_absence.snps')
        v2 = pymummer.variant.Variant(pymummer.snp.Snp('14\tC\tA\t14\tx\tx\t42\t42\tx\tx\tpresence_absence\tcontig1'))
        v3 = pymummer.variant.Variant(pymummer.snp.Snp('15\tG\tC\t15\tx\tx\t42\t42\tx\tx\tpresence_absence\tcontig1'))

        ref_nucmer_coords = {
            'contig1': [pyfastaq.intervals.Interval(0, 30)],
            'contig2': [pyfastaq.intervals.Interval(10, 41)],
        }

        ctg_nucmer_coords = {
            'contig1': [pyfastaq.intervals.Interval(0, 30)],
            'contig2': [pyfastaq.intervals.Interval(10, 41)],
        }

        expected = {
            'contig1': [
               (4, 'p', 'A5D', 'NONSYN', [v2, v3], {meta4}, set()),
               (None, 'p', None, None, None, {meta1}, set()),
               (None, 'p', None, None, None, {meta3}, set()),
            ],
            'contig2': [
               (None, 'p', None, None, None, {meta3}, set()),
               (None, 'p', None, None, None, {meta4}, set()),
               (None, 'p', None, None, None, {meta5}, set()),
            ],
        }

        a_variants = assembly_variants.AssemblyVariants(refdata, nucmer_snp_file)
        got = a_variants.get_variants('presence_absence', ctg_nucmer_coords, ref_nucmer_coords)
        self.assertEqual(expected, got)


    def test_get_variants_variants_only(self):
        '''test get_variants variants only'''
        meta1 = sequence_metadata.SequenceMetadata('variants_only\t1\t0\tD2E\tid1\tref has wild type D (GAT=D, GAA=E)')
        meta2 = sequence_metadata.SequenceMetadata('variants_only\t1\t0\tS3R\tid1\tref has variant type R (AGA=R, AGT=S)')
        meta3 = sequence_metadata.SequenceMetadata('variants_only\t1\t0\tD4E\tid1\tref has variant type E (GAA=E, GAC=D)')

        metadata_tsv = 'tmp.test_get_variants_variants_only.metadata.tsv'
        with open(metadata_tsv, 'w') as f:
            print(meta1, file=f)
            print(meta2, file=f)
            print(meta3, file=f)

        fasta_in = os.path.join(data_dir, 'assembly_variants_test_get_variants_variants_only.fa')
        refdata = reference_data.ReferenceData([fasta_in], [metadata_tsv])
        os.unlink(metadata_tsv)

        nucmer_snp_file = os.path.join(data_dir, 'assembly_variants_test_get_variants_variants_only.snps')
        v2 = pymummer.variant.Variant(pymummer.snp.Snp('14\tC\tA\t14\tx\tx\t42\t42\tx\tx\tvariants_only\tcontig1'))
        v3 = pymummer.variant.Variant(pymummer.snp.Snp('15\tG\tC\t15\tx\tx\t42\t42\tx\tx\tvariants_only\tcontig1'))

        ctg_nucmer_coords = {
            'contig1': [pyfastaq.intervals.Interval(0, 41)],
            'contig2': [pyfastaq.intervals.Interval(10, 41)],
        }

        ref_nucmer_coords = {
            'contig1': [pyfastaq.intervals.Interval(0, 41)],
            'contig2': [pyfastaq.intervals.Interval(10, 41)],
        }
        expected = {
            'contig1': [
                (4, 'p', 'A5D', 'NONSYN', [v2, v3], set(), set()),
                (None, 'p', None, None, None, {meta1}, set()),
                (None, 'p', None, None, None, {meta3}, set()),
            ],
            'contig2': [(None, 'p', None, None, None, {meta3}, set())],
        }

        a_variants = assembly_variants.AssemblyVariants(refdata, nucmer_snp_file)
        got = a_variants.get_variants('variants_only', ctg_nucmer_coords, ref_nucmer_coords)
        self.assertEqual(expected, got)

