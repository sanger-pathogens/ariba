import unittest
import os
import filecmp
import pyfastaq
from ariba import samtools_variants

modules_dir = os.path.dirname(os.path.abspath(samtools_variants.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


def file2lines(filename):
    f = pyfastaq.utils.open_file_read(filename)
    lines = f.readlines()
    pyfastaq.utils.close(f)
    return lines


class TestSamtoolsVariants(unittest.TestCase):
    def test_make_vcf_and_depths_files(self):
        '''test _make_vcf_and_read_depths_files'''
        ref = os.path.join(data_dir, 'samtools_variants_make_vcf_and_depths_files.asmbly.fa')
        bam = os.path.join(data_dir, 'samtools_variants_make_vcf_and_depths_files.bam')
        expected_vcf = os.path.join(data_dir, 'samtools_variants_make_vcf_and_depths_files.expect.vcf')
        expected_depths = os.path.join(data_dir, 'samtools_variants_make_vcf_and_depths_files.expect.depths.gz')
        expected_coverage = os.path.join(data_dir, 'samtools_variants_make_vcf_and_depths_files.expect.cov')
        tmp_prefix = 'tmp.test_make_vcf_and_depths_files'
        sv = samtools_variants.SamtoolsVariants(
            ref,
            bam,
            tmp_prefix,
        )
        sv._make_vcf_and_read_depths_files()

        def get_vcf_call_lines(fname):
            with open(fname) as f:
                lines = [x for x in f.readlines() if not x.startswith('#')]
            return lines

        expected_lines = get_vcf_call_lines(expected_vcf)
        got_lines = get_vcf_call_lines(sv.vcf_file)

        # need to check that the vcf lines look the same. Column 8 is ;-delimited
        # and can be an any order.
        self.assertEqual(len(expected_lines), len(got_lines))

        for i in range(len(expected_lines)):
            expected = expected_lines[i].split('\t')
            got = got_lines[i].split('\t')
            self.assertEqual(len(expected), len(got))
            self.assertEqual(expected[:7], got[:7])
            self.assertEqual(expected[-2:], got[-2:])
            exp_set = set(expected[7].split(';'))
            got_set = set(got[7].split(';'))
            self.assertEqual(exp_set, got_set)


        # samtools-1.2 and 1.3 output not xonsistent in final column, so
        # ignore those.
        expected_lines = file2lines(expected_depths)
        got_lines = file2lines(sv.read_depths_file)
        self.assertEqual(len(expected_lines), len(got_lines))

        for i in range(len(expected_lines)):
            self.assertEqual(expected_lines[i].split('\t')[:-1], got_lines[i].split('\t')[:-1])

        os.unlink(sv.vcf_file)
        os.unlink(sv.read_depths_file)
        os.unlink(sv.read_depths_file + '.tbi')

        self.assertTrue(filecmp.cmp(expected_coverage, tmp_prefix + '.contig_depths', shallow=False))
        os.unlink(tmp_prefix + '.contig_depths')


    def test_get_read_depths(self):
        '''test _get_read_depths'''
        read_depths_file = os.path.join(data_dir, 'samtools_variants_test_get_read_depths.gz')

        tests = [
            ( ('ref1', 42), None ),
            ( ('ref2', 1), None ),
            ( ('ref1', 0), ('G', 1, '1') ),
            ( ('ref1', 2), ('T,A', 3, '2,1') ),
            ( ('ref1', 3), ('C,A,G', 42, '21,11,10') ),
            ( ('ref1', 4), ('C,AC', 41, '0,42') )
        ]

        for (name, position), expected in tests:
            self.assertEqual(expected, samtools_variants.SamtoolsVariants._get_read_depths(read_depths_file, name, position))


    def test_get_variant_positions_from_vcf(self):
        '''test _get_variant_positions_from_vcf'''
        vcf_file = os.path.join(data_dir, 'samtools_variants_test_get_variant_positions_from_vcf.vcf')

        expected = [
            ('16__cat_2_M35190.scaffold.1', 92),
            ('16__cat_2_M35190.scaffold.1', 179),
            ('16__cat_2_M35190.scaffold.1', 263),
            ('16__cat_2_M35190.scaffold.6', 93)
        ]
        self.assertEqual(expected, samtools_variants.SamtoolsVariants._get_variant_positions_from_vcf(vcf_file))


    def test_get_variants(self):
        '''test _get_variants'''
        vcf_file = os.path.join(data_dir, 'samtools_variants_test_get_variants.vcf')
        read_depths_file = os.path.join(data_dir, 'samtools_variants_test_get_variants.read_depths.gz')
        positions = [
            ('16__cat_2_M35190.scaffold.1', 92),
            ('16__cat_2_M35190.scaffold.1', 179),
            ('16__cat_2_M35190.scaffold.1', 263),
            ('16__cat_2_M35190.scaffold.6', 93)
        ]
        expected = {
            '16__cat_2_M35190.scaffold.1': {
                92: ('T,A',123, '65,58'),
                179: ('A,T', 86, '41,45'),
                263: ('G,C', 97, '53,44'),
            },
            '16__cat_2_M35190.scaffold.6': {
                93: ('T,G', 99, '56,43')
            }
        }

        got = samtools_variants.SamtoolsVariants._get_variants(vcf_file, read_depths_file, positions=positions)
        self.assertEqual(expected, got)


    def test_total_depth_per_contig(self):
        '''test total_depth_per_contig'''
        infile = os.path.join(data_dir, 'samtools_variants_test_total_depth_per_contig')
        expected = {'scaff1': 42, 'scaff64738': 11}
        got = samtools_variants.SamtoolsVariants.total_depth_per_contig(infile)
        self.assertEqual(expected, got)


    def test_variants_in_coords(self):
        '''test variants_in_coords'''
        vcf_file = os.path.join(data_dir, 'samtools_variants_test_variants_in_coords.vcf')

        nucmer_hits = {
            'scaff1': [pyfastaq.intervals.Interval(0, 41)]
        }

        expected = {'scaff1': {22}}
        got = samtools_variants.SamtoolsVariants.variants_in_coords(nucmer_hits, vcf_file)
        self.assertEqual(expected, got)


    def test_get_depths_at_position(self):
        '''test get_depths_at_position'''
        bam = os.path.join(data_dir, 'samtools_variants_test_get_depths_at_position.bam')
        ref_fa = os.path.join(data_dir, 'samtools_variants_test_get_depths_at_position.ref.fa')
        tmp_prefix = 'tmp.test_get_depths_at_position'
        samtools_vars = samtools_variants.SamtoolsVariants(
            ref_fa,
            bam,
            tmp_prefix,
        )
        samtools_vars.run()
        tests = [
            (('ref', 425), ('C,T', 31, '18,13')),
            (('not_a_ref', 10), ('ND', 'ND', 'ND')),
            (('ref', 1000000000), ('ND', 'ND', 'ND'))
        ]
        for (ref, pos), expected in tests:
            got = samtools_vars.get_depths_at_position(ref, pos)
            self.assertEqual(expected, got)

        os.unlink(samtools_vars.vcf_file)
        os.unlink(samtools_vars.read_depths_file)
        os.unlink(samtools_vars.read_depths_file + '.tbi')
        os.unlink(samtools_vars.contig_depths_file)
