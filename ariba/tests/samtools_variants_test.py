import unittest
import os
from ariba import samtools_variants

modules_dir = os.path.dirname(os.path.abspath(samtools_variants.__file__))


class TestSamtoolsVariants(unittest.TestCase):
    def test_make_vcf_and_read_depths_files(self):
        '''test _make_vcf_and_read_depths_files'''
        cluster_dir = os.path.join(data_dir, 'cluster_test_generic')
        clean_cluster_dir(cluster_dir)
        c = cluster.Cluster(cluster_dir, 'name')
        c.final_assembly_fa = os.path.join(data_dir, 'cluster_test_make_assembly_vcf.assembly.fa')
        c.final_assembly_bam = os.path.join(data_dir, 'cluster_test_make_assembly_vcf.assembly.bam')
        expected_vcf = os.path.join(data_dir, 'cluster_test_make_assembly_vcf.assembly.vcf')
        expected_depths = os.path.join(data_dir, 'cluster_test_make_assembly_vcf.assembly.read_depths.gz')
        c._make_assembly_vcf()

        def get_vcf_call_lines(fname):
            with open(fname) as f:
                lines = [x for x in f.readlines() if not x.startswith('#')]
            return lines

        expected_lines = get_vcf_call_lines(expected_vcf)
        got_lines = get_vcf_call_lines(c.final_assembly_vcf)
        self.assertEqual(expected_lines, got_lines)
        self.assertEqual(file2lines(expected_depths), file2lines(c.final_assembly_read_depths))
        clean_cluster_dir(cluster_dir)


    def test_get_read_depths(self):
        '''test _get_read_depths'''
        cluster_dir = os.path.join(data_dir, 'cluster_test_generic')
        clean_cluster_dir(cluster_dir)
        c = cluster.Cluster(cluster_dir, 'name')
        c.final_assembly_read_depths = os.path.join(data_dir, 'cluster_test_get_assembly_read_depths.gz')
        tests = [
            ( ('ref1', 42), None ),
            ( ('ref2', 1), None ),
            ( ('ref1', 0), ('G', '.', 1, '1') ),
            ( ('ref1', 2), ('T', 'A', 3, '2,1') ),
            ( ('ref1', 3), ('C', 'A,G', 42, '21,11,10') ),
            ( ('ref1', 4), ('C', 'AC', 41, '0,42') )
        ]

        for t in tests:
            self.assertEqual(c._get_assembly_read_depths(t[0][0], t[0][1]), t[1])


    def test_get_variant_positions_from_vcf(self):
        '''test _get_variant_positions_from_vcf'''
        cluster_dir = os.path.join(data_dir, 'cluster_test_generic')
        clean_cluster_dir(cluster_dir)
        c = cluster.Cluster(cluster_dir, 'name')
        c.final_assembly_vcf = os.path.join(data_dir, 'cluster_test_get_samtools_variant_positions.vcf')
        expected = [
            ('16__cat_2_M35190.scaffold.1', 92),
            ('16__cat_2_M35190.scaffold.1', 179),
            ('16__cat_2_M35190.scaffold.1', 263),
            ('16__cat_2_M35190.scaffold.6', 93)
        ]
        self.assertEqual(expected, c._get_samtools_variant_positions())


    def test_get_variants(self):
        '''test _get_variants'''
        cluster_dir = os.path.join(data_dir, 'cluster_test_generic')
        clean_cluster_dir(cluster_dir)
        c = cluster.Cluster(cluster_dir, 'name')
        c.final_assembly_vcf = os.path.join(data_dir, 'cluster_test_get_samtools_variants.vcf')
        c.final_assembly_read_depths = os.path.join(data_dir, 'cluster_test_get_samtools_variants.read_depths.gz')
        positions = [
            ('16__cat_2_M35190.scaffold.1', 92),
            ('16__cat_2_M35190.scaffold.1', 179),
            ('16__cat_2_M35190.scaffold.1', 263),
            ('16__cat_2_M35190.scaffold.6', 93)
        ]
        expected = {
            '16__cat_2_M35190.scaffold.1': {
                92: ('T', 'A', 123, '65,58'),
                179: ('A', 'T', 86, '41,45'),
                263: ('G', 'C', 97, '53,44'),
            },
            '16__cat_2_M35190.scaffold.6': {
                93: ('T', 'G', 99, '56,43')
            }
        }

        got = c._get_samtools_variants(positions)
        self.assertEqual(expected, got)


    def test_variants_in_coords(self):
        '''test _get_variants_in_coords'''
        cluster_dir = os.path.join(data_dir, 'cluster_test_generic')
        clean_cluster_dir(cluster_dir)
        c = cluster.Cluster(cluster_dir, 'name')
        hit = ['1', '42', '1', '42', '42', '42', '100.00', '1000', '1000', '1', '1', 'gene', 'scaff1']
        c.nucmer_hits = {
            'scaff1': [pymummer.alignment.Alignment('\t'.join(hit))]
        }

        c.final_assembly_vcf = os.path.join(data_dir, 'cluster_test_get_vcf_variant_counts.vcf')
        c._get_vcf_variant_counts()
        expected = {'scaff1': 1}
        self.assertEqual(expected, c.vcf_variant_counts)
        clean_cluster_dir(cluster_dir)


