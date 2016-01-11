import unittest
from ariba import assembly_variants

modules_dir = os.path.dirname(os.path.abspath(assembly_variants.__file__))


class TestAssemblyVariants(unittest.TestCase):
    def test_get_codon_start(self):
        '''test _get_codon_start'''
        cluster_dir = os.path.join(data_dir, 'cluster_test_generic')
        clean_cluster_dir(cluster_dir)
        c = cluster.Cluster(cluster_dir, 'name')
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
        for t in tests:
            self.assertEqual(c._get_codon_start(t[0], t[1]), t[2])
        clean_cluster_dir(cluster_dir)


    def test_get_mummer_variants(self):
        '''test _get_mummer_variants'''
        cluster_dir = os.path.join(data_dir, 'cluster_test_generic')
        clean_cluster_dir(cluster_dir)
        c = cluster.Cluster(cluster_dir, 'name')
        snp_file = os.path.join(data_dir, 'cluster_test_get_mummer_variants.none.snps')
        shutil.copyfile(snp_file, c.assembly_vs_gene_coords + '.snps')
        c._get_mummer_variants()
        self.assertEqual(c.mummer_variants, {})

        clean_cluster_dir(cluster_dir)
        snp_file = os.path.join(data_dir, 'cluster_test_get_mummer_variants.snp.snps')
        shutil.copyfile(snp_file, c.assembly_vs_gene_coords + '.snps')
        v1 = pymummer.variant.Variant(pymummer.snp.Snp('42\tA\tG\t42\t42\t42\t500\t500\t1\t1\tgene\tcontig1'))
        v2 = pymummer.variant.Variant(pymummer.snp.Snp('42\tA\tG\t42\t42\t42\t500\t500\t1\t1\tgene\tcontig2'))
        v3 = pymummer.variant.Variant(pymummer.snp.Snp('40\tT\tC\t40\t42\t42\t500\t500\t1\t1\tgene\tcontig1'))
        v4 = pymummer.variant.Variant(pymummer.snp.Snp('2\tC\tG\t2\t42\t42\t500\t500\t1\t1\tgene\tcontig1'))
        expected = {
            'contig1': [[v4], [v3, v1]],
            'contig2': [[v2]]
        }
        shutil.copyfile(snp_file, c.assembly_vs_gene_coords + '.snps')
        c._get_mummer_variants()
        self.assertEqual(c.mummer_variants, expected)
        clean_cluster_dir(cluster_dir)


    def test_filter_mummer_variants(self):
        '''test filter_mummer_variants'''
        cluster_dir = os.path.join(data_dir, 'cluster_test_generic')
        clean_cluster_dir(cluster_dir)
        c = cluster.Cluster(cluster_dir, 'name')
        c.gene = pyfastaq.sequences.Fasta('gene', 'GATCGCGAAGCGATGACCCATGAAGCGACCGAACGCTGA')
        v1 = pymummer.variant.Variant(pymummer.snp.Snp('4\tC\tT\t4\tx\tx\t39\t39\tx\tx\tgene\tcontig'))
        v2 = pymummer.variant.Variant(pymummer.snp.Snp('6\tC\tA\t6\tx\tx\t39\t39\tx\tx\tgene\tcontig'))
        v3 = pymummer.variant.Variant(pymummer.snp.Snp('12\tG\tT\t12\tx\tx\t39\t39\tx\tx\tgene\tcontig'))
        c.mummer_variants = {'contig': [[v1, v2], v3]}
        c._filter_mummer_variants()
        expected = {'contig': [[v1, v2]]}
        self.assertEqual(expected, c.mummer_variants)
        clean_cluster_dir(cluster_dir)


    def test_get_variant_effect(self):
        '''test _get_variant_effect'''
        cluster_dir = os.path.join(data_dir, 'cluster_test_generic')
        clean_cluster_dir(cluster_dir)
        c = cluster.Cluster(cluster_dir, 'name')
        c.gene = pyfastaq.sequences.Fasta('gene', 'GATCGCGAAGCGATGACCCATGAAGCGACCGAACGCTGA')
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
            ([v1], ('SYN', '.')),
            ([v2], ('NONSYN', 'R2S')),
            ([v2, v1], ('NONSYN', 'R2S')),
            ([v3, v4], ('TRUNC', 'R2trunc')),
            ([v5], ('FSHIFT', 'R2fs')),
            ([v6], ('FSHIFT', 'R2fs')),
            ([v7], ('INS', 'R2_E3insD')),
            ([v8], ('TRUNC', 'R2trunc')),
            ([v9], ('INS', 'R2_E3insIP')),
            ([v10], ('DEL', 'R2del')),
            ([v11], ('DEL', 'R2_E3del')),
        ]

        for t in variants:
            self.assertEqual(t[1], c._get_variant_effect(t[0]))

        clean_cluster_dir(cluster_dir)

