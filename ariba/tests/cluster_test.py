import unittest
import os
import shutil
import filecmp
from ariba import cluster, common, reference_data

modules_dir = os.path.dirname(os.path.abspath(cluster.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')
cluster.unittest = True


def clean_cluster_dir(d, exclude=None):
    if not os.path.exists(d):
        return

    '''Cleans up all files made except original ones in a cluster directory'''
    keep = set(['genes.fa', 'reads_1.fq', 'reads_2.fq'])
    if exclude is not None:
        for f in exclude:
            keep.add(f)

    for name in os.listdir(d):
        if name not in keep:
            full_path = os.path.join(d, name)
            if os.path.isdir(full_path):
                common.rmtree(full_path)
            else:
                os.unlink(full_path)


class TestCluster(unittest.TestCase):
    def test_init_fail_files_missing(self):
        '''test init_fail_files_missing'''
        refdata_fa = os.path.join(data_dir, 'cluster_test_init_refdata.fa')
        meatadata_tsv = os.path.join(data_dir, 'cluster_test_init_refdata.tsv')
        refdata = reference_data.ReferenceData([refdata_fa], [meatadata_tsv])

        dirs = [
            'cluster_test_init_no_refs_fa',
            'cluster_test_init_no_reads_1',
            'cluster_test_init_no_reads_2',
        ]
        dirs = [os.path.join(data_dir, d) for d in dirs]
        for d in dirs:
            tmpdir = 'tmp.cluster_test_init_fail_files_missing'
            common.rmtree(tmpdir)
            shutil.copytree(d, tmpdir)
            with self.assertRaises(cluster.Error):
                cluster.Cluster(tmpdir, 'name', refdata=refdata, total_reads=42, total_reads_bases=4242)
            common.rmtree(tmpdir)


    def test_number_of_reads_for_assembly(self):
        '''Test _number_of_reads_for_assembly'''
        tests = [
            (50, 1000, 10, 20, 40),
            (50, 999, 10, 20, 42),
            (50, 1000, 10, 10, 20),
            (50, 1000, 10, 5, 10),
        ]

        for insert, bases, reads, coverage, expected in tests:
            self.assertEqual(expected, cluster.Cluster._number_of_reads_for_assembly(100, insert, bases, reads, coverage))


    def test_make_reads_for_assembly_proper_sample(self):
        '''Test _make_reads_for_assembly when sampling from reads'''
        reads_in1 = os.path.join(data_dir, 'cluster_test_make_reads_for_assembly.in1.fq')
        reads_in2 = os.path.join(data_dir, 'cluster_test_make_reads_for_assembly.in2.fq')
        expected_out1 = os.path.join(data_dir, 'cluster_test_make_reads_for_assembly.out1.fq')
        expected_out2 = os.path.join(data_dir, 'cluster_test_make_reads_for_assembly.out2.fq')
        reads_out1 = 'tmp.test_make_reads_for_assembly.reads.out1.fq'
        reads_out2 = 'tmp.test_make_reads_for_assembly.reads.out2.fq'
        reads_written = cluster.Cluster._make_reads_for_assembly(10, 20, reads_in1, reads_in2, reads_out1, reads_out2, random_seed=42)
        self.assertEqual(14, reads_written)
        self.assertTrue(filecmp.cmp(expected_out1, reads_out1, shallow=False))
        self.assertTrue(filecmp.cmp(expected_out2, reads_out2, shallow=False))
        os.unlink(reads_out1)
        os.unlink(reads_out2)


    def test_make_reads_for_assembly_symlinks(self):
        '''Test _make_reads_for_assembly when just makes symlinks'''
        reads_in1 = os.path.join(data_dir, 'cluster_test_make_reads_for_assembly.in1.fq')
        reads_in2 = os.path.join(data_dir, 'cluster_test_make_reads_for_assembly.in2.fq')
        reads_out1 = 'tmp.test_make_reads_for_assembly.reads.out1.fq'
        reads_out2 = 'tmp.test_make_reads_for_assembly.reads.out2.fq'
        reads_written = cluster.Cluster._make_reads_for_assembly(20, 20, reads_in1, reads_in2, reads_out1, reads_out2, random_seed=42)
        self.assertEqual(20, reads_written)
        self.assertTrue(os.path.islink(reads_out1))
        self.assertTrue(os.path.islink(reads_out2))
        self.assertEqual(os.readlink(reads_out1), reads_in1)
        self.assertEqual(os.readlink(reads_out2), reads_in2)
        os.unlink(reads_out1)
        os.unlink(reads_out2)


    def test_full_run_no_reads_after_filtering(self):
        '''test complete run of cluster when filtering removes all reads'''
        fasta_in = os.path.join(data_dir, 'cluster_test_full_run_no_reads_after_filtering.in.fa')
        tsv_in = os.path.join(data_dir, 'cluster_test_full_run_no_reads_after_filtering.in.tsv')
        refdata = reference_data.ReferenceData([fasta_in], [tsv_in])
        tmpdir = 'tmp.test_full_run_no_reads_after_filtering'
        common.rmtree(tmpdir)
        shutil.copytree(os.path.join(data_dir, 'cluster_test_full_run_no_reads_after_filtering'), tmpdir)

        c = cluster.Cluster(tmpdir, 'cluster_name', refdata, total_reads=0, total_reads_bases=0)
        c.run()

        expected = '\t'.join(['.', '.', '.', '.', '64', '0', 'cluster_name'] + ['.'] * 24)
        self.assertEqual([expected], c.report_lines)
        self.assertFalse(c.status_flag.has('ref_seq_choose_fail'))
        self.assertTrue(c.status_flag.has('assembly_fail'))
        common.rmtree(tmpdir)


    def test_full_run_choose_ref_fail(self):
        '''test complete run of cluster when choosing ref seq fails'''
        fasta_in = os.path.join(data_dir, 'cluster_test_full_run_choose_ref_fail.in.fa')
        tsv_in = os.path.join(data_dir, 'cluster_test_full_run_choose_ref_fail.in.tsv')
        refdata = reference_data.ReferenceData([fasta_in], [tsv_in])
        tmpdir = 'tmp.test_full_run_choose_ref_fail'
        common.rmtree(tmpdir)
        shutil.copytree(os.path.join(data_dir, 'cluster_test_full_run_choose_ref_fail'), tmpdir)

        c = cluster.Cluster(tmpdir, 'cluster_name', refdata, total_reads=2, total_reads_bases=108)
        c.run()

        expected = '\t'.join(['.', '.', '.', '.', '1024', '2', 'cluster_name'] + ['.'] * 24)
        self.assertEqual([expected], c.report_lines)
        self.assertTrue(c.status_flag.has('ref_seq_choose_fail'))
        self.assertFalse(c.status_flag.has('assembly_fail'))
        common.rmtree(tmpdir)


    def test_full_run_ref_not_in_cluster(self):
        '''test complete run of cluster when nearest ref is outside cluster'''
        fasta_in = os.path.join(data_dir, 'cluster_test_full_run_ref_not_in_cluster.in.fa')
        tsv_in = os.path.join(data_dir, 'cluster_test_full_run_ref_not_in_cluster.in.tsv')
        refdata = reference_data.ReferenceData([fasta_in], [tsv_in])
        tmpdir = 'tmp.test_full_run_ref_not_in_cluster'
        all_refs_fa =  os.path.join(data_dir, 'cluster_test_full_run_ref_not_in_cluster.all_refs.fa')
        common.rmtree(tmpdir)
        shutil.copytree(os.path.join(data_dir, 'cluster_test_full_run_ref_not_in_cluster'), tmpdir)

        c = cluster.Cluster(tmpdir, 'cluster_name', refdata, total_reads=72, total_reads_bases=3600, all_ref_seqs_fasta=all_refs_fa)
        c.run()

        expected = '\t'.join(['.', '.', '.', '.', '1024', '72', 'cluster_name'] + ['.'] * 24)
        self.assertEqual([expected], c.report_lines)
        self.assertTrue(c.status_flag.has('ref_seq_choose_fail'))
        self.assertFalse(c.status_flag.has('assembly_fail'))
        common.rmtree(tmpdir)


    def test_full_run_assembly_fail(self):
        '''test complete run of cluster when assembly fails'''
        fasta_in = os.path.join(data_dir, 'cluster_test_full_run_assembly_fail.in.fa')
        tsv_in = os.path.join(data_dir, 'cluster_test_full_run_assembly_fail.in.tsv')
        refdata = reference_data.ReferenceData([fasta_in], [tsv_in])
        tmpdir = 'tmp.test_full_run_assembly_fail'
        common.rmtree(tmpdir)
        shutil.copytree(os.path.join(data_dir, 'cluster_test_full_run_assembly_fail'), tmpdir)

        c = cluster.Cluster(tmpdir, 'cluster_name', refdata, total_reads=4, total_reads_bases=304)
        c.run()

        expected = '\t'.join(['.', '.', '.', '.', '64', '4', 'cluster_name'] + ['.'] * 24)
        self.assertEqual([expected], c.report_lines)
        self.assertFalse(c.status_flag.has('ref_seq_choose_fail'))
        self.assertTrue(c.status_flag.has('assembly_fail'))
        common.rmtree(tmpdir)


    def test_full_run_ok_non_coding(self):
        '''test complete run of cluster on a noncoding sequence'''
        fasta_in = os.path.join(data_dir, 'cluster_test_full_run_ok_non_coding.fa')
        tsv_in = os.path.join(data_dir, 'cluster_test_full_run_ok_non_coding.metadata.tsv')
        refdata = reference_data.ReferenceData([fasta_in], [tsv_in])
        tmpdir = 'tmp.test_full_run_ok_non_coding'
        common.rmtree(tmpdir)
        shutil.copytree(os.path.join(data_dir, 'cluster_test_full_run_ok_non_coding'), tmpdir)

        c = cluster.Cluster(tmpdir, 'cluster_name', refdata, total_reads=72, total_reads_bases=3600)
        c.run()

        self.maxDiff=None
        expected = [
            'noncoding1\tnoncoding1\t0\t0\t531\t72\tcluster_name\t120\t120\t95.87\tcluster_name.l15.c30.ctg.1\t234\t15.4\t1\tSNP\tn\tA14T\t1\tA14T\tSNP\t14\t14\tA\t74\t74\tT\t19\tT\t19\tnoncoding1:0:0:A14T:.:ref has wild type, reads has variant so should report\tgeneric description of noncoding1',
            'noncoding1\tnoncoding1\t0\t0\t531\t72\tcluster_name\t120\t120\t95.87\tcluster_name.l15.c30.ctg.1\t234\t15.4\t0\t.\tn\t.\t0\tG61T\tSNP\t61\t61\tG\t121\t121\tT\t24\tT\t24\t.\tgeneric description of noncoding1',
            'noncoding1\tnoncoding1\t0\t0\t531\t72\tcluster_name\t120\t120\t95.87\tcluster_name.l15.c30.ctg.1\t234\t15.4\t0\t.\tn\t.\t0\t.82C\tINS\t82\t82\tA\t143\t143\tC\t23\tC\t23\t.\tgeneric description of noncoding1',
            'noncoding1\tnoncoding1\t0\t0\t531\t72\tcluster_name\t120\t120\t95.87\tcluster_name.l15.c30.ctg.1\t234\t15.4\t0\t.\tn\t.\t0\tT108.\tDEL\t108\t108\tT\t168\t168\tC\t17\tC\t17\t.\tgeneric description of noncoding1',
            'noncoding1\tnoncoding1\t0\t0\t531\t72\tcluster_name\t120\t120\t95.87\tcluster_name.l15.c30.ctg.1\t234\t15.4\t1\tSNP\tn\tA6G\t1\t.\t.\t6\t6\tG\t66\t66\tG\t19\tG\t19\tnoncoding1:0:0:A6G:.:variant in ref and reads so should report\tgeneric description of noncoding1',
            'noncoding1\tnoncoding1\t0\t0\t531\t72\tcluster_name\t120\t120\t95.87\tcluster_name.l15.c30.ctg.1\t234\t15.4\t1\tSNP\tn\tG9T\t0\t.\t.\t9\t9\tG\t69\t69\tG\t19\tG\t19\tnoncoding1:0:0:G9T:.:wild type in ref and reads\tgeneric description of noncoding1'
        ]

        self.assertEqual(expected, c.report_lines)
        common.rmtree(tmpdir)


    def test_full_run_ok_presence_absence(self):
        '''test complete run of cluster on a presence absence gene'''
        fasta_in = os.path.join(data_dir, 'cluster_test_full_run_ok_presence_absence.fa')
        tsv_in = os.path.join(data_dir, 'cluster_test_full_run_ok_presence_absence.metadata.tsv')
        refdata = reference_data.ReferenceData([fasta_in], [tsv_in])
        tmpdir = 'tmp.cluster_test_full_run_ok_presence_absence'
        common.rmtree(tmpdir)
        shutil.copytree(os.path.join(data_dir, 'cluster_test_full_run_ok_presence_absence'), tmpdir)

        c = cluster.Cluster(tmpdir, 'cluster_name', refdata, total_reads=64, total_reads_bases=3200)
        c.run()

        expected = [
            'presence_absence1\tpresence_absence1\t1\t0\t539\t64\tcluster_name\t96\t96\t97.92\tcluster_name.l15.c30.ctg.1\t213\t15.0\t1\tSNP\tp\tA10V\t1\tA10V\tNONSYN\t28\t30\tGCG\t83\t85\tGTG\t22;22;21\tG;T;G\t22;22;21\tpresence_absence1:1:0:A10V:.:Ref has wild, reads have variant so report\tGeneric description of presence_absence1',
            'presence_absence1\tpresence_absence1\t1\t0\t539\t64\tcluster_name\t96\t96\t97.92\tcluster_name.l15.c30.ctg.1\t213\t15.0\t0\t.\tp\t.\t0\t.\tSYN\t52\t54\tATT\t107\t109\tATC\t31;31;32\tA;T;C\t31;31;32\t.\tGeneric description of presence_absence1',
            'presence_absence1\tpresence_absence1\t1\t0\t539\t64\tcluster_name\t96\t96\t97.92\tcluster_name.l15.c30.ctg.1\t213\t15.0\t1\tSNP\tp\tR3S\t0\t.\t.\t7\t9\tCGC\t62\t64\tCGC\t18;17;17\tC;G;C\t18;17;17\tpresence_absence1:1:0:R3S:.:Ref and assembly have wild type\tGeneric description of presence_absence1',
            'presence_absence1\tpresence_absence1\t1\t0\t539\t64\tcluster_name\t96\t96\t97.92\tcluster_name.l15.c30.ctg.1\t213\t15.0\t1\tSNP\tp\tI5A\t1\t.\t.\t13\t15\tGCG\t68\t70\tGCG\t18;20;20\tG;C;G\t18;20;20\tpresence_absence1:1:0:I5A:.:Ref and reads have variant so report\tGeneric description of presence_absence1',
        ]

        self.assertEqual(expected, c.report_lines)
        common.rmtree(tmpdir)


    def test_full_run_ok_variants_only_variant_not_present(self):
        '''test complete run of cluster on a variants only gene when variant not present'''
        fasta_in = os.path.join(data_dir, 'cluster_test_full_run_ok_variants_only.fa')
        tsv_in = os.path.join(data_dir, 'cluster_test_full_run_ok_variants_only.not_present.metadata.tsv')
        refdata = reference_data.ReferenceData([fasta_in], [tsv_in])
        tmpdir = 'tmp.cluster_test_full_run_ok_variants_only.not_present'
        common.rmtree(tmpdir)
        shutil.copytree(os.path.join(data_dir, 'cluster_test_full_run_ok_variants_only'), tmpdir)

        c = cluster.Cluster(tmpdir, 'cluster_name', refdata, total_reads=66, total_reads_bases=3300)
        c.run()
        expected = [
            'variants_only1\tvariants_only1\t1\t1\t27\t66\tcluster_name\t96\t96\t100.0\tcluster_name.l15.c30.ctg.1\t215\t15.3\t1\tSNP\tp\tR3S\t0\t.\t.\t7\t9\tCGC\t65\t67\tCGC\t18;18;19\tC;G;C\t18;18;19\tvariants_only1:1:1:R3S:.:Ref and assembly have wild type, so do not report\tGeneric description of variants_only1'
        ]
        self.assertEqual(expected, c.report_lines)
        common.rmtree(tmpdir)


    def test_full_run_ok_variants_only_variant_not_present_always_report(self):
        '''test complete run of cluster on a variants only gene when variant not present but always report variant'''
        fasta_in = os.path.join(data_dir, 'cluster_test_full_run_ok_variants_only.fa')
        tsv_in = os.path.join(data_dir, 'cluster_full_run_varonly.not_present.always_report.tsv')
        refdata = reference_data.ReferenceData([fasta_in], [tsv_in])
        tmpdir = 'tmp.cluster_full_run_varonly.not_present.always_report'
        common.rmtree(tmpdir)
        shutil.copytree(os.path.join(data_dir, 'cluster_test_full_run_ok_variants_only'), tmpdir)

        c = cluster.Cluster(tmpdir, 'cluster_name', refdata, total_reads=66, total_reads_bases=3300)
        c.run()
        expected = [
            'variants_only1\tvariants_only1\t1\t1\t27\t66\tcluster_name\t96\t96\t100.0\tcluster_name.l15.c30.ctg.1\t215\t15.3\t1\tSNP\tp\tR3S\t0\t.\t.\t7\t9\tCGC\t65\t67\tCGC\t18;18;19\tC;G;C\t18;18;19\tvariants_only1:1:1:R3S:.:Ref and assembly have wild type, but always report anyway\tGeneric description of variants_only1'
        ]
        self.assertEqual(expected, c.report_lines)
        common.rmtree(tmpdir)


    def test_full_run_ok_variants_only_variant_is_present(self):
        '''test complete run of cluster on a variants only gene when variant is present'''
        fasta_in = os.path.join(data_dir, 'cluster_test_full_run_ok_variants_only.fa')
        tsv_in = os.path.join(data_dir, 'cluster_test_full_run_ok_variants_only.present.metadata.tsv')
        refdata = reference_data.ReferenceData([fasta_in], [tsv_in])
        tmpdir = 'tmp.cluster_test_full_run_ok_variants_only.present'
        common.rmtree(tmpdir)
        shutil.copytree(os.path.join(data_dir, 'cluster_test_full_run_ok_variants_only'), tmpdir)

        c = cluster.Cluster(tmpdir, 'cluster_name', refdata, total_reads=66, total_reads_bases=3300)
        c.run()

        expected = [
            'variants_only1\tvariants_only1\t1\t1\t27\t66\tcluster_name\t96\t96\t100.0\tcluster_name.l15.c30.ctg.1\t215\t15.3\t1\tSNP\tp\tR3S\t0\t.\t.\t7\t9\tCGC\t65\t67\tCGC\t18;18;19\tC;G;C\t18;18;19\tvariants_only1:1:1:R3S:.:Ref and assembly have wild type\tGeneric description of variants_only1',
            'variants_only1\tvariants_only1\t1\t1\t27\t66\tcluster_name\t96\t96\t100.0\tcluster_name.l15.c30.ctg.1\t215\t15.3\t1\tSNP\tp\tI5A\t1\t.\t.\t13\t15\tGCG\t71\t73\tGCG\t17;17;17\tG;C;G\t17;17;17\tvariants_only1:1:1:I5A:.:Ref and reads have variant so report\tGeneric description of variants_only1',
        ]
        self.assertEqual(expected, c.report_lines)
        common.rmtree(tmpdir)


    def test_full_run_ok_gene_start_mismatch(self):
        '''test complete run where gene extended because too different at end for full nucmer match'''
        fasta_in = os.path.join(data_dir, 'cluster_test_full_run_ok_gene_start_mismatch.fa')
        tsv_in = os.path.join(data_dir, 'cluster_test_full_run_ok_gene_start_mismatch.metadata.tsv')
        refdata = reference_data.ReferenceData([fasta_in], [tsv_in])
        tmpdir = 'tmp.cluster_test_full_run_ok_gene_start_mismatch'
        common.rmtree(tmpdir)
        shutil.copytree(os.path.join(data_dir, 'cluster_test_full_run_ok_gene_start_mismatch'), tmpdir)
        c = cluster.Cluster(tmpdir, 'cluster_name', refdata, total_reads=112, total_reads_bases=1080)
        c.run()
        expected = [
            'gene\tgene\t1\t0\t27\t112\tcluster_name\t96\t96\t100.0\tcluster_name.l6.c30.ctg.1\t362\t27.8\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\tGeneric description of gene'
        ]
        self.assertEqual(expected, c.report_lines)
        common.rmtree(tmpdir)


    def test_full_run_smtls_snp_presabs_gene(self):
        '''test complete run where samtools calls a snp in a presence/absence gene'''
        fasta_in = os.path.join(data_dir, 'cluster_full_run_smtls_snp_presabs_gene.fa')
        tsv_in = os.path.join(data_dir, 'cluster_full_run_smtls_snp_presabs_gene.tsv')
        refdata = reference_data.ReferenceData([fasta_in], [tsv_in])
        tmpdir = 'tmp.cluster_test_full_run_ok_samtools_snp_pres_abs_gene'
        common.rmtree(tmpdir)
        shutil.copytree(os.path.join(data_dir, 'cluster_full_run_smtls_snp_presabs_gene'), tmpdir)
        c = cluster.Cluster(tmpdir, 'cluster_name', refdata, total_reads=148, total_reads_bases=13320)
        c.run()
        expected = [
            'ref_gene\tref_gene\t1\t0\t155\t148\tcluster_name\t96\t96\t100.0\tcluster_name.l15.c30.ctg.1\t335\t39.8\t0\tHET\t.\t.\t.\tG18A\t.\t18\t18\tG\t137\t137\tG\t63\tG,A\t32,31\t.\tGeneric description of ref_gene'
        ]
        self.assertEqual(expected, c.report_lines)
        common.rmtree(tmpdir)


    def test_full_run_smtls_snp_varonly_gene_2(self):
        '''test complete run where samtools calls a snp in a variant only gene'''
        # _2 because I think test_full_run_smtls_snp_varonly_gene tests the asame functionality.
        # ... but let's leave both tests in anyway
        fasta_in = os.path.join(data_dir, 'cluster_full_run_smtls_snp_varonly_gene_2.fa')
        tsv_in = os.path.join(data_dir, 'cluster_full_run_smtls_snp_varonly_gene_2.tsv')
        refdata = reference_data.ReferenceData([fasta_in], [tsv_in])
        tmpdir = 'tmp.cluster_full_run_smtls_snp_varonly_gene_2'
        common.rmtree(tmpdir)
        shutil.copytree(os.path.join(data_dir, 'cluster_full_run_smtls_snp_varonly_gene_2'), tmpdir)
        c = cluster.Cluster(tmpdir, 'cluster_name', refdata, total_reads=148, total_reads_bases=13320)
        c.run()
        expected = [
            'ref_gene\tref_gene\t1\t1\t155\t148\tcluster_name\t96\t96\t100.0\tcluster_name.l15.c30.ctg.1\t335\t39.8\t0\tHET\t.\t.\t.\tG18A\t.\t18\t18\tG\t137\t137\tG\t63\tG,A\t32,31\t.\tGeneric description of ref_gene'
        ]
        self.assertEqual(expected, c.report_lines)
        common.rmtree(tmpdir)


    def test_full_run_known_smtls_snp_presabs_gene(self):
        '''test complete run where samtools calls a snp at a known snp location in a presence/absence gene'''
        fasta_in = os.path.join(data_dir, 'cluster_full_run_known_smtls_snp_presabs_gene.fa')
        tsv_in = os.path.join(data_dir, 'cluster_full_run_known_smtls_snp_presabs_gene.tsv')
        refdata = reference_data.ReferenceData([fasta_in], [tsv_in])
        tmpdir = 'tmp.cluster_test_full_run_ok_samtools_snp_known_position_pres_abs_gene'
        common.rmtree(tmpdir)
        shutil.copytree(os.path.join(data_dir, 'cluster_full_run_known_smtls_snp_presabs_gene'), tmpdir)
        c = cluster.Cluster(tmpdir, 'cluster_name', refdata, total_reads=148, total_reads_bases=13320)
        c.run()

        # We shouldn't get an extra 'HET' line because we already know about the snp, so
        # included in the report of the known snp
        expected = [
            'ref_gene\tref_gene\t1\t0\t155\t148\tcluster_name\t96\t96\t100.0\tcluster_name.l15.c30.ctg.1\t335\t39.8\t1\tSNP\tp\tM6I\t0\t.\t.\t16\t18\tATG\t135\t137\tATG\t65;64;63\tA;T;G,A\t65;64;32,31\tref_gene:1:0:M6I:.:Description of M6I snp\t.'
        ]
        self.assertEqual(expected, c.report_lines)
        common.rmtree(tmpdir)


    def test_full_run_smtls_snp_varonly_gene_no_snp(self):
        '''test complete run where samtools calls a snp at a known snp location in a variant only gene, gene does not have variant'''
        fasta_in = os.path.join(data_dir, 'cluster_full_run_smtls_snp_varonly_gene_no_snp.fa')
        tsv_in = os.path.join(data_dir, 'cluster_full_run_smtls_snp_varonly_gene_no_snp.tsv')
        refdata = reference_data.ReferenceData([fasta_in], [tsv_in])
        tmpdir = 'tmp.cluster_test_full_run_smtls_snp_varonly_gene_no_snp'
        common.rmtree(tmpdir)
        shutil.copytree(os.path.join(data_dir, 'cluster_full_run_smtls_snp_varonly_gene_no_snp'), tmpdir)
        c = cluster.Cluster(tmpdir, 'cluster_name', refdata, total_reads=148, total_reads_bases=13320)
        c.run()

        # We shouldn't get an extra 'HET' line because we already know about the snp, so
        # included in the report of the known snp
        expected = [
            'ref_gene\tref_gene\t1\t1\t155\t148\tcluster_name\t96\t96\t100.0\tcluster_name.l15.c30.ctg.1\t335\t39.8\t1\tSNP\tp\tM6I\t0\t.\t.\t16\t18\tATG\t135\t137\tATG\t65;64;63\tA;T;G,A\t65;64;32,31\tref_gene:1:1:M6I:.:Description of M6I snp\t.'
        ]
        self.assertEqual(expected, c.report_lines)
        common.rmtree(tmpdir)


    def test_full_run_smtls_snp_varonly_gene(self):
        '''test complete run where samtools calls a snp at a known snp location in a variant only gene, gene does have variant'''
        fasta_in = os.path.join(data_dir, 'cluster_full_run_smtls_snp_varonly_gene.fa')
        tsv_in = os.path.join(data_dir, 'cluster_full_run_smtls_snp_varonly_gene.tsv')
        refdata = reference_data.ReferenceData([fasta_in], [tsv_in])
        tmpdir = 'tmp.cluster_test_full_run_ok_samtools_snp_known_position_var_only_gene_does_have_var'
        common.rmtree(tmpdir)
        shutil.copytree(os.path.join(data_dir, 'cluster_full_run_smtls_snp_varonly_gene'), tmpdir)
        c = cluster.Cluster(tmpdir, 'cluster_name', refdata, total_reads=148, total_reads_bases=13320)
        c.run()

        # We shouldn't get an extra 'HET' line because we already know about the snp, so
        # included in the report of the known snp
        expected = [
            'ref_gene\tref_gene\t1\t1\t155\t148\tcluster_name\t96\t96\t100.0\tcluster_name.l15.c30.ctg.1\t335\t39.8\t1\tSNP\tp\tI6M\t1\t.\t.\t16\t18\tATG\t135\t137\tATG\t65;64;63\tA;T;G,A\t65;64;32,31\tref_gene:1:1:I6M:.:Description of I6M snp\t.'
        ]
        self.assertEqual(expected, c.report_lines)
        common.rmtree(tmpdir)


    def test_full_run_smtls_snp_presabs_nonc(self):
        '''test complete run where samtools calls a snp in a presence/absence noncoding sequence'''
        fasta_in = os.path.join(data_dir, 'cluster_full_run_smtls_snp_presabs_nonc.fa')
        tsv_in = os.path.join(data_dir, 'cluster_full_run_smtls_snp_presabs_nonc.tsv')
        refdata = reference_data.ReferenceData([fasta_in], [tsv_in])
        tmpdir = 'tmp.cluster_test_full_run_smtls_snp_presabs_nonc'
        common.rmtree(tmpdir)
        shutil.copytree(os.path.join(data_dir, 'cluster_full_run_smtls_snp_presabs_nonc'), tmpdir)
        c = cluster.Cluster(tmpdir, 'cluster_name', refdata, total_reads=148, total_reads_bases=13320)
        c.run()
        expected = [
            'ref_seq\tref_seq\t0\t0\t147\t148\tcluster_name\t96\t96\t100.0\tcluster_name.l15.c30.ctg.1\t335\t39.8\t0\tHET\t.\t.\t.\tG18A\t.\t18\t18\tG\t137\t137\tG\t63\tG,A\t32,31\t.\tGeneric description of ref_seq'
        ]
        self.assertEqual(expected, c.report_lines)
        common.rmtree(tmpdir)


    def test_full_run_smtls_known_snp_presabs_nonc(self):
        '''test complete run where samtools calls a snp in a presence/absence noncoding sequence at a known snp position'''
        fasta_in = os.path.join(data_dir, 'cluster_full_run_smtls_known_snp_presabs_nonc.fa')
        tsv_in = os.path.join(data_dir, 'cluster_full_run_smtls_known_snp_presabs_nonc.tsv')
        refdata = reference_data.ReferenceData([fasta_in], [tsv_in])
        tmpdir = 'tmp.cluster_test_full_run_smtls_known_snp_presabs_nonc'
        common.rmtree(tmpdir)
        shutil.copytree(os.path.join(data_dir, 'cluster_full_run_smtls_known_snp_presabs_nonc'), tmpdir)
        c = cluster.Cluster(tmpdir, 'cluster_name', refdata, total_reads=148, total_reads_bases=13320)
        c.run()
        expected = [
            'ref_seq\tref_seq\t0\t0\t147\t148\tcluster_name\t96\t96\t100.0\tcluster_name.l15.c30.ctg.1\t335\t39.8\t1\tSNP\tn\tG18A\t0\t.\t.\t18\t18\tG\t137\t137\tG\t63\tG,A\t32,31\tref_seq:0:0:G18A:.:Description of G18A\tGeneric description of ref_seq'
        ]
        self.assertEqual(expected, c.report_lines)
        common.rmtree(tmpdir)


    def test_full_run_smtls_snp_varonly_nonc(self):
        '''test complete run where samtools calls a snp in a presence/absence noncoding sequence'''
        fasta_in = os.path.join(data_dir, 'cluster_full_run_smtls_snp_varonly_nonc.fa')
        tsv_in = os.path.join(data_dir, 'cluster_full_run_smtls_snp_varonly_nonc.tsv')
        refdata = reference_data.ReferenceData([fasta_in], [tsv_in])
        tmpdir = 'tmp.cluster_full_run_smtls_snp_varonly_nonc'
        common.rmtree(tmpdir)
        shutil.copytree(os.path.join(data_dir, 'cluster_full_run_smtls_snp_varonly_nonc'), tmpdir)
        c = cluster.Cluster(tmpdir, 'cluster_name', refdata, total_reads=148, total_reads_bases=13320)
        c.run()
        expected = [
            'ref_seq\tref_seq\t0\t1\t147\t148\tcluster_name\t96\t96\t100.0\tcluster_name.l15.c30.ctg.1\t335\t39.8\t0\tHET\t.\t.\t.\tG18A\t.\t18\t18\tG\t137\t137\tG\t63\tG,A\t32,31\t.\tGeneric description of ref_seq'
        ]
        self.assertEqual(expected, c.report_lines)
        common.rmtree(tmpdir)


    def test_full_run_known_smtls_snp_presabs_nonc(self):
        '''test complete run where samtools calls a snp at a known snp location in a presence/absence noncoding'''
        fasta_in = os.path.join(data_dir, 'cluster_full_run_known_smtls_snp_presabs_nonc.fa')
        tsv_in = os.path.join(data_dir, 'cluster_full_run_known_smtls_snp_presabs_nonc.tsv')
        refdata = reference_data.ReferenceData([fasta_in], [tsv_in])
        tmpdir = 'tmp.cluster_test_full_run_ok_samtools_snp_known_position_pres_abs_noncoding'
        common.rmtree(tmpdir)
        shutil.copytree(os.path.join(data_dir, 'cluster_full_run_known_smtls_snp_presabs_nonc'), tmpdir)
        c = cluster.Cluster(tmpdir, 'cluster_name', refdata, total_reads=148, total_reads_bases=13320)
        c.run()

        # We shouldn't get an extra 'HET' line because we already know about the snp, so
        # included in the report of the known snp
        expected = [
            'ref_seq\tref_seq\t0\t0\t147\t148\tcluster_name\t96\t96\t100.0\tcluster_name.l15.c30.ctg.1\t335\t39.8\t1\tSNP\tn\tG18A\t0\t.\t.\t18\t18\tG\t137\t137\tG\t63\tG,A\t32,31\tref_seq:0:0:G18A:.:Description of G18A snp\t.'
        ]
        self.assertEqual(expected, c.report_lines)
        common.rmtree(tmpdir)


    def test_full_run_smtls_snp_varonly_nonc_no_snp(self):
        '''test complete run where samtools calls a snp at a known snp location in a presence/absence noncoding and sample does not have the var'''
        fasta_in = os.path.join(data_dir, 'cluster_full_run_smtls_snp_varonly_nonc_no_snp.fa')
        tsv_in = os.path.join(data_dir, 'cluster_full_run_smtls_snp_varonly_nonc_no_snp.tsv')
        refdata = reference_data.ReferenceData([fasta_in], [tsv_in])
        tmpdir = 'tmp.cluster_test_full_run_ok_samtools_snp_known_position_var_only_noncoding'
        common.rmtree(tmpdir)
        shutil.copytree(os.path.join(data_dir, 'cluster_full_run_smtls_snp_varonly_nonc_no_snp'), tmpdir)
        c = cluster.Cluster(tmpdir, 'cluster_name', refdata, total_reads=148, total_reads_bases=13320)
        c.run()

        # We shouldn't get an extra 'HET' line because we already know about the snp, so
        # included in the report of the known snp
        expected = [
            'ref_seq\tref_seq\t0\t1\t147\t148\tcluster_name\t96\t96\t100.0\tcluster_name.l15.c30.ctg.1\t335\t39.8\t1\tSNP\tn\tG18A\t0\t.\t.\t18\t18\tG\t137\t137\tG\t63\tG,A\t32,31\tref_seq:0:1:G18A:.:Description of G18A snp\t.'
        ]
        self.assertEqual(expected, c.report_lines)
        common.rmtree(tmpdir)


    def test_full_run_cluster_test_full_run_smtls_snp_varonly_nonc(self):
        '''test complete run where samtools calls a snp at a known snp location in a presence/absence noncoding and sample has the var'''
        fasta_in = os.path.join(data_dir, 'cluster_test_full_run_smtls_snp_varonly_nonc.fa')
        tsv_in = os.path.join(data_dir, 'cluster_test_full_run_smtls_snp_varonly_nonc.tsv')
        refdata = reference_data.ReferenceData([fasta_in], [tsv_in])
        tmpdir = 'tmp.cluster_test_full_run_ok_samtools_snp_known_position_var_only_noncoding'
        common.rmtree(tmpdir)
        shutil.copytree(os.path.join(data_dir, 'cluster_test_full_run_smtls_snp_varonly_nonc'), tmpdir)
        c = cluster.Cluster(tmpdir, 'cluster_name', refdata, total_reads=148, total_reads_bases=13320)
        c.run()

        # We shouldn't get an extra 'HET' line because we already know about the snp, so
        # included in the report of the known snp
        expected = [
            'ref_seq\tref_seq\t0\t1\t147\t148\tcluster_name\t96\t96\t100.0\tcluster_name.l15.c30.ctg.1\t335\t39.8\t1\tSNP\tn\tA18G\t1\t.\t.\t18\t18\tG\t137\t137\tG\t63\tG,A\t32,31\tref_seq:0:1:A18G:.:Description of A18G snp\t.'
        ]
        self.assertEqual(expected, c.report_lines)
        common.rmtree(tmpdir)


    def test_full_run_partial_assembly(self):
        '''Test complete run where only part of the ref gene is present in the reads'''
        fasta_in = os.path.join(data_dir, 'cluster_test_full_run_partial_asmbly.fa')
        tsv_in = os.path.join(data_dir, 'cluster_test_full_run_partial_asmbly.tsv')
        refdata = reference_data.ReferenceData([fasta_in], [tsv_in])
        tmpdir = 'tmp.cluster_test_full_run_partial_assembly'
        common.rmtree(tmpdir)
        shutil.copytree(os.path.join(data_dir, 'cluster_test_full_run_partial_asmbly'), tmpdir)
        c = cluster.Cluster(tmpdir, 'cluster_name', refdata, total_reads=278, total_reads_bases=15020)
        c.run()

        expected = [
            'presence_absence1\tpresence_absence1\t1\t0\t19\t278\tcluster_name\t96\t77\t100.0\tcluster_name.l15.c17.ctg.1\t949\t20.5\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\tGeneric description of presence_absence1'
        ]
        self.assertEqual(expected, c.report_lines)
        common.rmtree(tmpdir)


    def test_full_run_multiple_vars_in_codon(self):
        '''Test complete run where there is a codon with a SNP and an indel'''
        fasta_in = os.path.join(data_dir, 'cluster_test_full_run_multiple_vars.fa')
        tsv_in = os.path.join(data_dir, 'cluster_test_full_run_multiple_vars.tsv')
        refdata = reference_data.ReferenceData([fasta_in], [tsv_in])
        tmpdir = 'tmp.cluster_test_full_run_multiple_vars'
        common.rmtree(tmpdir)
        shutil.copytree(os.path.join(data_dir, 'cluster_test_full_run_multiple_vars'), tmpdir)
        c = cluster.Cluster(tmpdir, 'cluster_name', refdata, total_reads=292, total_reads_bases=20900)
        c.run()

        expected = [
            'presence_absence1\tpresence_absence1\t1\t0\t539\t292\tcluster_name\t96\t96\t96.91\tcluster_name.l15.c30.ctg.1\t1074\t20.4\t0\t.\tp\t.\t0\t.\tMULTIPLE\t25\t26\tGA\t487\t489\tCAT\t27;26;25\tC;A;T\t27;26;25\t.\tGeneric description of presence_absence1',
            'presence_absence1\tpresence_absence1\t1\t0\t539\t292\tcluster_name\t96\t96\t96.91\tcluster_name.l15.c30.ctg.1\t1074\t20.4\t0\t.\tp\t.\t0\tA10fs\tFSHIFT\t28\t28\tG\t491\t491\tG\t26\tG\t26\t.\tGeneric description of presence_absence1',
        ]
        self.assertEqual(expected, c.report_lines)
        common.rmtree(tmpdir)


    def test_full_run_delete_codon(self):
        '''Test complete run where there is a deleted codon'''
        fasta_in = os.path.join(data_dir, 'cluster_test_full_run_delete_codon.fa')
        tsv_in = os.path.join(data_dir, 'cluster_test_full_run_delete_codon.tsv')
        refdata = reference_data.ReferenceData([fasta_in], [tsv_in])
        tmpdir = 'tmp.cluster_test_full_delete_codon'
        common.rmtree(tmpdir)
        shutil.copytree(os.path.join(data_dir, 'cluster_test_full_run_delete_codon'), tmpdir)
        c = cluster.Cluster(tmpdir, 'cluster_name', refdata, total_reads=292, total_reads_bases=20900)
        c.run()

        expected = [
            'presence_absence1\tpresence_absence1\t1\t0\t539\t292\tcluster_name\t117\t117\t92.31\tcluster_name.l15.c30.ctg.1\t1104\t20.0\t0\t.\tp\t.\t0\tR25_A26del\tDEL\t73\t73\tA\t553\t553\tA\t27\tA\t27\t.\tGeneric description of presence_absence1',
        ]
        self.assertEqual(expected, c.report_lines)
        common.rmtree(tmpdir)


    def test_full_run_insert_codon(self):
        '''Test complete run where there is a inserted codon'''
        fasta_in = os.path.join(data_dir, 'cluster_test_full_run_insert_codon.fa')
        tsv_in = os.path.join(data_dir, 'cluster_test_full_run_insert_codon.tsv')
        refdata = reference_data.ReferenceData([fasta_in], [tsv_in])
        tmpdir = 'tmp.cluster_test_full_insert_codon'
        common.rmtree(tmpdir)
        shutil.copytree(os.path.join(data_dir, 'cluster_test_full_run_insert_codon'), tmpdir)
        c = cluster.Cluster(tmpdir, 'cluster_name', refdata, total_reads=292, total_reads_bases=20900)
        c.run()

        expected = [
            'presence_absence1\tpresence_absence1\t1\t0\t539\t292\tcluster_name\t108\t108\t92.31\tcluster_name.l15.c30.ctg.1\t1115\t19.9\t0\t.\tp\t.\t0\tS25_M26insELI\tINS\t73\t73\tA\t554\t554\tG\t24\tG\t24\t.\tGeneric description of presence_absence1'
        ]
        self.assertEqual(expected, c.report_lines)
        common.rmtree(tmpdir)
