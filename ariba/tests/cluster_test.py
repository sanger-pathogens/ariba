import unittest
import os
import copy
import shutil
import filecmp
import pyfastaq
import pysam
import pymummer
from ariba import cluster, flag, reference_data

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
                shutil.rmtree(full_path)
            else:
                os.unlink(full_path)


class TestCluster(unittest.TestCase):
    def test_init_fail_files_missing(self):
        '''test init_fail_files_missing'''
        refdata_fa = os.path.join(data_dir, 'cluster_test_init_refdata.fa')
        refdata = reference_data.ReferenceData(presence_absence_fa=refdata_fa)

        dirs = [
            'cluster_test_init_no_refs_fa',
            'cluster_test_init_no_reads_1',
            'cluster_test_init_no_reads_2',
        ]
        dirs = [os.path.join(data_dir, d) for d in dirs]
        for d in dirs:
            tmpdir = 'tmp.cluster_test_init_fail_files_missing'
            shutil.copytree(d, tmpdir)
            with self.assertRaises(cluster.Error):
                c = cluster.Cluster(tmpdir, 'name', refdata=refdata)
            shutil.rmtree(tmpdir)

        with self.assertRaises(cluster.Error):
            c = cluster.Cluster('directorydoesnotexistshouldthrowerror', 'name', refdata=refdata)


    def test_count_reads(self):
        '''test _count_reads pass'''
        reads1 = os.path.join(data_dir, 'cluster_test_count_reads_1.fq')
        reads2 = os.path.join(data_dir, 'cluster_test_count_reads_2.fq')
        self.assertEqual(4, cluster.Cluster._count_reads(reads1, reads2))


    def test_full_run_choose_ref_fail(self):
        '''test complete run of cluster when choosing ref seq fails'''
        refdata = reference_data.ReferenceData(
            presence_absence_fa=os.path.join(data_dir, 'cluster_test_full_run_choose_ref_fail.presence_absence.fa')
        )
        tmpdir = 'tmp.test_full_run_choose_ref_fail'
        shutil.copytree(os.path.join(data_dir, 'cluster_test_full_run_choose_ref_fail'), tmpdir)

        c = cluster.Cluster(tmpdir, 'cluster_name', refdata)
        c.run()

        expected = '\t'.join(['.', '.', '1088', '2', 'cluster_name'] + ['.'] * 22)
        self.assertEqual([expected], c.report_lines)
        self.assertTrue(c.status_flag.has('ref_seq_choose_fail'))
        self.assertTrue(c.status_flag.has('assembly_fail'))
        shutil.rmtree(tmpdir)


    def test_full_run_assembly_fail(self):
        '''test complete run of cluster when assembly fails'''
        refdata = reference_data.ReferenceData(
            non_coding_fa=os.path.join(data_dir, 'cluster_test_full_run_assembly_fail.noncoding.fa')
        )
        tmpdir = 'tmp.test_full_run_assembly_fail'
        shutil.copytree(os.path.join(data_dir, 'cluster_test_full_run_assembly_fail'), tmpdir)

        c = cluster.Cluster(tmpdir, 'cluster_name', refdata)
        c.run()

        expected = '\t'.join(['noncoding_ref_seq', 'non_coding', '64', '4', 'cluster_name'] + ['.'] * 22)
        self.assertEqual([expected], c.report_lines)
        self.assertFalse(c.status_flag.has('ref_seq_choose_fail'))
        self.assertTrue(c.status_flag.has('assembly_fail'))
        shutil.rmtree(tmpdir)


    def test_full_run_ok_non_coding(self):
        '''test complete run of cluster on a noncoding sequence'''
        refdata = reference_data.ReferenceData(
            non_coding_fa=os.path.join(data_dir, 'cluster_test_full_run_ok_non_coding.fa'),
            metadata_tsv=os.path.join(data_dir, 'cluster_test_full_run_ok_non_coding.metadata.tsv')
        )

        tmpdir = 'tmp.test_full_run_ok_non_coding'
        shutil.copytree(os.path.join(data_dir, 'cluster_test_full_run_ok_non_coding'), tmpdir)

        c = cluster.Cluster(tmpdir, 'cluster_name', refdata, spades_other_options='--only-assembler')
        c.run()

        expected = [
            'noncoding1\tnon_coding\t19\t72\tcluster_name\t120\t120\t95.87\tnoncoding1.scaffold.1\t234\tT\tSNP\tn\tT\tA14T\tSNP\t13\t13\tA\t73\t73\tT\t19\t.\t19\tnoncoding1_n_A14T_N_ref has wild type, reads has variant so should report\tgeneric description of noncoding1',
            'noncoding1\tnon_coding\t19\t72\tcluster_name\t120\t120\t95.87\tnoncoding1.scaffold.1\t234\tF\t.\tn\tF\tG61T\tSNP\t60\t60\tG\t120\t120\tT\t24\t.\t24\t.\tgeneric description of noncoding1',
            'noncoding1\tnon_coding\t19\t72\tcluster_name\t120\t120\t95.87\tnoncoding1.scaffold.1\t234\tF\t.\tn\tF\t.82C\tINS\t81\t81\t.\t142\t142\tC\t23\t.\t23\t.\tgeneric description of noncoding1',
            'noncoding1\tnon_coding\t19\t72\tcluster_name\t120\t120\t95.87\tnoncoding1.scaffold.1\t234\tF\t.\tn\tF\tT108.\tDEL\t107\t107\tT\t167\t167\t.\t17\t.\t17\t.\tgeneric description of noncoding1',
            'noncoding1\tnon_coding\t19\t72\tcluster_name\t120\t120\t95.87\tnoncoding1.scaffold.1\t234\tT\tSNP\tn\tT\t.\t.\t6\t6\tG\t66\t66\tG\t19\t.\t19\tnoncoding1_n_A6G_N_variant in ref and reads so should report\tgeneric description of noncoding1',
        ]

        self.assertEqual(expected, c.report_lines)
        shutil.rmtree(tmpdir)


    def test_full_run_ok_presence_absence(self):
        '''test complete run of cluster on a presence absence gene'''
        refdata = reference_data.ReferenceData(
            presence_absence_fa=os.path.join(data_dir, 'cluster_test_full_run_ok_presence_absence.fa'),
            metadata_tsv=os.path.join(data_dir, 'cluster_test_full_run_ok_presence_absence.metadata.tsv'),
        )

        tmpdir = 'tmp.cluster_test_full_run_ok_presence_absence'
        shutil.copytree(os.path.join(data_dir, 'cluster_test_full_run_ok_presence_absence'), tmpdir)

        c = cluster.Cluster(tmpdir, 'cluster_name', refdata, spades_other_options='--only-assembler')
        c.run()

        expected = [
            'presence_absence1\tpresence_absence\t27\t64\tcluster_name\t96\t96\t97.92\tpresence_absence1.scaffold.1\t213\tT\tSNP\tp\tT\tA10V\tNONSYN\t28\t28\tC\t83\t83\tT\t22\t.\t22\tpresence_absence1_p_A10V_N_Ref has wild, reads have variant so report\tGeneric description of presence_absence1',
            'presence_absence1\tpresence_absence\t27\t64\tcluster_name\t96\t96\t97.92\tpresence_absence1.scaffold.1\t213\tF\t.\tp\tF\t.\tSYN\t53\t53\tT\t108\t108\tC\t32\t.\t32\t.\tGeneric description of presence_absence1',
            'presence_absence1\tpresence_absence\t27\t64\tcluster_name\t96\t96\t97.92\tpresence_absence1.scaffold.1\t213\tT\tSNP\tp\tT\t.\t.\t13\t15\tG;C;G\t68\t70\tG;C;G\t18;20;20\t.;.;.\t18;20;20\tpresence_absence1_p_I5A_N_Ref and reads have variant so report\tGeneric description of presence_absence1',
        ]

        self.assertEqual(expected, c.report_lines)
        shutil.rmtree(tmpdir)


    def test_full_run_ok_variants_only_gene_not_present(self):
        '''test complete run of cluster on a variants only gene when variant not present'''
        refdata = reference_data.ReferenceData(
            variants_only_fa=os.path.join(data_dir, 'cluster_test_full_run_ok_variants_only.fa'),
            metadata_tsv=os.path.join(data_dir, 'cluster_test_full_run_ok_variants_only.not_present.metadata.tsv'),
        )

        tmpdir = 'tmp.cluster_test_full_run_ok_variants_only.not_present'
        shutil.copytree(os.path.join(data_dir, 'cluster_test_full_run_ok_variants_only'), tmpdir)

        c = cluster.Cluster(tmpdir, 'cluster_name', refdata, spades_other_options='--only-assembler')
        c.run()
        expected = [
            'variants_only1\tvariants_only\t27\t66\tcluster_name\t96\t96\t100.0\tvariants_only1.scaffold.1\t215\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\t.\tGeneric description of variants_only1'
        ]
        self.assertEqual(expected, c.report_lines)
        shutil.rmtree(tmpdir)


    def test_full_run_ok_variants_only_gene_is_present(self):
        '''test complete run of cluster on a variants only gene when variant is present'''
        refdata = reference_data.ReferenceData(
            variants_only_fa=os.path.join(data_dir, 'cluster_test_full_run_ok_variants_only.fa'),
            metadata_tsv=os.path.join(data_dir, 'cluster_test_full_run_ok_variants_only.present.metadata.tsv'),
        )

        tmpdir = 'tmp.cluster_test_full_run_ok_variants_only.present'
        shutil.copytree(os.path.join(data_dir, 'cluster_test_full_run_ok_variants_only'), tmpdir)

        c = cluster.Cluster(tmpdir, 'cluster_name', refdata, spades_other_options='--only-assembler')
        c.run()

        expected = [
            'variants_only1\tvariants_only\t27\t66\tcluster_name\t96\t96\t100.0\tvariants_only1.scaffold.1\t215\tT\tSNP\tp\tT\t.\t.\t13\t15\tG;C;G\t71\t73\tG;C;G\t17;17;17\t.;.;.\t17;17;17\tvariants_only1_p_I5A_N_Ref and reads have variant so report\tGeneric description of variants_only1'
        ]

        self.assertEqual(expected, c.report_lines)
        shutil.rmtree(tmpdir)
