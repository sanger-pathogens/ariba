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
                c = cluster.Cluster(tmpdir, 'name', refdata=refdata, total_reads=42, total_reads_bases=4242)
            shutil.rmtree(tmpdir)


    def test_number_of_reads_for_assembly(self):
        '''Test _number_of_reads_for_assembly'''
        # ref is 100bp long
        ref_fa = os.path.join(data_dir, 'cluster_test_number_of_reads_for_assembly.ref.fa')
        tests = [
            (50, 1000, 10, 20, 40),
            (50, 999, 10, 20, 42),
            (50, 1000, 10, 10, 20),
            (50, 1000, 10, 5, 10),
        ]

        for insert, bases, reads, coverage, expected in tests:
            self.assertEqual(expected, cluster.Cluster._number_of_reads_for_assembly(ref_fa, insert, bases, reads, coverage))


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
        expected_out1 = os.path.join(data_dir, 'cluster_test_make_reads_for_assembly.out1.fq')
        expected_out2 = os.path.join(data_dir, 'cluster_test_make_reads_for_assembly.out2.fq')
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


    def test_full_run_choose_ref_fail(self):
        '''test complete run of cluster when choosing ref seq fails'''
        refdata = reference_data.ReferenceData(
            presence_absence_fa=os.path.join(data_dir, 'cluster_test_full_run_choose_ref_fail.presence_absence.fa')
        )
        tmpdir = 'tmp.test_full_run_choose_ref_fail'
        shutil.copytree(os.path.join(data_dir, 'cluster_test_full_run_choose_ref_fail'), tmpdir)

        c = cluster.Cluster(tmpdir, 'cluster_name', refdata, total_reads=2, total_reads_bases=108)
        c.run()

        expected = '\t'.join(['.', '.', '1088', '2', 'cluster_name'] + ['.'] * 24)
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

        c = cluster.Cluster(tmpdir, 'cluster_name', refdata, total_reads=4, total_reads_bases=304)
        c.run()

        expected = '\t'.join(['noncoding_ref_seq', 'non_coding', '64', '4', 'cluster_name'] + ['.'] * 24)
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

        c = cluster.Cluster(tmpdir, 'cluster_name', refdata, spades_other_options='--only-assembler', total_reads=72, total_reads_bases=3600)
        c.run()

        expected = [
            'noncoding1\tnon_coding\t531\t72\tcluster_name\t120\t120\t95.87\tnoncoding1.scaffold.1\t234\t15.4\t1\tSNP\tn\tA14T\t1\tA14T\tSNP\t13\t13\tA\t73\t73\tT\t19\t.\t19\tnoncoding1_n_A14T_ref has wild type, reads has variant so should report\tgeneric description of noncoding1',
            'noncoding1\tnon_coding\t531\t72\tcluster_name\t120\t120\t95.87\tnoncoding1.scaffold.1\t234\t15.4\t0\t.\tn\t.\t0\tG61T\tSNP\t60\t60\tG\t120\t120\tT\t24\t.\t24\t.\tgeneric description of noncoding1',
            'noncoding1\tnon_coding\t531\t72\tcluster_name\t120\t120\t95.87\tnoncoding1.scaffold.1\t234\t15.4\t0\t.\tn\t.\t0\t.82C\tINS\t81\t81\t.\t142\t142\tC\t23\t.\t23\t.\tgeneric description of noncoding1',
            'noncoding1\tnon_coding\t531\t72\tcluster_name\t120\t120\t95.87\tnoncoding1.scaffold.1\t234\t15.4\t0\t.\tn\t.\t0\tT108.\tDEL\t107\t107\tT\t167\t167\t.\t17\t.\t17\t.\tgeneric description of noncoding1',
            'noncoding1\tnon_coding\t531\t72\tcluster_name\t120\t120\t95.87\tnoncoding1.scaffold.1\t234\t15.4\t1\tSNP\tn\tA6G\t1\t.\t.\t6\t6\tG\t66\t66\tG\t19\t.\t19\tnoncoding1_n_A6G_variant in ref and reads so should report\tgeneric description of noncoding1',
            'noncoding1\tnon_coding\t531\t72\tcluster_name\t120\t120\t95.87\tnoncoding1.scaffold.1\t234\t15.4\t1\tSNP\tn\tG9T\t0\t.\t.\t9\t9\tG\t69\t69\tG\t19\t.\t19\tnoncoding1_n_G9T_wild type in ref and reads\tgeneric description of noncoding1'
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

        c = cluster.Cluster(tmpdir, 'cluster_name', refdata, spades_other_options='--only-assembler', total_reads=64, total_reads_bases=3200)
        c.run()

        expected = [
            'presence_absence1\tpresence_absence\t539\t64\tcluster_name\t96\t96\t97.92\tpresence_absence1.scaffold.1\t213\t15.0\t1\tSNP\tp\tA10V\t1\tA10V\tNONSYN\t28\t28\tC\t83\t83\tT\t22\t.\t22\tpresence_absence1_p_A10V_Ref has wild, reads have variant so report\tGeneric description of presence_absence1',
            'presence_absence1\tpresence_absence\t539\t64\tcluster_name\t96\t96\t97.92\tpresence_absence1.scaffold.1\t213\t15.0\t0\t.\tp\t.\t0\t.\tSYN\t53\t53\tT\t108\t108\tC\t32\t.\t32\t.\tGeneric description of presence_absence1',

            'presence_absence1\tpresence_absence\t539\t64\tcluster_name\t96\t96\t97.92\tpresence_absence1.scaffold.1\t213\t15.0\t1\tSNP\tp\tR3S\t0\t.\t.\t7\t9\tC;G;C\t62\t64\tC;G;C\t18;17;17\t.;.;.\t18;17;17\tpresence_absence1_p_R3S_Ref and assembly have wild type\tGeneric description of presence_absence1',

            'presence_absence1\tpresence_absence\t539\t64\tcluster_name\t96\t96\t97.92\tpresence_absence1.scaffold.1\t213\t15.0\t1\tSNP\tp\tI5A\t1\t.\t.\t13\t15\tG;C;G\t68\t70\tG;C;G\t18;20;20\t.;.;.\t18;20;20\tpresence_absence1_p_I5A_Ref and reads have variant so report\tGeneric description of presence_absence1',
        ]

        self.assertEqual(expected, c.report_lines)
        shutil.rmtree(tmpdir)


    def test_full_run_ok_variants_only_variant_not_present(self):
        '''test complete run of cluster on a variants only gene when variant not present'''
        refdata = reference_data.ReferenceData(
            variants_only_fa=os.path.join(data_dir, 'cluster_test_full_run_ok_variants_only.fa'),
            metadata_tsv=os.path.join(data_dir, 'cluster_test_full_run_ok_variants_only.not_present.metadata.tsv'),
        )

        tmpdir = 'tmp.cluster_test_full_run_ok_variants_only.not_present'
        shutil.copytree(os.path.join(data_dir, 'cluster_test_full_run_ok_variants_only'), tmpdir)

        c = cluster.Cluster(tmpdir, 'cluster_name', refdata, spades_other_options='--only-assembler', total_reads=66, total_reads_bases=3300)
        c.run()
        expected = [
            'variants_only1\tvariants_only\t27\t66\tcluster_name\t96\t96\t100.0\tvariants_only1.scaffold.1\t215\t15.3\t1\tSNP\tp\tR3S\t0\t.\t.\t7\t9\tC;G;C\t65\t67\tC;G;C\t18;18;19\t.;.;.\t18;18;19\tvariants_only1_p_R3S_Ref and assembly have wild type, so do not report\tGeneric description of variants_only1'
        ]
        self.assertEqual(expected, c.report_lines)
        shutil.rmtree(tmpdir)


    def test_full_run_ok_variants_only_variant_not_present_always_report(self):
        '''test complete run of cluster on a variants only gene when variant not present but always report variant'''
        refdata = reference_data.ReferenceData(
            variants_only_fa=os.path.join(data_dir, 'cluster_test_full_run_ok_variants_only.fa'),
            metadata_tsv=os.path.join(data_dir, 'cluster_test_full_run_ok_variants_only.not_present.always_report.metadata.tsv'),
        )

        tmpdir = 'tmp.cluster_test_full_run_ok_variants_only.not_present.always_report'
        shutil.copytree(os.path.join(data_dir, 'cluster_test_full_run_ok_variants_only'), tmpdir)

        c = cluster.Cluster(tmpdir, 'cluster_name', refdata, spades_other_options='--only-assembler', total_reads=66, total_reads_bases=3300)
        c.run()
        expected = [
            'variants_only1\tvariants_only\t27\t66\tcluster_name\t96\t96\t100.0\tvariants_only1.scaffold.1\t215\t15.3\t1\tSNP\tp\tR3S\t0\t.\t.\t7\t9\tC;G;C\t65\t67\tC;G;C\t18;18;19\t.;.;.\t18;18;19\tvariants_only1_p_R3S_Ref and assembly have wild type, but always report anyway\tGeneric description of variants_only1'
        ]
        self.assertEqual(expected, c.report_lines)
        shutil.rmtree(tmpdir)


    def test_full_run_ok_variants_only_variant_is_present(self):
        '''test complete run of cluster on a variants only gene when variant is present'''
        refdata = reference_data.ReferenceData(
            variants_only_fa=os.path.join(data_dir, 'cluster_test_full_run_ok_variants_only.fa'),
            metadata_tsv=os.path.join(data_dir, 'cluster_test_full_run_ok_variants_only.present.metadata.tsv'),
        )

        tmpdir = 'tmp.cluster_test_full_run_ok_variants_only.present'
        shutil.copytree(os.path.join(data_dir, 'cluster_test_full_run_ok_variants_only'), tmpdir)

        c = cluster.Cluster(tmpdir, 'cluster_name', refdata, spades_other_options='--only-assembler', total_reads=66, total_reads_bases=3300)
        c.run()

        expected = [
            'variants_only1\tvariants_only\t27\t66\tcluster_name\t96\t96\t100.0\tvariants_only1.scaffold.1\t215\t15.3\t1\tSNP\tp\tR3S\t0\t.\t.\t7\t9\tC;G;C\t65\t67\tC;G;C\t18;18;19\t.;.;.\t18;18;19\tvariants_only1_p_R3S_Ref and assembly have wild type\tGeneric description of variants_only1',
            'variants_only1\tvariants_only\t27\t66\tcluster_name\t96\t96\t100.0\tvariants_only1.scaffold.1\t215\t15.3\t1\tSNP\tp\tI5A\t1\t.\t.\t13\t15\tG;C;G\t71\t73\tG;C;G\t17;17;17\t.;.;.\t17;17;17\tvariants_only1_p_I5A_Ref and reads have variant so report\tGeneric description of variants_only1',
        ]
        self.assertEqual(expected, c.report_lines)
        shutil.rmtree(tmpdir)
