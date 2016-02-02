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

        expected = '\t'.join(['.', '.', '1088', '2', 'cluster_name'] + ['.'] * 18)
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

        expected = '\t'.join(['noncoding_ref_seq', 'non_coding', '64', '4', 'cluster_name'] + ['.'] * 18)
        self.assertEqual([expected], c.report_lines)
        self.assertFalse(c.status_flag.has('ref_seq_choose_fail'))
        self.assertTrue(c.status_flag.has('assembly_fail'))
        shutil.rmtree(tmpdir)


    def test_full_run_ok_non_coding(self):
        '''test complete run of cluster on a noncoding sequence'''
        refdata = reference_data.ReferenceData(
            non_coding_fa=os.path.join(data_dir, 'cluster_test_full_run_ok_non_coding.fa')
        )

        tmpdir = 'tmp.test_full_run_ok_non_coding'
        shutil.copytree(os.path.join(data_dir, 'cluster_test_full_run_ok_non_coding'), tmpdir)

        c = cluster.Cluster(tmpdir, 'cluster_name', refdata, spades_other_options='--only-assembler')
        c.run()

        expected = ['noncoding_ref1\tnon_coding\t19\t154\tcluster_name\t400\t' +  '\t'.join(['.'] * 17)]
        self.assertEqual(expected, c.report_lines)
        shutil.rmtree(tmpdir)

    #FIXME more tests with full runs and variants!
