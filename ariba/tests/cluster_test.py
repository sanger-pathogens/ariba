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


def file2lines(filename):
    f = pyfastaq.utils.open_file_read(filename)
    lines = f.readlines()
    pyfastaq.utils.close(f)
    return lines


def load_gene(filename):
    file_reader = pyfastaq.sequences.file_reader(filename)
    seq = None
    for seq in file_reader:
        pass
    return seq


class TestCluster(unittest.TestCase):
    def test_init_fail_files_missing(self):
        '''test init_fail_files_missing'''
        dirs = [
            'cluster_test_directorynotexist'
            'cluster_test_init_no_genes_fa',
            'cluster_test_init_no_reads_1',
            'cluster_test_init_no_reads_2',
        ]
        dirs = [os.path.join(data_dir, d) for d in dirs]
        for d in dirs:
            clean_cluster_dir(d)
            with self.assertRaises(cluster.Error):
                c = cluster.Cluster(d, 'name')
            clean_cluster_dir(d)


    def test_count_reads(self):
        '''test _count_reads pass'''
        reads1 = os.path.join(data_dir, 'cluster_test_count_reads_1.fq')
        reads2 = os.path.join(data_dir, 'cluster_test_count_reads_2.fq')
        self.assertEqual(4, cluster.Cluster._count_reads(reads1, reads2))


    def test_make_report_lines_nonsynonymous(self):
        '''test _make_report_lines'''
        cluster_dir = os.path.join(data_dir, 'cluster_test_generic')
        clean_cluster_dir(cluster_dir)
        c = cluster.Cluster(cluster_dir, 'cluster_name')
        c.gene = pyfastaq.sequences.Fasta('gene', 'GATCGCGAAGCGATGACCCATGAAGCGACCGAACGCTGA')
        v1 = pymummer.variant.Variant(pymummer.snp.Snp('8\tA\tG\t8\tx\tx\t39\t39\tx\tx\tgene\tcontig'))

        nucmer_hit = ['1', '10', '1', '10', '10', '10', '90.00', '1000', '1000', '1', '1', 'gene', 'contig']
        c.nucmer_hits = {'contig': [pymummer.alignment.Alignment('\t'.join(nucmer_hit))]}
        c.mummer_variants = {'contig': [[v1]]}
        c.percent_identities = {'contig': 92.42}
        c.status_flag.set_flag(42)
        c.assembled_ok = True
        c.final_assembly_read_depths = os.path.join(data_dir, 'cluster_test_make_report_lines.read_depths.gz')
        c._make_report_lines()
        expected = [[
            'gene',
            554,
            2,
            'cluster_name',
            39,
            10,
            92.42,
            'SNP',
            'NONSYN',
            'E3G',
            8,
            8,
            'A',
            'contig',
            39,
            8,
            8,
            'G',
            '.',
            '.',
            '.'
        ]]
        self.assertEqual(expected, c.report_lines)
        clean_cluster_dir(cluster_dir)


    def test_make_report_lines_synonymous(self):
        '''test _make_report_lines'''
        cluster_dir = os.path.join(data_dir, 'cluster_test_generic')
        clean_cluster_dir(cluster_dir)
        c = cluster.Cluster(cluster_dir, 'cluster_name')
        c.gene = pyfastaq.sequences.Fasta('gene', 'GATCGCGAAGCGATGACCCATGAAGCGACCGAACGCTGA')
        v1 = pymummer.variant.Variant(pymummer.snp.Snp('6\tC\tT\t6\tx\tx\t39\t39\tx\tx\tgene\tcontig'))

        nucmer_hit = ['1', '10', '1', '10', '10', '10', '90.00', '1000', '1000', '1', '1', 'gene', 'contig']
        c.nucmer_hits = {'contig': [pymummer.alignment.Alignment('\t'.join(nucmer_hit))]}
        c.mummer_variants = {'contig': [[v1]]}
        c.percent_identities = {'contig': 92.42}
        c.status_flag.set_flag(42)
        c.assembled_ok = True
        c.final_assembly_read_depths = os.path.join(data_dir, 'cluster_test_make_report_lines.read_depths.gz')
        c._make_report_lines()
        expected = [[
            'gene',
            42,
            2,
            'cluster_name',
            39,
            10,
            92.42,
            'SNP',
            'SYN',
            '.',
            6,
            6,
            'C',
            'contig',
            39,
            6,
            6,
            'T',
            42,
            'G',
            '22,20'
        ]]
        self.assertEqual(expected, c.report_lines)
        clean_cluster_dir(cluster_dir)


    def test_make_report_lines_assembly_fail(self):
        '''test _make_report_lines when assembly fails'''
        cluster_dir = os.path.join(data_dir, 'cluster_test_generic')
        clean_cluster_dir(cluster_dir)
        refdata = reference_data.ReferenceData(presence_absence_fa=os.path.join(cluster_dir, 'genes.fa'))
        c = cluster.Cluster(cluster_dir, 'cluster_name', refdata=refdata)
        c.gene = pyfastaq.sequences.Fasta('gene', 'GATCGCGAAGCGATGACCCATGAAGCGACCGAACGCTGA')
        c.status_flag.set_flag(64)
        c.assembled_ok = False
        c._make_report_lines()
        expected = [
            [
                'gene',
                64,
                2,
                'cluster_name',
                39,
            ] + ['.'] * 16
        ]
        self.assertEqual(expected, c.report_lines)
        clean_cluster_dir(cluster_dir)


# to test:
# count_reads()

