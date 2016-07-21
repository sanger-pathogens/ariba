import unittest
import os
import shutil
import filecmp
from ariba import external_progs, ref_preparer

modules_dir = os.path.dirname(os.path.abspath(ref_preparer.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


class TestRefPreparer(unittest.TestCase):
    def test_run(self):
        '''test run'''
        fasta_in = [
            os.path.join(data_dir, 'ref_preparer_test_run.in.1.fa'),
            os.path.join(data_dir, 'ref_preparer_test_run.in.2.fa'),
            os.path.join(data_dir, 'ref_preparer_test_run.in.3.fa'),
        ]
        tsv_in = [
            os.path.join(data_dir, 'ref_preparer_test_run.in.1.tsv'),
            os.path.join(data_dir, 'ref_preparer_test_run.in.2.tsv'),
        ]

        extern_progs = external_progs.ExternalProgs()
        refprep = ref_preparer.RefPreparer(fasta_in, tsv_in, extern_progs, genetic_code=1)
        tmp_out = 'tmp.ref_preparer_test_run'
        refprep.run(tmp_out)
        expected_outdir = os.path.join(data_dir, 'ref_preparer_test_run.out')

        test_files = [
            '01.filter.check_metadata.tsv',
            '01.filter.check_genes.log',
            '01.filter.check_metadata.log',
            '02.cdhit.all.fa',
            '02.cdhit.clusters.tsv',
            '02.cdhit.gene.fa',
            '02.cdhit.gene.varonly.fa',
            '02.cdhit.noncoding.fa',
            '02.cdhit.noncoding.varonly.fa',
        ]

        for filename in test_files:
            expected = os.path.join(expected_outdir, filename)
            got = os.path.join(tmp_out, filename)
            self.assertTrue(filecmp.cmp(expected, got, shallow=False))

        shutil.rmtree(tmp_out)

