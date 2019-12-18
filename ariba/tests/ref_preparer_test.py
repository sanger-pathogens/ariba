import unittest
import os
import filecmp
from ariba import common, external_progs, ref_preparer

modules_dir = os.path.dirname(os.path.abspath(ref_preparer.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


class TestRefPreparer(unittest.TestCase):
    def test_fasta_to_metadata(self):
        '''test _fasta_to_metadata'''
        infile = os.path.join(data_dir, 'ref_preparer_test_fasta_to_metadata.fa')
        tmp_out = 'tmp.test_fasta_to_metadata.tsv'
        expected_coding = os.path.join(data_dir, 'ref_preparer_test_fasta_to_metadata.coding.tsv')
        expected_noncoding = os.path.join(data_dir, 'ref_preparer_test_fasta_to_metadata.noncoding.tsv')

        with open(tmp_out, 'w') as f:
            ref_preparer.RefPreparer._fasta_to_metadata(infile, f, True)
        self.assertTrue(filecmp.cmp(expected_coding, tmp_out, shallow=False))

        with open(tmp_out, 'w') as f:
            ref_preparer.RefPreparer._fasta_to_metadata(infile, f, False)
        self.assertTrue(filecmp.cmp(expected_noncoding, tmp_out, shallow=False))

        os.unlink(tmp_out)


    def test_rename_clusters(self):
        '''test _rename_clusters'''
        clusters_in = {
           '0': {'no_dot_in_name'},
           '1': {'another_no_dot_in_name'},
           '2': {'foo.blah_blah_blah', 'foo.xyz'},
           '3': {'foo.abc', 'foo.def'},
           '4': {'pre1.abc', 'pre2.abc'},
           '5': {'pre1.def', 'pre2.pqr', 'pre2.zxy'},
           '6': {'prefix1.abc', 'prefix1.def', 'something_else.abc'},
           '7': {'prefix1.fgh', 'prefix1.ijk', 'something_else_again.abc'},
           '8': {'xyz.1', 'xyz.2', 'abcdefgh'},
           '9': {'a.foo', 'a.bar'},
           '10': {'abc_1.1'},
           '11': {'abc.2'},
           '12': {'abc.3'},
           '13': {'abc.4'},
           '14': {'def.1'},
           '15': {'def_1.2'},
           '16': {'def_2.3'},
           '17': {'def.4'},
           '18': {'def.5'},
           '19': {'x_1.foo'},
           '20': {'x_1.bar'},
           '21': {'x_1.baz'},
           '22': {'x_1_2.abc'},
           '23': {'x_1_2.def'},
           '24': {'y_1.foo'},
           '25': {'y_1_2.def'},
           '26': {'y_1.bar'},
           '27': {'y_1.baz'},
           '28': {'y_1_2.abc'},
        }

        expected = {
            'cluster': {'no_dot_in_name'},
            'cluster_1': {'another_no_dot_in_name'},
            'foo': {'foo.blah_blah_blah', 'foo.xyz'},
            'foo_1': {'foo.abc', 'foo.def'},
            'pre-': {'pre1.abc', 'pre2.abc'},
            'pre-_1': {'pre1.def', 'pre2.pqr', 'pre2.zxy'},
            'prefix1+': {'prefix1.abc', 'prefix1.def', 'something_else.abc'},
            'prefix1+_1': {'prefix1.fgh', 'prefix1.ijk', 'something_else_again.abc'},
            'xyz+': {'xyz.1', 'xyz.2', 'abcdefgh'},
            'cluster_2': {'a.foo', 'a.bar'},
            'abc_1': {'abc_1.1'},
            'abc': {'abc.2'},
            'abc_2': {'abc.3'},
            'abc_3': {'abc.4'},
            'def_1': {'def_1.2'},
            'def_2': {'def_2.3'},
            'def': {'def.1'},
            'def_3': {'def.4'},
            'def_4': {'def.5'},
            'x_1': {'x_1.foo'},
            'x_1_1': {'x_1.bar'},
            'x_1_2_1': {'x_1_2.abc'},
            'x_1_2_2': {'x_1_2.def'},
            'x_1_2': {'x_1.baz'},
            'y_1': {'y_1.foo'},
            'y_1_1': {'y_1.bar'},
            'y_1_2_1': {'y_1_2.abc'},
            'y_1_2': {'y_1_2.def'},
            'y_1_3': {'y_1.baz'},
        }

        got = ref_preparer.RefPreparer._rename_clusters(clusters_in)
        self.assertEqual(expected, got)


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
        refprep = ref_preparer.RefPreparer(fasta_in, extern_progs, metadata_tsv_files=tsv_in, genetic_code=1)
        tmp_out = 'tmp.ref_preparer_test_run'
        refprep.run(tmp_out)
        expected_outdir = os.path.join(data_dir, 'ref_preparer_test_run.out')

        test_files = [
            '01.filter.check_metadata.tsv',
            '01.filter.check_genes.log',
            '01.filter.check_noncoding.log',
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

        common.rmtree(tmp_out)


    def test_run_all_noncoding(self):
        '''test run with no metadata input, all sequences are noncoding'''
        fasta_in = [
            os.path.join(data_dir, 'ref_preparer_test_run.in.1.fa'),
            os.path.join(data_dir, 'ref_preparer_test_run.in.2.fa'),
            os.path.join(data_dir, 'ref_preparer_test_run.in.3.fa'),
        ]

        extern_progs = external_progs.ExternalProgs()
        refprep = ref_preparer.RefPreparer(fasta_in, extern_progs, all_coding='no', genetic_code=1)
        tmp_out = 'tmp.ref_preparer_test_run'
        refprep.run(tmp_out)
        expected_outdir = os.path.join(data_dir, 'ref_preparer_test_run_all_noncoding.out')

        test_files = [
            '00.auto_metadata.tsv',
            '01.filter.check_metadata.tsv',
            '01.filter.check_genes.log',
            '01.filter.check_noncoding.log',
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

        common.rmtree(tmp_out)

    def test_run_noncoding_checks(self):
        '''test run with noncoding sequences that are outside of the allowed size range'''
        fasta_in = [
            os.path.join(data_dir, 'ref_preparer_test_run.in.4.fa')
        ]
        tsv_in = [
            os.path.join(data_dir, 'ref_preparer_test_run.in.4.tsv')
        ]

        extern_progs = external_progs.ExternalProgs()
        refprep = ref_preparer.RefPreparer(
            fasta_in, extern_progs, min_noncoding_length=6, max_noncoding_length=20,
            metadata_tsv_files=tsv_in, genetic_code=1)
        tmp_out = 'tmp.ref_preparer_test_run_noncoding_checks'
        refprep.run(tmp_out)
        expected_outdir = os.path.join(data_dir, 'ref_preparer_test_run_noncoding_checks.out')

        test_files = [
            '01.filter.check_metadata.tsv',
            '01.filter.check_genes.log',
            '01.filter.check_noncoding.log',
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

        common.rmtree(tmp_out)


