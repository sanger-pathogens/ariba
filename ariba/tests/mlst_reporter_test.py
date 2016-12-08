import unittest
import os
import filecmp
from ariba import mlst_reporter

modules_dir = os.path.dirname(os.path.abspath(mlst_reporter.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


class TestMlstReporter(unittest.TestCase):
    def test_all_present_perfect(self):
        '''test when all alleles 100% there no variants'''
        profile_in = os.path.join(data_dir, 'mlst_reporter.profile.in.tsv')
        report_in = os.path.join(data_dir, 'mlst_reporter.all_present_perfect.report.in.tsv')
        expected_out_simple = os.path.join(data_dir, 'mlst_reporter.all_present_perfect.report.out.tsv')
        expected_out_all = os.path.join(data_dir, 'mlst_reporter.all_present_perfect.report.out.details.tsv')
        tmp_out = 'tmp.mlst_reporter.test_all_present_perfect'
        got_simple = tmp_out + '.tsv'
        got_details = tmp_out + '.details.tsv'
        reporter = mlst_reporter.MlstReporter(report_in, profile_in, tmp_out)
        reporter.run()
        self.assertTrue(filecmp.cmp(expected_out_simple, got_simple, shallow=False))
        self.assertTrue(filecmp.cmp(expected_out_all, got_details, shallow=False))
        os.unlink(got_simple)
        os.unlink(got_details)


    def test_new_set(self):
        '''test when allele combination is new'''
        profile_in = os.path.join(data_dir, 'mlst_reporter.profile.in.tsv')
        report_in = os.path.join(data_dir, 'mlst_reporter.new_st.report.in.tsv')
        expected_out_simple = os.path.join(data_dir, 'mlst_reporter.new_st.report.out.tsv')
        expected_out_all = os.path.join(data_dir, 'mlst_reporter.new_st.report.out.details.tsv')
        tmp_out = 'tmp.mlst_reporter.test_new_st'
        got_simple = tmp_out + '.tsv'
        got_details = tmp_out + '.details.tsv'
        reporter = mlst_reporter.MlstReporter(report_in, profile_in, tmp_out)
        reporter.run()
        self.assertTrue(filecmp.cmp(expected_out_simple, got_simple, shallow=False))
        self.assertTrue(filecmp.cmp(expected_out_all, got_details, shallow=False))
        os.unlink(got_simple)
        os.unlink(got_details)


    def test_one_gene_missing(self):
        '''test when one gene missing'''
        profile_in = os.path.join(data_dir, 'mlst_reporter.profile.in.tsv')
        report_in = os.path.join(data_dir, 'mlst_reporter.one_gene_missing.in.tsv')
        expected_out_simple = os.path.join(data_dir, 'mlst_reporter.one_gene_missing.out.tsv')
        expected_out_all = os.path.join(data_dir, 'mlst_reporter.one_gene_missing.out.details.tsv')
        tmp_out = 'tmp.mlst_reporter.test_one_gene_missing'
        got_simple = tmp_out + '.tsv'
        got_details = tmp_out + '.details.tsv'
        reporter = mlst_reporter.MlstReporter(report_in, profile_in, tmp_out)
        reporter.run()
        self.assertTrue(filecmp.cmp(expected_out_simple, got_simple, shallow=False))
        self.assertTrue(filecmp.cmp(expected_out_all, got_details, shallow=False))
        os.unlink(got_simple)
        os.unlink(got_details)


    def test_het_snp(self):
        '''test when one gene has a het snp'''
        profile_in = os.path.join(data_dir, 'mlst_reporter.profile.in.tsv')
        report_in = os.path.join(data_dir, 'mlst_reporter.het_snp.in.tsv')
        expected_out_simple = os.path.join(data_dir, 'mlst_reporter.het_snp.out.tsv')
        expected_out_all = os.path.join(data_dir, 'mlst_reporter.het_snp.out.details.tsv')
        tmp_out = 'tmp.mlst_reporter.test_het_snp'
        got_simple = tmp_out + '.tsv'
        got_details = tmp_out + '.details.tsv'
        reporter = mlst_reporter.MlstReporter(report_in, profile_in, tmp_out)
        reporter.run()
        self.assertTrue(filecmp.cmp(expected_out_simple, got_simple, shallow=False))
        self.assertTrue(filecmp.cmp(expected_out_all, got_details, shallow=False))
        os.unlink(got_simple)
        os.unlink(got_details)


    def test_het_snps(self):
        '''test when one gene has two het snps'''
        profile_in = os.path.join(data_dir, 'mlst_reporter.profile.in.tsv')
        report_in = os.path.join(data_dir, 'mlst_reporter.het_snps.in.tsv')
        expected_out_simple = os.path.join(data_dir, 'mlst_reporter.het_snps.out.tsv')
        expected_out_all = os.path.join(data_dir, 'mlst_reporter.het_snps.out.details.tsv')
        tmp_out = 'tmp.mlst_reporter.test_het_snps'
        got_simple = tmp_out + '.tsv'
        got_details = tmp_out + '.details.tsv'
        reporter = mlst_reporter.MlstReporter(report_in, profile_in, tmp_out)
        reporter.run()
        self.assertTrue(filecmp.cmp(expected_out_simple, got_simple, shallow=False))
        self.assertTrue(filecmp.cmp(expected_out_all, got_details, shallow=False))
        os.unlink(got_simple)
        os.unlink(got_details)
