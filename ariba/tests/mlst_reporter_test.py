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
        expected_out_all = os.path.join(data_dir, 'mlst_reporter.all_present_perfect.report.out.all.tsv')
        tmp_out = 'tmp.mlst_reporter.test_all_present_perfect'
        got_simple = tmp_out + '.tsv'
        got_all = tmp_out + '.all.tsv'
        reporter = mlst_reporter.MlstReporter(report_in, profile_in, tmp_out)
        reporter.run()
        self.assertTrue(filecmp.cmp(expected_out_simple, got_simple, shallow=False))
        self.assertTrue(filecmp.cmp(expected_out_all, got_all, shallow=False))
        os.unlink(got_simple)
        os.unlink(got_all)
