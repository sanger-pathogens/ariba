import unittest
import os
import filecmp
from ariba import report_flag_expander

modules_dir = os.path.dirname(os.path.abspath(report_flag_expander.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


class TestReportFlagExpander(unittest.TestCase):
    def test_run(self):
        '''test run'''
        infile = os.path.join(data_dir, 'report_flag_expander.run.in.tsv')
        expected = os.path.join(data_dir, 'report_flag_expander.run.out.tsv')
        tmp_out = 'tmp.report_flag_expander.out.tsv'
        expander = report_flag_expander.ReportFlagExpander(infile, tmp_out)
        expander.run()
        self.assertTrue(filecmp.cmp(expected, tmp_out, shallow=False))
        os.unlink(tmp_out)

