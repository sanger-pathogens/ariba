import os
import pyfastaq
from ariba import report

class Error (Exception): pass

class ReportFilter:
    def __init__(self, infile):
        self.report = self._load_report(infile)


    @classmethod
    def _report_line_to_dict(cls, line):
        '''Takes report line string as input. Returns a dict of column name -> value in line'''
        data = line.split('\t')
        if len(data) != len(report.columns):
            return None

        return dict(zip(report.columns, data))


    @staticmethod
    def _load_report(infile):
        '''Loads report file into a dictionary. Key=refrence name.
        Value = list of report lines for that reference'''
        report_dict = {}
        f = pyfastaq.utils.open_file_read(infile)
        first_line = True

        for line in f:
            line = line.rstrip()

            if first_line:
                expected_first_line = '#' + '\t'.join(report.columns)
                if line != expected_first_line:
                    pyfastaq.utils.close(f)
                    raise Error('Error reading report file. Expected first line of file is\n' + expected_first_line + '\nbut got:\n' + line)
                first_line = False
            else:
                line_dict = ReportFilter._report_line_to_dict(line)
                if line_dict is None:
                    pyfastaq.utils.close(f)
                    raise Error('Error reading report file. Expected ' + str(len(report.columns)) + ' columns but got ' + str(len(data)) + ' columns at this line:\n' + line)
                ref_name = line_dict['ref_name']
                if ref_name not in report_dict:
                    report_dict[ref_name] = []
                report_dict[ref_name].append(line_dict)

        pyfastaq.utils.close(f)
        return report_dict

