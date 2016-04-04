import os
import pyfastaq
from ariba import report, flag

class Error (Exception): pass

class ReportFilter:
    def __init__(self,
            infile=None,
            min_pc_ident=90,
            min_ref_base_assembled=1,
            ignore_not_has_known_variant=True,
        ):

        if infile is not None:
            self.report = self._load_report(infile)
        else:
            self.report = {}

        self.min_pc_ident = min_pc_ident
        self.min_ref_base_assembled = min_ref_base_assembled
        self.ignore_not_has_known_variant = ignore_not_has_known_variant


    @classmethod
    def _report_line_to_dict(cls, line):
        '''Takes report line string as input. Returns a dict of column name -> value in line'''
        data = line.split('\t')
        if len(data) != len(report.columns):
            return None

        d = dict(zip(report.columns, data))
        for key in report.int_columns:
            d[key] = int(d[key])

        for key in report.float_columns:
            d[key] = float(d[key])

        d['flag'] = flag.Flag(int(d['flag']))
        return d


    @classmethod
    def _dict_to_report_line(cls, report_dict):
        '''Takes a report_dict as input and returns a report line'''
        return '\t'.join([str(report_dict[x]) for x in report.columns])


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


    @staticmethod
    def _report_dict_passes_known_variant_filter(ignore_not_has_known_variant, report_dict):
        if ignore_not_has_known_variant:
            return report_dict['has_known_var'] == '1'
        else:
            return True


    def _report_dict_passes_filters(self, report_dict):
        return report_dict['pc_ident'] >= self.min_pc_ident \
                   and report_dict['ref_base_assembled'] >= self.min_ref_base_assembled \
                   and self._report_dict_passes_known_variant_filter(self.ignore_not_has_known_variant, report_dict)


    def _filter_dicts(self):
        '''Filters out all the report_dicts that do not pass the cutoffs. If any ref sequence
           loses all of its report_dicts, then it is completely removed.'''
        keys_to_remove = set()

        for ref_name in self.report:
            self.report[ref_name] = [x for x in self.report[ref_name] if self._report_dict_passes_filters(x)]
            if len(self.report[ref_name]) == 0:
                keys_to_remove.add(ref_name)

        for key in keys_to_remove:
            del self.report[key]


    def _write_report_tsv(self, outfile):
        f = pyfastaq.utils.open_file_write(outfile)
        print('#' + '\t'.join(report.columns), file=f)

        for key, report_dicts in sorted(self.report.items()):
            for d in report_dicts:
                print(ReportFilter._dict_to_report_line(d), file=f)

        pyfastaq.utils.close(f)


    def run(self, outfile):
        self._filter_dicts()
        self._write_report_tsv(outfile)
