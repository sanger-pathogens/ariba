import pyfastaq
from ariba import report, flag

class Error (Exception): pass

class ReportFilter:
    def __init__(self,
            infile=None,
            min_pc_ident=90,
            min_ref_base_assembled=1,
            ignore_not_has_known_variant=False,
            remove_synonymous_snps=True,
            exclude_flags=None,
        ):

        if infile is not None:
            self.report = self._load_report(infile)
        else:
            self.report = {}

        self.min_pc_ident = min_pc_ident
        self.min_ref_base_assembled = min_ref_base_assembled
        self.ignore_not_has_known_variant = ignore_not_has_known_variant
        self.remove_synonymous_snps = remove_synonymous_snps

        if exclude_flags is None:
            self.exclude_flags = ['assembly_fail', 'ref_seq_choose_fail']
        else:
            self.exclude_flags = exclude_flags


    @classmethod
    def _report_line_to_dict(cls, line):
        '''Takes report line string as input. Returns a dict of column name -> value in line'''
        data = line.split('\t')
        if len(data) != len(report.columns):
            return None

        d = dict(zip(report.columns, data))
        for key in report.int_columns:
            try:
                d[key] = int(d[key])
            except:
                assert d[key] == '.'

        for key in report.float_columns:
            try:
                d[key] = float(d[key])
            except:
                assert d[key] == '.'

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
                    raise Error('Error reading report file at this line:\n' + line)
                ref_name = line_dict['ref_name']
                ctg_name = line_dict['ctg']
                if ref_name not in report_dict:
                    report_dict[ref_name] = {}
                if ctg_name not in report_dict[ref_name]:
                    report_dict[ref_name][ctg_name] = []

                report_dict[ref_name][ctg_name].append(line_dict)

        pyfastaq.utils.close(f)
        return report_dict


    @staticmethod
    def _flag_passes_filter(flag, exclude_flags):
        for f in exclude_flags:
            if flag.has(f):
                return False
        return True


    def _report_dict_passes_non_essential_filters(self, report_dict):
        # known_var == '.' iff this line is not reporting a variant. Which means it passes all non-essential filters
        if report_dict.get('known_var', '.') == '.':
            return True

        if self.remove_synonymous_snps and report_dict.get('ref_ctg_effect', None) == 'SYN':
            return False

        if self.ignore_not_has_known_variant and report_dict['known_var'] == '1' and report_dict['has_known_var'] == '0':
            return False

        return True


    def _report_dict_passes_essential_filters(self, report_dict):
        return ReportFilter._flag_passes_filter(report_dict['flag'], self.exclude_flags) \
                   and report_dict['pc_ident'] >= self.min_pc_ident \
                   and report_dict['ref_base_assembled'] >= self.min_ref_base_assembled


    def _filter_list_of_dicts(self, dicts_list):
        if len(dicts_list) == 0:
            return []

        pass_dicts = []
        essential_dicts = []
        fail_dicts = []

        for d in dicts_list:
            if self._report_dict_passes_essential_filters(d):
                if self._report_dict_passes_non_essential_filters(d):
                    pass_dicts.append(d)
                else:
                    essential_dicts.append(d)
            else:
                fail_dicts.append(d)

        if len(pass_dicts) == 0:
            assert len(fail_dicts) + len(essential_dicts) > 0
            if len(essential_dicts) > 0:
                new_d = essential_dicts[0]
                for key in report.var_columns:
                    new_d[key] = '.'
                pass_dicts.append(new_d)

        return ReportFilter._remove_all_after_first_frameshift(pass_dicts)


    @staticmethod
    def _remove_all_after_first_frameshift(dicts_list):
        fshift_starts = [int(d['ref_start']) for d in dicts_list if d.get('ref_ctg_effect', None) == 'FSHIFT']
        if len(fshift_starts) == 0:
            return dicts_list

        first_start = min(fshift_starts)

        return [d for d in dicts_list if d['ref_start'] == '.' or d['ref_start'] <= first_start]


    def _filter_dicts(self):
        '''Filters out all the report_dicts that do not pass the cutoffs. If any ref sequence
           loses all of its report_dicts, then it is completely removed.'''
        keys_to_remove = set()

        for ref_name in self.report:
            for ctg_name in self.report[ref_name]:
                self.report[ref_name][ctg_name] = self._filter_list_of_dicts(self.report[ref_name][ctg_name])
                if len(self.report[ref_name][ctg_name]) == 0:
                    keys_to_remove.add((ref_name, ctg_name))

        refs_to_remove = set()

        for ref_name, ctg_name in keys_to_remove:
            del self.report[ref_name][ctg_name]
            if len(self.report[ref_name]) == 0:
                refs_to_remove.add(ref_name)

        for ref_name in refs_to_remove:
            del self.report[ref_name]


    def _write_report_tsv(self, outfile):
        f = pyfastaq.utils.open_file_write(outfile)
        print('#' + '\t'.join(report.columns), file=f)

        for ref_name in sorted(self.report):
            for ctg_name, report_dicts in sorted(self.report[ref_name].items()):
                for d in report_dicts:
                    print(ReportFilter._dict_to_report_line(d), file=f)

        pyfastaq.utils.close(f)


    def run(self, outfile):
        self._filter_dicts()
        self._write_report_tsv(outfile)

