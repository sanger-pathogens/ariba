import os
import openpyxl
import pyfastaq
from ariba import flag, common, reference_data, report

class Error (Exception): pass


int_columns = [
    'reads',
    'ref_len',
    'ref_base_assembled',
    'ctg_len',
    'ref_start',
    'ref_end',
    'ctg_start',
    'ctg_end',
]


float_columns = ['pc_ident']


class Summary:
    def __init__(
      self,
      outfile,
      filenames=None,
      fofn=None,
      filter_output=True,
      js_candy_prefix=None,
      min_id=90.0
    ):
        if filenames is None and fofn is None:
            raise Error('Error! Must supply filenames or fofn to Summary(). Cannot continue')

        if filenames is None:
            self.filenames = []
        else:
            self.filenames = filenames

        if fofn is not None:
            self.filenames.extend(self._load_fofn(fofn))

        self.filter_output = filter_output
        self.min_id = min_id
        self.outfile = outfile
        self.js_candy_prefix = js_candy_prefix


    def _load_fofn(self, fofn):
        f = pyfastaq.utils.open_file_read(fofn)
        filenames = [x.rstrip() for x in f.readlines()]
        pyfastaq.utils.close(f)
        return filenames


    def _check_files_exist(self):
        for fname in self.filenames:
            if not os.path.exists(fname):
                raise Error('File not found: "' + fname + '". Cannot continue')


    @classmethod
    def _line2dict(cls, line):
        data = line.rstrip().split('\t')
        if len(data) != len(report.columns):
            raise Error('Wrong number of columns in the following line. Expected ' + str(len(report.columns)) + ' but got ' + str(len(data)) + '\n' + line)
        d = {report.columns[i]: data[i] for i in range(len(data))}
        d['flag'] = flag.Flag(int(d['flag']) )
        for key in int_columns:
            try:
                d[key] = int(d[key])
            except:
                assert d[key] == '.'

        for key in float_columns:
            try:
                d[key] = float(d[key])
            except:
                assert d[key] == '.'

        return d


    @classmethod
    def _dict2key(cls, d):
        if d['var_type'] == '.':
            return d['ref_name'], '', ''
        elif d['known_var_change'] == d['ref_ctg_change'] == '.':
            raise Error('Unexpected data in ariba summary... \n' + d + '\n... known_var_change and ref_ctg_change both equal to ".", but var_type was not a ".". Cannot continue')
        else:
            if '.' not in [d['known_var_change'], d['ref_ctg_change']] and d['known_var_change'] != d['ref_ctg_change']:
                raise Error('Unexpected data in ariba summary... \n' + d + '\n... known_var_change != ref_ctg_change. Cannot continue')
            if d['known_var_change'] != '.':
                change = d['known_var_change']
            else:
                change = d['ref_ctg_change']

            return d['ref_name'], d['var_seq_type'], change


    @classmethod
    def _load_file(cls, filename):
        f = pyfastaq.utils.open_file_read(filename)
        d = {}

        for line in f:
            if line.startswith('#'):
                if line.rstrip()[1:].split('\t') != report.columns:
                    pyfastaq.utils.close(f)
                    raise Error('Error parsing the following line.\n' + line)
                continue
            data = Summary._line2dict(line)
            key = Summary._dict2key(data)
            if key[0] not in d:
                d[key[0]] = {}
            d[key[0]][key] = data

        pyfastaq.utils.close(f)
        return d


    @classmethod
    def _pc_id_of_longest(self, data_dict, seq_name):
        longest = 0
        identity = 0
        assert seq_name in data_dict

        for d in data_dict[seq_name].values():
            if d['ref_base_assembled'] > longest:
                longest = d['ref_base_assembled']
                identity = d['pc_ident']

        return identity


    def _to_summary_number_for_seq(self, l):
        f = l[0]['flag']
        if f.has('assembly_fail') or not f.has('assembled') or self._pc_id_of_longest(l) <= self.min_id:
            return 0

        if f.has('hit_both_strands') or (not f.has('complete_orf')):
            return 1

        if f.has('unique_contig') and f.has('assembled_into_one_contig') and f.has('complete_orf'):
            if f.has('has_nonsynonymous_variants'):
                return 3
            else:
                return 4
        else:
            return 2





    def _gather_output_rows(self):
        self.data = {filename: self._load_file(filename) for filename in self.filenames}

        all_genes = set()
        for l in self.data.values():
            all_genes.update(set(l.keys()))
        all_genes = list(all_genes)
        all_genes.sort()

        self.rows_out = []
        self.rows_out.append(['filename'] + all_genes)

        for filename in self.filenames:
            new_row = [filename]
            for gene in all_genes:
                if gene not in self.data[filename]:
                    new_row.append(0)
                else:
                    new_row.append(self._to_summary_number(self.data[filename][gene]))

            self.rows_out.append(new_row)


    def _filter_output_rows(self):
        if not self.filter_output:
            return

        # remove rows that are all zeros
        self.rows_out = [x for x in self.rows_out if x[1:] != [0]*(len(x)-1)]

        # remove columns that are all zeros
        to_remove = []
        for i in range(1, len(self.rows_out[0])):
            if sum([x[i] for x in self.rows_out[1:]]) == 0:
                to_remove.append(i)

        for i in range(len(self.rows_out)):
            self.rows_out[i] = [self.rows_out[i][j] for j in range(len(self.rows_out[i])) if j not in to_remove]



    def _write_tsv(self):
        f = pyfastaq.utils.open_file_write(self.outfile)
        print('#', end='', file=f)
        for row in self.rows_out:
            print('\t'.join([str(x) for x in row]), file=f)
        pyfastaq.utils.close(f)


    def _write_js_candy_csv(self, outfile):
        f = pyfastaq.utils.open_file_write(outfile)
        # js candy needs the "name" column.
        # Names must match those used in the tree file
        print('name', *self.rows_out[0][1:], sep=',', file=f)
        for row in self.rows_out[1:]:
            print(*row, sep=',', file=f)
        pyfastaq.utils.close(f)


    def _write_xls(self):
        workbook = openpyxl.Workbook()
        worksheet = workbook.worksheets[0]
        worksheet.title = 'ARIBA_summary'
        for row in self.rows_out:
            worksheet.append(row)
        workbook.save(self.outfile)


    @staticmethod
    def _distance_score_between_values(value1, value2):
        if value1 != value2 and 0 in [value1, value2]:
            return 1
        else:
            return 0


    @classmethod
    def _distance_score_between_lists(cls, scores1, scores2):
        assert len(scores1) == len(scores2)
        return sum([cls._distance_score_between_values(scores1[i], scores2[i]) for i in range(1, len(scores1))])


    def _write_distance_matrix(self, outfile):
        if len(self.rows_out) < 3:
            raise Error('Cannot calculate distance matrix to make tree for js_candy.\n' +
                        'Only one sample present.')

        if len(self.rows_out[0]) < 2:
            raise Error('Cannot calculate distance matrix to make tree for js_candy.\n' +
                        'No genes present in output')

        with open(outfile, 'w') as f:
            sample_names = [x[0] for x in self.rows_out]
            print(*sample_names[1:], sep='\t', file=f)

            for i in range(1,len(self.rows_out)):
                scores = []
                for j in range(2, len(self.rows_out)):
                    scores.append(self._distance_score_between_lists(self.rows_out[i], self.rows_out[j]))
                print(self.rows_out[i][0], *scores, sep='\t', file=f)


    @staticmethod
    def _newick_from_dist_matrix(distance_file, outfile):
        r_script = outfile + '.tmp.R'

        with open(r_script, 'w') as f:
            print('library(ape)', file=f)
            print('a=read.table("', distance_file, '", header=TRUE, row.names=1)', sep='', file=f)
            print('h=hclust(dist(a))', file=f)
            print('write.tree(as.phylo(h), file="', outfile, '")', sep='', file=f)

        common.syscall('R CMD BATCH --no-save ' + r_script)
        os.unlink(r_script + 'out')
        os.unlink(r_script)


    def _write_js_candy_files(self, outprefix):
        distance_file = outprefix + '.distance_matrix'
        tree_file = outprefix + '.tre'
        csv_file = outprefix + '.csv'
        self._write_distance_matrix(distance_file)
        self._newick_from_dist_matrix(distance_file, tree_file)
        os.unlink(distance_file)
        self._write_js_candy_csv(csv_file)


    def run(self):
        self._check_files_exist()
        self._gather_output_rows()
        self._filter_output_rows()
        if self.outfile.endswith('.xls'):
            self._write_xls()
        else:
            self._write_tsv()

        if self.js_candy_prefix is not None:
            self._write_js_candy_files(self.js_candy_prefix)
