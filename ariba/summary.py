import os
import re
import openpyxl
import pyfastaq
from ariba import flag, common, report, summary_cluster

class Error (Exception): pass

class Summary:
    def __init__(
      self,
      outfile,
      filenames=None,
      fofn=None,
      filter_output=True,
      phandango_prefix=None,
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
        self.phandango_prefix = phandango_prefix
        self.clusters = {}


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
    def _dict2key(cls, d):
        if d['var_type'] == '.':
            return d['ref_name'], '', ''
        elif d['known_var_change'] == d['ref_ctg_change'] == '.':
            raise Error('Unexpected data in ariba summary... \n' + str(d) + '\n... known_var_change and ref_ctg_change both equal to ".", but var_type was not a ".". Cannot continue')
        else:
            if '.' not in [d['known_var_change'], d['ref_ctg_change']] and d['known_var_change'] != d['ref_ctg_change']:
                raise Error('Unexpected data in ariba summary... \n' + str(d) + '\n... known_var_change != ref_ctg_change. Cannot continue')
            if d['known_var_change'] != '.':
                change = d['known_var_change']
            else:
                change = d['ref_ctg_change']

            return d['ref_name'], d['var_seq_type'], change


    @classmethod
    def _load_file(cls, filename):
        f = pyfastaq.utils.open_file_read(filename)
        clusters = {}

        for line in f:
            if line.startswith('#'):
                if line.rstrip()[1:].split('\t') != report.columns:
                    pyfastaq.utils.close(f)
                    raise Error('Error parsing the following line.\n' + line)
                continue
            data = summary_cluster.SummaryCluster.line2dict(line)
            cluster = data['cluster']

            if cluster not in clusters:
                clusters[cluster] = summary_cluster.SummaryCluster()

            clusters[cluster].add_data_dict(data)

        pyfastaq.utils.close(f)
        return clusters


    @classmethod
    def _gather_output_rows(cls, filenames, min_id):
        data = {filename: Summary._load_file(filename) for filename in filenames}

        all_column_tuples = set()

        for filename, data_dict in data.items():
            for seq_name, seq_data_dict in data_dict.items():
                all_column_tuples.update(set(seq_data_dict.keys()))
                all_column_tuples.add((seq_name, '', ''))



        all_column_tuples = list(all_column_tuples)
        all_column_tuples.sort()
        rows = [['filename']]
        for t in all_column_tuples:
            if t[1] == t[2] == '':
                rows[0].append(t[0])
            else:
                rows[0].append(t[0] + ';' + 'var.' + t[1] + '.' + t[2])

        for filename in filenames:
            new_row = [filename]
            for column_tuple in all_column_tuples:
                if column_tuple[0] not in data[filename]:
                    new_row.append(0)
                elif column_tuple[1] == '':
                    new_row.append(Summary._to_summary_number_for_seq(data[filename], column_tuple[0], min_id))
                elif column_tuple in data[filename][column_tuple[0]]:
                    new_row.append(Summary._to_summary_number_for_variant(data[filename][column_tuple[0]][column_tuple]))
                else:
                    new_row.append(0)

            rows.append(new_row)

        return rows


    @classmethod
    def _filter_output_rows(cls, rows):
        # remove rows that are all zeros
        rows = [x for x in rows if x[1:] != [0]*(len(x)-1)]

        # remove columns that are all zeros
        to_remove = []
        for i in range(1, len(rows[0])):
            if sum([x[i] for x in rows[1:]]) == 0:
                to_remove.append(i)

        for i in range(len(rows)):
            rows[i] = [rows[i][j] for j in range(len(rows[i])) if j not in to_remove]

        return rows


    @classmethod
    def _write_tsv(cls, rows, outfile):
        f = pyfastaq.utils.open_file_write(outfile)
        print('#', end='', file=f)
        for row in rows:
            print('\t'.join([str(x) for x in row]), file=f)
        pyfastaq.utils.close(f)


    @classmethod
    def _write_phandango_csv(cls, rows, outfile):
        # phandango needs the "name" column.
        # Names must match those used in the tree file.
        # we also need to add suffixes like :z1 to make phandango colour
        # the columns consistently. We want to colour just sequence
        # columns the same, and then all variant columns the same
        header_line = ['name']
        var_regex = re.compile('^.*;var\.[np.]\.\S+$')
        var_columns = []

        for i in range(1, len(rows[0]), 1):
            heading = rows[0][i]
            if var_regex.search(heading) is not None:
                var_columns.append(i)

            header_line.append(heading + ':o1')

        f = pyfastaq.utils.open_file_write(outfile)
        print(*header_line, sep=',', file=f)
        for row in rows[1:]:
            to_print = list(row)
            for i in var_columns:
                to_print[i] *= 3
            print(*to_print, sep=',', file=f)
        pyfastaq.utils.close(f)


    @classmethod
    def _write_xls(cls, rows, outfile):
        workbook = openpyxl.Workbook()
        worksheet = workbook.worksheets[0]
        worksheet.title = 'ARIBA_summary'
        for row in rows:
            worksheet.append(row)
        workbook.save(outfile)


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


    @classmethod
    def _write_distance_matrix(cls, rows, outfile):
        if len(rows) < 3:
            raise Error('Cannot calculate distance matrix to make tree for phandango.\n' +
                        'Only one sample present.')

        if len(rows[0]) < 2:
            raise Error('Cannot calculate distance matrix to make tree for phandango.\n' +
                        'No genes present in output')

        with open(outfile, 'w') as f:
            sample_names = [x[0] for x in rows]
            print(*sample_names[1:], sep='\t', file=f)

            for i in range(1,len(rows)):
                scores = []
                for j in range(2, len(rows)):
                    scores.append(Summary._distance_score_between_lists(rows[i], rows[j]))
                print(rows[i][0], *scores, sep='\t', file=f)


    @classmethod
    def _newick_from_dist_matrix(cls, distance_file, outfile):
        r_script = outfile + '.tmp.R'

        with open(r_script, 'w') as f:
            print('library(ape)', file=f)
            print('a=read.table("', distance_file, '", header=TRUE, row.names=1, comment.char="")', sep='', file=f)
            print('h=hclust(dist(a))', file=f)
            print('write.tree(as.phylo(h), file="', outfile, '")', sep='', file=f)

        common.syscall('R CMD BATCH --no-save ' + r_script)
        if os.path.exists(r_script + 'out'):
            os.unlink(r_script + 'out')
        os.unlink(r_script)


    @classmethod
    def _write_phandango_files(cls, rows, outprefix):
        distance_file = outprefix + '.distance_matrix'
        tree_file = outprefix + '.tre'
        csv_file = outprefix + '.csv'
        Summary._write_distance_matrix(rows, distance_file)
        Summary._newick_from_dist_matrix(distance_file, tree_file)
        os.unlink(distance_file)
        Summary._write_phandango_csv(rows, csv_file)


    def run(self):
        self._check_files_exist()
        rows = Summary._gather_output_rows(self.filenames, self.min_id)

        if self.filter_output:
            rows = Summary._filter_output_rows(rows)

        if self.outfile.endswith('.xls'):
            Summary._write_xls(rows, self.outfile)
        else:
            Summary._write_tsv(rows, self.outfile)

        if self.phandango_prefix is not None:
            Summary._write_phandango_files(rows, self.phandango_prefix)
