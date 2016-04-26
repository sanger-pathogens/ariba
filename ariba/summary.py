import os
import re
import openpyxl
import pyfastaq
from ariba import flag, common, report, summary_cluster, summary_sample

class Error (Exception): pass

class Summary:
    def __init__(
      self,
      outfile,
      filenames=None,
      fofn=None,
      include_all_variant_columns=False,
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

        self.include_all_variant_columns = include_all_variant_columns
        self.min_id = min_id
        self.outfile = outfile
        self.phandango_prefix = phandango_prefix


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
    def _load_input_files(cls, filenames, min_id):
        samples = {}
        for filename in filenames:
            samples[filename] = summary_sample.SummarySample(filename, min_pc_id=min_id)
            samples[filename].run()
        return samples


    @classmethod
    def _get_all_cluster_names(cls, samples_dict):
        '''Input should be output of _load_input_files'''
        cluster_names = set()
        for filename, sample in samples_dict.items():
            cluster_names.update(set(sample.clusters.keys()))
        return cluster_names


    @classmethod
    def _get_all_variant_columns(cls, samples_dict):
        '''Input should be output of _load_input_files'''
        columns = {}
        for filename, sample in samples_dict.items():
            for cluster in sample.column_summary_data:
                if sample.column_summary_data[cluster]['assembled'] == 'yes':
                    for key, tuple_set in sample.variant_column_names_tuples.items():
                        for t in tuple_set:
                            if key not in columns:
                                columns[key] = set()
                            columns[key].add(t)
        return columns


    def _gather_output_rows(self):
        all_cluster_names = Summary._get_all_cluster_names(self.samples)
        all_var_columns = Summary._get_all_variant_columns(self.samples)
        rows = {}

        for filename, sample in self.samples.items():
            rows[filename] = {}

            for cluster in all_cluster_names:
                rows[filename][cluster] = {}

                if cluster in sample.column_summary_data and sample.column_summary_data[cluster]['assembled'].startswith('yes'):
                    rows[filename][cluster] = sample.column_summary_data[cluster]
                else:
                    rows[filename][cluster] = {
                        'assembled': 'no',
                        'ref_seq': 'NA',
                        'any_var': 'NA',
                        'pct_id': 'NA'
                    }

                if self.include_all_variant_columns and cluster in all_var_columns:
                    for (ref_name, variant) in all_var_columns[cluster]:
                        key = ref_name + '.' + variant
                        if rows[filename][cluster]['assembled'] == 'no':
                            rows[filename][cluster][key] = 'NA'
                        elif cluster in sample.variant_column_names_tuples and (ref_name, variant) in sample.variant_column_names_tuples[cluster]:
                            rows[filename][cluster][key] = 'yes'
                        else:
                            rows[filename][cluster][key] = 'no'

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
        self.samples = self._load_input_files(self.filenames, self.min_id)

        if self.filter_output:
            rows = Summary._filter_output_rows(rows)

        if self.outfile.endswith('.xls'):
            Summary._write_xls(rows, self.outfile)
        else:
            Summary._write_tsv(rows, self.outfile)

        if self.phandango_prefix is not None:
            Summary._write_phandango_files(rows, self.phandango_prefix)
