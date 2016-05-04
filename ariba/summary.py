import os
import re
import sys
import openpyxl
import pyfastaq
from ariba import flag, common, report, summary_cluster, summary_sample

class Error (Exception): pass

required_keys_for_difference = {'no', 'yes', 'yes_nonunique', 'fragmented'}

class Summary:
    def __init__(
      self,
      outprefix,
      filenames=None,
      fofn=None,
      include_all_known_variant_columns=True,
      include_all_novel_variant_columns=False,
      min_id=90.0,
      cluster_cols='assembled,has_res,ref_seq,pct_id,known_var,novel_var',
      verbose=False,
    ):
        if filenames is None and fofn is None:
            raise Error('Error! Must supply filenames or fofn to Summary(). Cannot continue')

        if filenames is None:
            self.filenames = []
        else:
            self.filenames = filenames

        if fofn is not None:
            self.filenames.extend(self._load_fofn(fofn))

        self.cluster_columns = self._determine_cluster_cols(cluster_cols)
        self.include_all_known_variant_columns = include_all_known_variant_columns
        self.include_all_novel_variant_columns = include_all_novel_variant_columns
        self.min_id = min_id
        self.outprefix = outprefix
        self.verbose = verbose


    @staticmethod
    def _determine_cluster_cols(cols_string):
        allowed_cols = {'assembled', 'has_res', 'ref_seq', 'pct_id', 'known_var', 'novel_var'}
        if cols_string == '' or cols_string is None:
            return {x: False for x in allowed_cols}
        wanted_cols = set(cols_string.split(','))
        if not wanted_cols.issubset(allowed_cols):
            raise Error('Error in cluster names. Allowed values are: ' + str(','.join(list(allowed_cols))) + '. Got: ' + cols_string)
        return {x: x in wanted_cols for x in allowed_cols}


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
    def _load_input_files(cls, filenames, min_id, verbose=False):
        samples = {}
        for filename in filenames:
            samples[filename] = summary_sample.SummarySample(filename, min_pc_id=min_id)
            samples[filename].run()
            if verbose:
                print('Loaded file', filename, flush=True)
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
                        'has_res': 'NA',
                        'ref_seq': 'NA',
                        'known_var': 'NA',
                        'novel_var': 'NA',
                        'pct_id': 'NA'
                    }

                wanted_var_types = set()
                if self.include_all_known_variant_columns:
                    wanted_var_types.add('known')
                if self.include_all_novel_variant_columns:
                    wanted_var_types.add('unknown')

                if len(wanted_var_types) and cluster in all_var_columns:
                    for (ref_name, variant, known_or_unknown) in all_var_columns[cluster]:
                        if known_or_unknown not in wanted_var_types:
                            continue

                        key = ref_name + '.' + variant
                        if rows[filename][cluster]['assembled'] == 'no':
                            rows[filename][cluster][key] = 'NA'
                        elif cluster in sample.variant_column_names_tuples and (ref_name, variant, known_or_unknown) in sample.variant_column_names_tuples[cluster]:
                            rows[filename][cluster][key] = 'yes'
                        else:
                            rows[filename][cluster][key] = 'no'

                for key, wanted in self.cluster_columns.items():
                    if not wanted:
                        del rows[filename][cluster][key]

        return rows


    @classmethod
    def _to_matrix(cls, filenames, rows, cluster_cols):
        '''rows = output from _gather_output_rows().
           filenames = self.filenames
           cluster_cols = self.cluster_columns'''
        matrix = []
        making_header_lines = True
        phandango_header = ['name']
        phandago_suffixes = {'assembled': ':o1', 'has_res': ':o1', 'ref': ':o2', 'pct_id': ':c1', 'known_var': ':o1', 'novel_var': 'o1'}
        csv_header = ['name']
        all_cluster_cols_in_order = ['assembled', 'has_res', 'ref_seq', 'pct_id', 'known_var', 'novel_var']
        all_cluster_cols_in_order_set = set(['assembled', 'has_res', 'ref_seq', 'pct_id', 'known_var', 'novel_var'])
        cluster_cols_in_order = [x for x in all_cluster_cols_in_order if cluster_cols[x]]
        cluster_cols_set = set(cluster_cols_in_order)

        for filename in filenames:
            assert filename in rows
            line = [filename]

            for cluster_name in sorted(rows[filename]):
                for col in cluster_cols_in_order:
                    if making_header_lines:
                        csv_header.append(cluster_name + '.' + col)
                        phandango_header.append(cluster_name + '.' + col + '.' + phandago_suffixes[col])

                    line.append(rows[filename][cluster_name][col])

                for col in sorted(rows[filename][cluster_name]):
                    if col in all_cluster_cols_in_order_set:
                        continue

                    if making_header_lines:
                        csv_header.append(cluster_name + '.' + col)
                        phandango_header.append(cluster_name + '.' + col + ':o1')

                    line.append(rows[filename][cluster_name][col])

            making_header_lines = False
            matrix.append(line)

        return phandango_header, csv_header, matrix


    @classmethod
    def _filter_matrix_rows(cls, matrix):
        '''matrix = output from _to_matrix'''
        indexes_to_keep = []

        for i in range(len(matrix)):
            keep_row = False
            for element in matrix[i]:
                if element not in {'NA', 'no'}:
                    keep_row = True
                    break
            if keep_row:
                indexes_to_keep.append(i)

        return [matrix[i] for i in indexes_to_keep]


    @classmethod
    def _filter_clusters(cls, rows):
        '''Removes any column where every sample has "no" or "NA".
           Returns tuple: (filtered rows, number of remaining columns)'''
        columns_to_keep = set()
        first_filename = True
        all_clusters = set()

        for filename in rows:
            for cluster in rows[filename]:
                if first_filename:
                    all_clusters.add(cluster)

                if cluster in columns_to_keep:
                    continue

                if rows[filename][cluster]['assembled'] == 'yes':
                    columns_to_keep.add(cluster)

            first_filename = False

        to_delete = all_clusters.difference(columns_to_keep)

        for filename in rows:
            for cluster in to_delete:
                del rows[filename][cluster]

        return rows, len(columns_to_keep)


    @classmethod
    def _write_csv(cls, filenames, rows, outfile, phandango=False):
        lines = []
        non_var_keys_list = ['assembled', 'ref_seq', 'pct_id', 'known_var', 'novel_var']
        non_var_keys_set = set(non_var_keys_list)
        making_header_line = True
        first_line = ['name']

        # loop over filenames not rows to preserve their order
        for filename in filenames:
            assert filename in rows
            line = [filename]

            for cluster_name in sorted(rows[filename]):
                if making_header_line:
                    first_line.extend([
                        cluster_name + '.assembled',
                        cluster_name + '.ref',
                        cluster_name + '.pct_id',
                        cluster_name + '.known_var',
                        cluster_name + '.novel_var',
                    ])
                    if phandango:
                        first_line[-5] += ':o1'
                        first_line[-4] += ':o2'
                        first_line[-3] += ':c1'
                        first_line[-2] += ':o1'
                        first_line[-1] += ':o1'

                d = rows[filename][cluster_name]
                line.extend([d[x] for x in non_var_keys_list])

                for key, value in sorted(d.items()):
                    if key in non_var_keys_set:
                        continue

                    if making_header_line:
                        if phandango:
                            first_line.append(cluster_name + '.' + key + ':o1')
                        else:
                            first_line.append(cluster_name + '.' + key)

                    line.append(value)

            making_header_line = False
            lines.append(line)

        f = pyfastaq.utils.open_file_write(outfile)
        print(','.join(first_line), file=f)
        for line in lines:
            print(*line, sep=',', file=f)
        pyfastaq.utils.close(f)
        return [first_line] + lines


    @staticmethod
    def _distance_score_between_values(value1, value2):
        value_set = {value1, value2}
        if value_set.isdisjoint(required_keys_for_difference) or value1 == value2 or value_set == {'NA', 'no'}:
            return 0
        else:
            return 1


    @classmethod
    def _distance_score_between_lists(cls, scores1, scores2):
        assert len(scores1) == len(scores2)
        return sum([cls._distance_score_between_values(scores1[i], scores2[i]) for i in range(1, len(scores1))])


    @classmethod
    def _write_distance_matrix(cls, lines, outfile):
        if len(lines) < 3:
            raise Error('Cannot calculate distance matrix to make tree for phandango.\n' +
                        'Only one sample present.')

        if len(lines[0]) < 2:
            raise Error('Cannot calculate distance matrix to make tree for phandango.\n' +
                        'No genes present in output')

        with open(outfile, 'w') as f:
            sample_names = [x[0] for x in lines]
            print(*sample_names[1:], sep='\t', file=f)

            for i in range(1,len(lines)):
                scores = []
                for j in range(2, len(lines)):
                    scores.append(Summary._distance_score_between_lists(lines[i], lines[j]))
                print(lines[i][0], *scores, sep='\t', file=f)


    @classmethod
    def _newick_from_dist_matrix(cls, distance_file, outfile):
        r_script = outfile + '.tmp.R'

        with open(r_script, 'w') as f:
            print('library(ape)', file=f)
            print('a=read.table("', distance_file, '", header=TRUE, row.names=1, comment.char="")', sep='', file=f)
            print('h=hclust(dist(a))', file=f)
            print('write.tree(as.phylo(h), file="', outfile, '")', sep='', file=f)

        common.syscall('Rscript --no-save ' + r_script)
        if os.path.exists(r_script + 'out'):
            os.unlink(r_script + 'out')
        os.unlink(r_script)


    def run(self):
        if self.verbose:
            print('Loading input files...', flush=True)
        self._check_files_exist()
        self.samples = self._load_input_files(self.filenames, self.min_id, verbose=self.verbose)
        if self.verbose:
            print('Generating output rows', flush=True)
        self.rows = self._gather_output_rows()

        if self.verbose:
            print('Filtering columns', flush=True)
        self.rows, remaining_clusters = Summary._filter_clusters(self.rows)

        if remaining_clusters == 0:
            print('No clusters found that are present in any sample. Will not write any output files', file=sys.stderr)
            sys.exit(1)

        csv_file = self.outprefix + '.csv'
        if self.verbose:
            print('Writing csv file', csv_file, flush=True)
        Summary._write_csv(self.filenames, self.rows, csv_file, phandango=False)

        if len(self.samples) > 1:
            if self.verbose:
                print('Making Phandango csv file', csv_file, flush=True)
            csv_file = self.outprefix + '.phandango.csv'
            lines = Summary._write_csv(self.filenames, self.rows, csv_file, phandango=True)
            dist_matrix_file = self.outprefix + '.phandango.distance_matrix'
            tree_file = self.outprefix + '.phandango.tre'

            if self.verbose:
                print('Making Phandango distance matrix', dist_matrix_file, flush=True)
            Summary._write_distance_matrix(lines, dist_matrix_file)

            if self.verbose:
                print('Making Phandango tree file', tree_file, flush=True)
            Summary._newick_from_dist_matrix(dist_matrix_file, tree_file)
            os.unlink(dist_matrix_file)
        else:
            print('Made csv file. Not making Phandango files because only one input file given', file=sys.stderr)

        if self.verbose:
            print('Finished', flush=True)
