import os
import copy
import sys
import pyfastaq
import dendropy
from ariba import summary_sample

class Error (Exception): pass

required_keys_for_difference = {'no', 'yes', 'yes_nonunique', 'fragmented'}

class Summary:
    def __init__(
      self,
      outprefix,
      filenames=None,
      fofn=None,
      filter_rows=True,
      filter_columns=True,
      min_id=90.0,
      cluster_cols='assembled,match,ref_seq,pct_id,known_var,novel_var',
      make_phandango_tree=True,
      only_clusters=None,
      show_var_groups=False,
      show_known_vars=False,
      show_novel_vars=False,
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
        self.filter_rows = filter_rows
        self.filter_columns = filter_columns
        self.min_id = min_id
        self.outprefix = outprefix
        self.make_phandango_tree = make_phandango_tree
        self.only_clusters = only_clusters
        self.show_var_groups = show_var_groups
        self.show_known_vars = show_known_vars
        self.show_novel_vars = show_novel_vars
        self.verbose = verbose


    @classmethod
    def _determine_cols(cls, cols_string, allowed_cols, error_string):
        if cols_string == '' or cols_string is None:
            return {x: False for x in allowed_cols}
        wanted_cols = set(cols_string.split(','))
        if not wanted_cols.issubset(allowed_cols):
            raise Error('Error in ' + error_string + '. Allowed values are: ' + str(','.join(list(allowed_cols))) + '. Got: ' + cols_string)
        return {x: x in wanted_cols for x in allowed_cols}


    @staticmethod
    def _determine_cluster_cols(cols_string):
        allowed_cols = {'assembled', 'match', 'ref_seq', 'pct_id', 'known_var', 'novel_var'}
        return Summary._determine_cols(cols_string, allowed_cols, 'cluster columns')


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
    def _load_input_files(cls, filenames, min_id, verbose=False, only_clusters=None):
        samples = {}
        for filename in filenames:
            samples[filename] = summary_sample.SummarySample(filename, min_pc_id=min_id, only_clusters=only_clusters)
            samples[filename].run()
            if verbose:
                print('Loaded file', filename, flush=True)
        return samples


    def _gather_unfiltered_output_data(self):
        self.all_potential_columns = {}
        self.all_data = {}

        for filename in sorted(self.samples):
            self.all_data[filename] = {}
            for cluster in self.samples[filename].clusters.values():
                self.all_data[filename][cluster.name] = {}
                if cluster.name not in self.all_potential_columns:
                    self.all_potential_columns[cluster.name] = {'summary' : set(), 'groups': set(), 'vars': set()}

                this_cluster_dict = {'groups': {}, 'vars': {}}

                if cluster.summary['assembled'] == 'no':
                    this_cluster_dict['summary'] = {
                            'assembled': 'no',
                            'known_var': 'NA',
                            'match': 'no',
                            'novel_var': 'NA',
                            'pct_id': 'NA',
                            'ref_seq': 'NA'
                    }
                else:
                    this_cluster_dict['summary'] = copy.copy(cluster.summary)
                    seen_groups = {}

                    for variant in cluster.variants:
                        if (self.show_known_vars and variant.known) or (self.show_novel_vars and not variant.known):
                            this_cluster_dict['vars'][variant.var_string] = 'yes' if (variant.het_percent is None or not variant.is_het) else 'het'
                            if variant.het_percent is not None:
                                this_cluster_dict['vars'][variant.var_string + '.%'] = variant.het_percent

                        if self.show_var_groups and variant.var_group != '.':
                            if variant.var_group not in seen_groups:
                                seen_groups[variant.var_group] = {'yes': 0, 'het': 0}

                            if variant.het_percent is not None:
                                this_cluster_dict['groups'][variant.var_group + '.%'] = variant.het_percent

                            if variant.is_het:
                                seen_groups[variant.var_group]['het'] += 1
                                this_cluster_dict['groups'][variant.var_group] = 'het'
                            else:
                                seen_groups[variant.var_group]['yes'] += 1
                                this_cluster_dict['groups'][variant.var_group] = 'yes'

                    for group, d in seen_groups.items():
                        if d['het'] > 0 and d['het'] + d['yes'] > 1:
                            this_cluster_dict['groups'][group] = 'yes_multi_het'
                            this_cluster_dict['groups'][group + '.%'] = 'NA'

                for x in this_cluster_dict:
                    self.all_potential_columns[cluster.name][x].update(set(this_cluster_dict[x].keys()))

                self.all_data[filename][cluster.name] = this_cluster_dict


    @classmethod
    def _to_matrix(cls, filenames, all_data, all_potential_columns, cluster_cols):
        matrix = []
        making_header_lines = True
        phandango_header = ['name']
        phandango_suffixes = {'assembled': ':o1', 'match': ':o1', 'ref_seq': ':o2', 'pct_id': ':c1', 'known_var': ':o1', 'novel_var': ':o1'}
        ref_seq_counter = 2
        csv_header = ['name']
        summary_cols_in_order = ['assembled', 'match', 'ref_seq', 'pct_id', 'known_var', 'novel_var']
        summary_cols_set = set(['assembled', 'match', 'ref_seq', 'pct_id', 'known_var', 'novel_var'])
        summary_cols_in_order = [x for x in summary_cols_in_order if cluster_cols[x]]

        for filename in filenames:
            line = [filename]

            for cluster_name in sorted(all_potential_columns):
                group_cols = sorted(list(all_potential_columns[cluster_name]['groups']))
                var_cols = sorted(list(all_potential_columns[cluster_name]['vars']))

                for col in summary_cols_in_order + group_cols + var_cols:
                    if making_header_lines:
                        csv_header.append(cluster_name + '.' + col)
                        if col == 'ref_seq':
                            phandango_suffixes[col] = ':o' + str(ref_seq_counter)
                            ref_seq_counter += 1
                            phandango_header.append(cluster_name + '.' + col + phandango_suffixes[col])
                        elif col in phandango_suffixes:
                            phandango_header.append(cluster_name + '.' + col + phandango_suffixes[col])
                        elif col.endswith('.%'):
                            phandango_header.append(cluster_name + '.' + col + ':c2')
                        else:
                            phandango_header.append(cluster_name + '.' + col + ':o1')

                    for col_type in ['summary', 'groups', 'vars']:
                        if cluster_name in all_data[filename] and col in all_data[filename][cluster_name][col_type]:
                            line.append(all_data[filename][cluster_name][col_type][col])
                            break
                    else:
                        if col in {'assembled', 'match'}:
                            line.append('no')
                        elif col in summary_cols_set:
                            line.append('NA')
                        elif cluster_name in all_data[filename] and all_data[filename][cluster_name]['summary'].get('assembled', 'no')  != 'no':
                            if col.endswith('.%'):
                                line.append('NA')
                            else:
                                line.append('no')
                        else:
                            line.append('NA')

            making_header_lines = False
            matrix.append(line)

        assert len(phandango_header) == len(csv_header)
        for line in matrix:
            assert len(line) == len(csv_header)
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
    def _filter_matrix_columns(cls, matrix, phandango_header, csv_header):
        '''phandango_header, csv_header, matrix = output from _to_matrix'''
        indexes_to_keep = set()

        for row in matrix:
            for i in range(len(row)):
                if row[i] not in {'NA', 'no'}:
                    indexes_to_keep.add(i)

        indexes_to_keep = sorted(list(indexes_to_keep))

        for i in range(len(matrix)):
            matrix[i] = [matrix[i][j] for j in indexes_to_keep]

        phandango_header = [phandango_header[i] for i in indexes_to_keep]
        csv_header = [csv_header[i] for i in indexes_to_keep]
        return phandango_header, csv_header, matrix


    @classmethod
    def _add_phandango_colour_columns(cls, header, matrix):
        header = copy.deepcopy(header)
        matrix = copy.deepcopy(matrix)
        cols_to_add_colour_col = [i for i in range(len(header)) if header[i].endswith(':o1')]
        field_to_col = {
            'yes': '#33a02c',
            'yes_nonunique': '#b2df8a',
            'no': '#fb9a99',
            'NA': '#ffffff',
            'het': '#fdbf6f',
            'fragmented': '#1f78b4',
            'interrupted': '#a6cee3',
            'partial': '#fdbf6f',
        }

        cols_to_add_colour_col.reverse()

        for col_index in cols_to_add_colour_col:
            header[col_index] = header[col_index][:-3]
            header.insert(col_index + 1, header[col_index] + ':colour')

            for row_index in range(len(matrix)):
                colour = field_to_col[matrix[row_index][col_index]]
                matrix[row_index].insert(col_index + 1, colour)

        return header, matrix


    @classmethod
    def _matrix_to_csv(cls, matrix, header, outfile, remove_nas=False):
        f = pyfastaq.utils.open_file_write(outfile)
        fixed_header = [x.replace(',', '/') for x in header]
        print(*fixed_header, sep=',', file=f)
        for line in matrix:
            if remove_nas:
                new_line = ['' if x=='NA' else x for x in line]
                print(*new_line, sep=',', file=f)
            else:
                print(*line, sep=',', file=f)
        pyfastaq.utils.close(f)


    @staticmethod
    def _distance_score_between_values(value1, value2):
        if value1 == 'partial':
            value1 = 'no'
        if value2 == 'partial':
            value2 = 'no'
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
        if len(lines) < 2:
            raise Error('Cannot calculate distance matrix to make tree for phandango.\n' +
                        'Only one sample present.')

        if len(lines[0]) < 2:
            raise Error('Cannot calculate distance matrix to make tree for phandango. Not enough columns')

        scores = [[0 for i in range(len(lines))] for j in range(len(lines))]

        for i in range(len(lines)):
            for j in range(i + 1, len(lines), 1):
                scores[i][j] = Summary._distance_score_between_lists(lines[i], lines[j])
                scores[j][i] = scores[i][j]

        with open(outfile, 'w') as f:
            sample_names = [''] + [x[0] for x in lines]
            print(*sample_names, sep='\t', file=f)
            for i in range(len(scores)):
                print(lines[i][0], *scores[i], sep='\t', file=f)


    @classmethod
    def _newick_from_dist_matrix(cls, distance_file, outfile):
        with open(distance_file) as f:
            pdm = dendropy.PhylogeneticDistanceMatrix.from_csv(src=f, delimiter='\t')
        upgma_tree = pdm.upgma_tree()
        with open(outfile, 'w') as f:
            print(upgma_tree.as_string("newick").replace("'", ''), end='', file=f)


    def run(self):
        if self.verbose:
            print('Loading input files...', flush=True)
        self._check_files_exist()
        self.samples = self._load_input_files(self.filenames, self.min_id, verbose=self.verbose, only_clusters=self.only_clusters)
        if self.verbose:
            print('Generating output rows', flush=True)
        self._gather_unfiltered_output_data()
        phandango_header, csv_header, matrix = Summary._to_matrix(self.filenames, self.all_data, self.all_potential_columns, self.cluster_columns)

        # sanity check same number of columns in headers and matrix
        lengths = {len(x) for x in matrix}
        assert len(lengths) == 1
        assert len(matrix[0]) == len(phandango_header) == len(csv_header)

        if self.filter_rows:
            if self.verbose:
                print('Filtering rows', flush=True)
            matrix = Summary._filter_matrix_rows(matrix)

        if len(matrix) == 0:
            print('No rows left after filtering rows. Cannot continue', file=sys.stderr)
            sys.exit(1)

        if self.filter_columns:
            if self.verbose:
                print('Filtering columns', flush=True)
            phandango_header, csv_header, matrix = Summary._filter_matrix_columns(matrix, phandango_header, csv_header)

        if len(matrix) == 0 or len(matrix[0]) == 0:
            print('No columns left after filtering columns. Cannot continue', file=sys.stderr)

        # sanity check same number of columns in headers and matrix
        lengths = {len(x) for x in matrix}
        assert len(lengths) == 1
        assert len(matrix[0]) == len(phandango_header) == len(csv_header)

        csv_file = self.outprefix + '.csv'
        if self.verbose:
            print('Writing csv file', csv_file, flush=True)
        Summary._matrix_to_csv(matrix, csv_header, csv_file)

        if len(matrix) > 1:
            if self.verbose:
                print('Making Phandango csv file', csv_file, flush=True)
            csv_file = self.outprefix + '.phandango.csv'
            phandango_header, phandango_matrix = Summary._add_phandango_colour_columns(phandango_header, matrix)
            Summary._matrix_to_csv(phandango_matrix, phandango_header, csv_file, remove_nas=True)

            if self.make_phandango_tree:
                dist_matrix_file = self.outprefix + '.phandango.distance_matrix'
                tree_file = self.outprefix + '.phandango.tre'

                if self.verbose:
                    print('Making Phandango distance matrix', dist_matrix_file, flush=True)
                Summary._write_distance_matrix(matrix, dist_matrix_file)

                if self.verbose:
                    print('Making Phandango tree file', tree_file, flush=True)
                Summary._newick_from_dist_matrix(dist_matrix_file, tree_file)
                os.unlink(dist_matrix_file)
            elif self.verbose:
                print('Skipping making tree because you asked me not to make it', flush=True)
        else:
            print('Made csv file. Not making Phandango files because only one sample remains after filtering', file=sys.stderr)

        if self.verbose:
            print('Finished', flush=True)

