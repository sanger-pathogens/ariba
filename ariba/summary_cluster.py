from ariba import flag, report

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

class SummaryCluster:
    def __init__(self, min_pc_id=90):
        self.min_pc_id = min_pc_id
        self.name = None
        self.ref_name = None
        self.flag = None
        self.data = []


    def __eq__(self, other):
       return type(other) is type(self) and self.__dict__ == other.__dict__


    @classmethod
    def line2dict(cls, line):
        data = line.rstrip().split('\t')
        if len(data) != len(report.columns):
            raise Error('Wrong number of columns in the following line. Expected ' + str(len(report.columns)) + ' but got ' + str(len(data)) + '\n' + line)
        d = {report.columns[i]: data[i] for i in range(len(data))}
        try:
            d['flag'] = flag.Flag(int(d['flag']) )
        except:
            raise Error('Error getting flag in the following line. Got "' + d['flag'] + '" for the flag.\n' + line)

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


    def add_data_dict(self, data_dict):
        if data_dict['pc_ident'] == '.' or data_dict['pc_ident'] < self.min_pc_id:
            return

        if self.name is None:
            assert self.ref_name is None and self.flag is None
            self.name = data_dict['cluster']
            self.ref_name = data_dict['ref_name']
            self.flag = data_dict['flag']

        if self.name != data_dict['cluster']:
            raise Error('Cannot add dict to SummaryCluster. Expected cluster name "' + self.name + '" but got "' + data_dict['cluster'] + '".')

        if self.ref_name != data_dict['ref_name']:
            raise Error('Cannot add dict to SummaryCluster. Expected ref_name "' + self.ref_name + '" but got "' + data_dict['ref_name'] + '".')

        if self.flag != data_dict['flag']:
            raise Error('Cannot add dict to SummaryCluster. Expected flag "' + str(self.flag) + '" but got "' + str(data_dict['flag']) + '".')
        self.data.append(data_dict)


    def pc_id_of_longest(self):
        longest = 0
        identity = 0

        for d in self.data:
            if d['ref_base_assembled'] > longest:
                longest = d['ref_base_assembled']
                identity = d['pc_ident']

        return identity


    def _to_cluster_summary_assembled(self):
        if self.data[0]['ref_type'] == 'non_coding':
            has_complete_orf = True
        else:
            has_complete_orf = self.flag.has('complete_orf')

        if self.flag.has('assembly_fail') or \
          (not self.flag.has('assembled')) or \
          self.flag.has('ref_seq_choose_fail'):
            return 'no'
        elif self.flag.has('assembled_into_one_contig') and has_complete_orf:
            if self.flag.has('unique_contig') and \
              (not self.flag.has('scaffold_graph_bad')) and \
              (not self.flag.has('variants_suggest_collapsed_repeat')) and \
              (not self.flag.has('hit_both_strands')) and \
              (not self.flag.has('region_assembled_twice')):
                return 'yes'
            else:
                return 'yes_nonunique'
        else:
            return 'fragmented'


    @classmethod
    def _has_nonsynonymous(cls, data_dict):
        return data_dict['ref_ctg_effect'] != 'SYN' and \
          (data_dict['has_known_var'] == '1' or (data_dict['known_var'] != '1' and data_dict['ref_ctg_change'] != '.'))


    def _has_any_nonsynonymous(self):
        for d in self.data:
            if self._has_nonsynonymous(d):
                return 'yes'
        return 'no'


    def _to_cluster_summary_has_nonsynonymous(self, assembled_summary):
        '''assembled_summary should be output of _to_cluster_summary_assembled'''
        if assembled_summary == 'no':
            return 'NA'
        else:
            return self._has_any_nonsynonymous()


    @staticmethod
    def _get_nonsynonymous_var(data_dict):
        '''if data_dict has a non synonymous variant, return string:
        ref_name.change. Otherwise return None'''
        has_nonsyn = SummaryCluster._has_nonsynonymous(data_dict)

        if not has_nonsyn:
            return None
        elif data_dict['known_var_change'] == data_dict['ref_ctg_change'] == '.':
            raise Error('Unexpected data in ariba summary... \n' + str(data_dict) + '\n... known_var_change and ref_ctg_change both equal to ".", but var_type was not a ".". Cannot continue')
        else:
            if '.' not in [data_dict['known_var_change'], data_dict['ref_ctg_change']] and \
              data_dict['known_var_change'] != data_dict['ref_ctg_change']:
                raise Error('Unexpected data in ariba summary... \n' + str(data_dict) + '\n... known_var_change != ref_ctg_change. Cannot continue')

            if data_dict['known_var_change'] != '.':
                return data_dict['known_var_change']
            else:
                return data_dict['ref_ctg_change']


    def column_summary_data(self):
        '''Returns a dictionary of column name -> value, for cluster-level results'''
        assembled_summary = self._to_cluster_summary_assembled()

        columns = {
            'assembled': self._to_cluster_summary_assembled(),
            'ref_seq': self.ref_name,
            'pct_id': str(self.pc_id_of_longest()),
            'any_var': self._to_cluster_summary_has_nonsynonymous(assembled_summary)
        }

        return columns


    def non_synon_variants(self):
        variants = {self._get_nonsynonymous_var(d) for d in self.data}
        variants.discard(None)
        return variants
