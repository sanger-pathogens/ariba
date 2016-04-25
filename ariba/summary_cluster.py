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
        if data_dict['pc_ident'] < self.min_pc_id:
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

