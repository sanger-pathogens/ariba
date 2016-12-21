import sys
import csv
import os

class Error (Exception): pass

class MlstProfile:
    def __init__(self, infile, duplicate_warnings=True):
        self.infile = infile
        self.duplicate_warnings = duplicate_warnings
        self.columns_to_ignore = ['clonal_complex', 'CC', 'Lineage', 'mlst_clade', 'species']

        if not os.path.exists(self.infile):
            raise Error('Error! Input file "' + self.infile + '" not found.')
        self._load_input_file()


    def _load_input_file(self):
        self.profile_to_type = {}

        with open(self.infile) as f:
            reader = csv.DictReader(f, delimiter='\t')
            if reader.fieldnames[0] != 'ST':
                raise Error('Error. Expected first column of profile file "' + self.infile + '" to be "ST"')

            self.genes_list = [column_name for column_name in reader.fieldnames[1:] if column_name not in self.columns_to_ignore]

            for row in reader:
                type_tuple = tuple(int(row[x]) for x in self.genes_list)
                if type_tuple in self.profile_to_type:
                    previous_st = self.profile_to_type[type_tuple]
                    this_st = int(row['ST'])
                    if previous_st != this_st:
                        new_st = min(previous_st, this_st)
                        if self.duplicate_warnings:
                            print('WARNING: Same profile found twice in input file, but two different STs. Going to use the ST with the smaller number (', new_st, ')', sep='', file=sys.stderr)
                            print(' ... STs are', previous_st, this_st, 'and alleles are', ', '.join([x + ':' + row[x] for x in self.genes_list]), file=sys.stderr)
                        self.profile_to_type[type_tuple] = new_st
                else:
                    self.profile_to_type[type_tuple] = int(row['ST'])

        self.genes_set = set(self.genes_list)


    def has_gene(self, gene):
        return gene in self.genes_set


    def get_sequence_type(self, type_dict):
        key = tuple(type_dict.get(x, 'ND') for x in self.genes_list)
        if 'ND' in key:
            return 'ND'
        else:
            return self.profile_to_type.get(key, 'Novel')


