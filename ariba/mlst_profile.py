import csv
import os

class Error (Exception): pass

class MlstProfile:
    def __init__(self, infile):
        self.infile = infile
        if not os.path.exists(self.infile):
            raise Error('Error! Input file "' + self.infile + '" not found.')
        self._load_input_file()


    def _load_input_file(self):
        self.profile_to_type = {}

        with open(self.infile) as f:
            reader = csv.DictReader(f, delimiter='\t')
            if reader.fieldnames[0] != 'ST':
                raise Error('Error. Expected first column of profile file "' + self.infile + '" to be "ST"')

            self.genes_list = reader.fieldnames[1:]
            if self.genes_list[-1] == 'clonal_complex':
                self.genes_list.pop()

            for row in reader:
                type_tuple = tuple(int(row[x]) for x in self.genes_list)
                assert type_tuple not in self.profile_to_type
                self.profile_to_type[type_tuple] = int(row['ST'])

        self.genes_set = set(self.genes_list)


    def has_gene(self, gene):
        return gene in self.genes_set


    def get_sequence_type(self, type_dict):
        key = tuple(type_dict.get(x, 'ND') for x in self.genes_list)
        return self.profile_to_type.get(key, 'ND')

