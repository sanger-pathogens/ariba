import re

class Error (Exception): pass


NUCLEOTIDE_SNP = 1
AMINO_ACID_SNP = 2

variant_type_strings = {
    'n': NUCLEOTIDE_SNP,
    'p': AMINO_ACID_SNP,
}

class Variant:
    def __init__(self, variant_type, variant_string):
        try:
            self.variant_type = variant_type_strings[variant_type]
        except:
            raise Error('Error! Variant type "' + variant_type + '" not recognised.\n' + \
                        'Must be one of:' + ', '.join(variant_type_strings.keys()))


        m = re.match('^([A-Z])([0-9]+)([A-Z])$', variant_string.upper())
        if m is None:
            raise Error('Unexpected format of variant string: ', variant_string)

        try:
            self.wild_value, self.position, self.variant_value = m.group(1, 2, 3)
        except:
            raise Error('Error getting amino acids and position of variant from', variant_string)

        self.position = int(self.position) - 1


    def __str__(self):
        return ''.join([self.wild_value, str(self.position + 1), self.variant_value])


    def sanity_check_against_seq(self, seq):
        return len(seq) >= self.position + 1 and seq[self.position].upper() in [self.wild_value, self.variant_value]
