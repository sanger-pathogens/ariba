import pyfastaq
import re

class Error (Exception): pass


allowed_variant_types = {'n', 'p'}

class Variant:
    def __init__(self, variant_type, variant_string, identifier):
        if variant_type not in allowed_variant_types:
            raise Error('Error! Variant type "' + variant_type + '" not recognised.\n' + \
                        'Must be one of:' + ', '.join(allowed_variant_types))

        self.variant_type = variant_type
        self.identifier = None if identifier == '.' else identifier


        m = re.match('^([A-Z])([0-9]+)([A-Z])$', variant_string.upper())
        if m is None:
            raise Error('Unexpected format of variant string: ', variant_string)

        try:
            self.wild_value, self.position, self.variant_value = m.group(1, 2, 3)
        except:
            raise Error('Error getting amino acids and position of variant from', variant_string)

        self.position = int(self.position) - 1


    def __eq__(self, other):
       return type(other) is type(self) and self.__dict__ == other.__dict__


    def __lt__(self, other):
        return self.position < other.position or \
            (self.position == other.position and self.variant_type < other.variant_type) or \
            (self.position == other.position and self.variant_type == other.variant_type and self.wild_value < other.wild_value) or \
            (self.position == other.position and self.variant_type == other.variant_type and self.wild_value == other.wild_value and self.variant_value < other.variant_value)


    def __str__(self):
        return ''.join([self.wild_value, str(self.position + 1), self.variant_value])


    def sanity_check_against_seq(self, seq, translate_seq=False):
        if translate_seq:
            seq = pyfastaq.sequences.Fasta('x', seq).translate().seq

        return len(seq) >= self.position + 1 and seq[self.position].upper() in [self.wild_value, self.variant_value]


    def has_variant(self, seq):
        if self.variant_type == 'p':
            test_seq = seq.translate()
        else:
            test_seq = seq

        assert self.position < len(test_seq)
        return test_seq[self.position] == self.variant_value


    def nucleotide_range(self):
        '''Returns the nucleotide (start, end) positions inclusive of this variant.
           start==end if it's an amino acid variant, otherwise start+2==end'''
        if self.variant_type == 'p':
            return 3 * self.position, 3 * self.position + 2
        else:
            return self.position, self.position

