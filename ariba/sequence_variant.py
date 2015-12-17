import pyfastaq
import re

class Error (Exception): pass


allowed_variant_types = {'n', 'p'}

class Variant:
    def __init__(self, variant_type, variant_string):
        if variant_type not in allowed_variant_types:
            raise Error('Error! Variant type "' + variant_type + '" not recognised.\n' + \
                        'Must be one of:' + ', '.join(allowed_variant_types))

            self.variant_type = variant_type


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


    def sanity_check_against_seq(self, seq, translate_seq=False):
        if translate_seq:
            seq = pyfastaq.sequences.Fasta('x', seq).translate().seq

        return len(seq) >= self.position + 1 and seq[self.position].upper() in [self.wild_value, self.variant_value]
