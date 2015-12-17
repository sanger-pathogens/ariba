import re

class Error (Exception): pass


class Variant:
    def __init__(self, variant_string):
        m = re.match('^([A-Z])([0-9]+)([A-Z])$', variant_string.upper())
        if m is None:
            raise Error('Unexpected format of variant string: ', variant_string)

        try:
            self.wild_aa, self.position, self.variant_aa = m.group(1, 2, 3)
        except:
            raise Error('Error getting amino acids and position of variant from', variant_string)

        self.position = int(self.position) - 1


    def __str__(self):
        return ''.join([self.wild_aa, str(self.position + 1), self.variant_aa])


    def agrees_with_protein_seq(self, protein_seq):
        return len(protein_seq) >= self.position + 1 and protein_seq[self.position].upper() in [self.wild_aa, self.variant_aa]
