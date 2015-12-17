from ariba import sequence_variant

class Error (Exception): pass


class SequenceMetadata:
    def __init__(self, line):
        try:
            self.name, variant_type, variant_string, *extra_columns = line.rstrip().split('\t')
        except:
            raise Error('Error parsing line of file:\n' + line)

        if len(extra_columns) == 0:
            self.free_text = None
        elif len(extra_columns) == 1:
            self.free_text = extra_columns[0]
        else:
            raise Error('Too many columns in this line:\n' + line)

        if variant_type == '.':
            self.variant = None
        else:
            self.variant = sequence_variant.Variant(variant_type, variant_string)

