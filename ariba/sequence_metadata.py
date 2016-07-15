from ariba import sequence_variant

class Error (Exception): pass


class SequenceMetadata:
    def __init__(self, line):
        try:
            self.name, seq_type, var_only, variant, variant_id, self.free_text = line.rstrip().split('\t')
        except:
            raise Error('Error parsing line of file:\n' + line)

        if seq_type not in {'0', '1'}:
            raise Error('Error. Second column must be "0" or "1". Cannot continue. Line was:\n' + line)

        self.seq_type = 'n' if seq_type == '0' else 'p'

        if var_only not in {'0', '1'}:
            raise Error('Error. Third column must be "0" or "1". Cannot continue. Line was:\n' + line)

        self.variant_only = var_only == '1'

        if variant == '.':
            self.variant = None
        else:
            self.variant = sequence_variant.Variant(self.seq_type, variant, variant_id)


    def __eq__(self, other):
       return type(other) is type(self) and self.name == other.name and self.seq_type == other.seq_type and self.variant_only == other.variant_only and self.variant == other.variant and self.free_text == other.free_text


    def __lt__(self, other):
        return self.name < other.name or (self.name == other.name and self.variant is not None and other.variant is not None and self.variant < other.variant)


    def __hash__(self):
        return hash((self.name, self.seq_type, str(self.variant), self.free_text))


    def __str__(self):
        return self.to_string()


    def to_string(self, separator='\t'):
        return separator.join([
            self.name,
            '1' if self.seq_type == 'p' else '0',
            '1' if self.variant_only else '0',
            '.' if self.variant is None else str(self.variant),
            '.' if (self.variant is None or self.variant.identifier is None) else self.variant.identifier,
            self.free_text
        ])


    def has_variant(self, seq):
        return self.variant is not None and self.variant.has_variant(seq)
