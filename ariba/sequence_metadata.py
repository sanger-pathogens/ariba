from ariba import sequence_variant

class Error (Exception): pass


class SequenceMetadata:
    def __init__(self, line):
        try:
            self.name, variant_type, variant_string, identifier, self.free_text = line.rstrip().split('\t')
        except:
            raise Error('Error parsing line of file:\n' + line)

        self.variant_type = variant_type

        if self.variant_type == '.':
            self.variant = None
        else:
            self.variant = sequence_variant.Variant(self.variant_type, variant_string, identifier)


    def __eq__(self, other):
       return type(other) is type(self) and self.name == other.name and self.variant == other.variant and self.variant_type == other.variant_type and self.free_text == other.free_text


    def __lt__(self, other):
        return self.name < other.name or (self.name == other.name and self.variant < other.variant)


    def __hash__(self):
        return hash((self.name, self.variant_type, str(self.variant), self.free_text))


    def __str__(self):
        return self.to_string()


    def to_string(self, separator='\t'):
        return separator.join([
            self.name,
            self.variant_type,
            '.' if self.variant is None else str(self.variant),
            '.' if (self.variant is None or self.variant.identifier is None) else self.variant.identifier,
            self.free_text
        ])


    def has_variant(self, seq):
        return self.variant is not None and self.variant.has_variant(seq)
