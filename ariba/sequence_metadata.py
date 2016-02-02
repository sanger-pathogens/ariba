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

        self.variant_type = variant_type

        if self.variant_type == '.':
            self.variant = None
        else:
            self.variant = sequence_variant.Variant(self.variant_type, variant_string)

        self.hashed = hash((self.name, self.variant_type, variant_string))

    def __eq__(self, other):
       return type(other) is type(self) and self.__dict__ == other.__dict__


    def __lt__(self, other):
        return self.name < other.name or (self.name == other.name and self.variant < other.variant)


    def __hash__(self):
        return self.hashed


    def __str__(self):
        return self.to_string()


    def to_string(self, separator='\t'):
        fields = [self.name, self.variant_type]
        if self.variant is None:
            fields.append('.')
        else:
            fields.append(str(self.variant))

        if self.free_text:
            return separator.join(fields + [self.free_text])
        else:
            return separator.join(fields)


    def has_variant(self, seq):
        return self.variant is not None and self.variant.has_variant(seq)
