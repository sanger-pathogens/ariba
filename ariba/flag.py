class Error (Exception): pass


flags_in_order = [
    'assembled',
    'assembled_into_one_contig',
    'region_assembled_twice',
    'complete_gene',
    'unique_contig',
    'scaffold_graph_bad',
    'assembly_fail',
    'variants_suggest_collapsed_repeat',
    'hit_both_strands',
    'has_variant',
    'ref_seq_choose_fail',
]


flag_bits = {flags_in_order[i]: pow(2, i) for i in range(len(flags_in_order))}


class Flag:
    def __init__(self, n=0):
        self.flags = {x: False for x in flags_in_order}
        self.set_flag(n)


    def set_flag(self, n):
        for f in self.flags:
            if flag_bits[f] & n != 0:
                self.flags[f] = True


    def add(self, f):
        self.flags[f] = True


    def to_number(self):
        n = 0
        for f in self.flags:
            if self.flags[f]:
                n += flag_bits[f]
        return n

    def __eq__(self, other):
       return type(other) is type(self) and self.__dict__ == other.__dict__


    def __str__(self):
        return str(self.to_number())


    def to_long_string(self):
        lines = []
        for f in flags_in_order:
            x_or_not = 'X' if self.flags[f] else ' '
            lines.append('[' + x_or_not + '] ' + f)
        return '\n'.join(lines)


    def has(self, s):
        return self.flags[s]

