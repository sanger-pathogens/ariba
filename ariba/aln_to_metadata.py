import os
import sys
import shutil
import pyfastaq
from ariba import sequence_variant

class Error (Exception): pass

class AlnToMetadata:
    def __init__(self,
      aln_file,
      vars_file,
      refs_are_coding
    ):
        self.padded_seqs = AlnToMetadata._load_aln_file(aln_file)
        self.variants = AlnToMetadata._load_vars_file(vars_file)
        self.refs_are_coding = refs_are_coding


    @classmethod
    def _load_aln_file(cls, aln_file):
        seqs = {}
        pyfastaq.tasks.file_to_dict(aln_file, seqs)
        return seqs


    @classmethod
    def _load_vars_file(cls, vars_file, refs_are_coding):
        var_type = 'p' if refs_are_coding else 'n'
        f = pyfastaq.utils.open_file_read(vars_file)
        variants = {}

        for line in f:
            try:
                ref_name, variant, identifier, description = line.rstrip().split('\t')
                variant = sequence_variant.Variant(var_type, variant, identifier)
            except:
                pyfastaq.utils.close(f)
                raise Error('Error in this line of variants file:\n' + line)

            if ref_name not in variants:
                variants[ref_name] = []

            variants[ref_name].append((variant, description))

        pyfastaq.utils.close(f)
        return variants


    @classmethod
    def _check_seq_lengths_same(cls, seqs):
        sequence_lengths = set([len(x) for x in seqs.values()])
        if len(sequence_lengths) > 1:
            raise Error('Input sequences must all be the same length. Cannot continue. Lengths found: ' + ','.join([str(x) for x in sequence_lengths]))
        return len(sequence_lengths) == 1


    def run(self, outprefix):
        pass
