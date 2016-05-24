import os
import re
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
    def _make_unpadded_seqs(cls, padded_seqs):
        unpadded_seqs = {}
        for seq in padded_seqs.values():
            unpadded_seqs[seq.id] = pyfastaq.sequences.Fasta(seq.id, seq.seq.replace('-', ''))
        return unpadded_seqs


    @classmethod
    def _check_seq_lengths_same(cls, seqs):
        sequence_lengths = set([len(x) for x in seqs.values()])
        if len(sequence_lengths) > 1:
            raise Error('Input sequences must all be the same length. Cannot continue. Lengths found: ' + ','.join([str(x) for x in sequence_lengths]))
        return len(sequence_lengths) == 1


    @classmethod
    def _insertion_coords(cls, sequence):
        insertions = []
        regex = re.compile('-+')
        for m in regex.finditer(sequence.seq):
             insertions.append(pyfastaq.intervals.Interval(m.span()[0], m.span()[1] - 1))
        return insertions


    @classmethod
    def _check_insertion_coords(cls, sequence):
        insertions = AlnToMetadata._insertion_coords(sequence)
        for coords in insertions:
            if coords.start % 3 !=0:
                raise Error('Insertion does not start in frame in sequence "' + sequence.id + '". Cannot continue')
            elif len(coords) % 3 != 0:
                raise Error('Insertion of length not a mulitple of 3 in sequence "' + sequence.id + '". Cannot continue')

        return True


    @classmethod
    def _check_coding_seq(cls, sequence, genetic_code=11):
        if len(sequence) % 3 != 0:
            raise Error('Length of sequence ' + sequence.id + ' is ' + str(len(sequence)) + ', which is not a multiple of 3. Cannot continue')

        pyfastaq.sequences.genetic_code = genetic_code
        protein_seq = sequence.translate()

        if sequence.seq[0:3].upper() not in pyfastaq.genetic_codes.starts[genetic_code]:
            raise Error('Sequence "' + sequence.id + '" does not start with a start codon. Cannot continue')
        elif protein_seq[-1] != '*':
            raise Error('Sequence "' + sequence.id + '" does not end with a stop codon. Cannot continue')
        elif '*' in protein_seq[:-1]:
            raise Error('Sequence "' + sequence.id + '" has an internal stop codon. Cannot continue')

        return True


    @classmethod
    def _check_sequences(cls, padded_sequences, unpadded_sequences, seqs_are_coding, genetic_code=11):
        AlnToMetadata._check_seq_lengths_same(padded_sequences)

        if seqs_are_coding:
            for sequence in unpadded_sequences.values():
                AlnToMetadata._check_insertion_coords(sequence)
                AlnToMetadata._check_coding_seq(sequence, genetic_code=genetic_code)

        return True


    @classmethod
    def _check_variants_match_sequences(cls, unpadded_sequences, variants, seqs_are_coding, genetic_code=11):
        pyfastaq.sequences.genetic_code = genetic_code
        for seqname, variants in variants.items():
            if seqname not in unpadded_sequences:
                raise Error('Sequence name "' + seqname + '" given in variants file, but sequence not found')
            for variant, description in variants:
                if not variant.sanity_check_against_seq(unpadded_sequences[seqname], translate_seq=seqs_are_coding):
                    raise Error('Variant "' + str(variant) + '" for sequence "' + seqname + '" does not match sequence. cannot continue')
        return True

    def run(self, outprefix):
        pass
