import re
import sys
import pyfastaq
from ariba import sequence_variant

class Error (Exception): pass

class AlnToMetadata:
    def __init__(self,
      aln_file,
      vars_file,
      refs_are_coding,
      refs_are_variant_only,
      genetic_code=11,
    ):
        self.padded_seqs = AlnToMetadata._load_aln_file(aln_file)
        self.refs_are_coding = refs_are_coding
        self.refs_are_variant_only = refs_are_variant_only
        self.variants = AlnToMetadata._load_vars_file(vars_file, self.refs_are_coding)
        self.genetic_code = genetic_code


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
    def _make_unpadded_insertion_coords(cls, unpadded_sequences):
        return {x.id: AlnToMetadata._insertion_coords(x) for x in unpadded_sequences.values()}


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

        original_code = pyfastaq.sequences.genetic_code
        pyfastaq.sequences.genetic_code = genetic_code
        protein_seq = sequence.translate()
        start_ok = sequence.seq[0:3].upper() in pyfastaq.genetic_codes.starts[genetic_code]
        pyfastaq.sequences.genetic_code = original_code

        if not start_ok:
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
        original_code = pyfastaq.sequences.genetic_code
        pyfastaq.sequences.genetic_code = genetic_code
        for seqname, variant_list in variants.items():
            if seqname not in unpadded_sequences:
                pyfastaq.sequences.genetic_code = original_code
                raise Error('Sequence name "' + seqname + '" given in variants file, but sequence not found')
            for variant, description in variant_list:
                if not variant.sanity_check_against_seq(unpadded_sequences[seqname], translate_seq=seqs_are_coding):
                    pyfastaq.sequences.genetic_code = original_code
                    raise Error('Variant "' + str(variant) + '" for sequence "' + seqname + '" does not match sequence. cannot continue')

        pyfastaq.sequences.genetic_code = original_code
        return True


    @classmethod
    def _variant_ids_are_unique(cls, variants):
        seen_variants = set()
        for variants_list in variants.values():
            for variant, description in variants_list:
                if variant.identifier in seen_variants:
                    raise Error('Variant identifier "' + variant.identifier + '" found more than once. Cannot continue')
                else:
                    seen_variants.add(variant.identifier)

        return True


    @classmethod
    def _unpadded_to_padded_nt_position(cls, position, insertions):
        if len(insertions) == 0:
            return position

        i = 0
        while i < len(insertions) and insertions[i].start <= position:
            position += len(insertions[i])
            i += 1

        return position


    @classmethod
    def _padded_to_unpadded_nt_position(cls, position, insertions):
        if len(insertions) == 0:
            return position

        i = 0
        total_gap_length = 0
        while i < len(insertions) and insertions[i].end < position:
            total_gap_length += len(insertions[i])
            i += 1

        if i < len(insertions) and insertions[i].distance_to_point(position) == 0:
            return None
        else:
            return position - total_gap_length


    @classmethod
    def _variants_to_tsv_lines(cls, variants, unpadded_sequences, padded_sequences, insertions, seqs_are_coding, seqs_are_var_only):
        if seqs_are_coding:
            unpadded_aa_sequences = {x: unpadded_sequences[x].translate() for x in unpadded_sequences}
            is_gene = '1'
        else:
            is_gene = '0'

        is_var_only = '1' if seqs_are_var_only else '0'

        lines = []
        for refname in sorted(variants):
            for variant, description in variants[refname]:
                if seqs_are_coding:
                    ref_unpadded_nt_position = 3 * variant.position
                else:
                    ref_unpadded_nt_position = variant.position

                padded_nt_position = AlnToMetadata._unpadded_to_padded_nt_position(ref_unpadded_nt_position, insertions[refname])
                lines.append('\t'.join([refname, is_gene, is_var_only, str(variant), variant.identifier, description]))

                for seqname, seq in sorted(padded_sequences.items()):
                    if seqname == refname:
                        continue

                    if seq[padded_nt_position] == '-':
                        print('Warning: position has a gap in sequence ', seqname, 'corresponding to variant', variant, '(' + variant.identifier + ') in sequence ', refname, '... Ignoring for ' + seqname, file=sys.stderr)
                        continue

                    unpadded_nt_position = AlnToMetadata._padded_to_unpadded_nt_position(padded_nt_position, insertions[seqname])
                    assert unpadded_nt_position is not None

                    if seqs_are_coding:
                        assert unpadded_nt_position % 3 == 0
                        unpadded_aa_position = unpadded_nt_position // 3

                        if unpadded_aa_sequences[seqname][unpadded_aa_position] in {variant.wild_value, variant.variant_value}:
                            variant_string = variant.wild_value
                        else:
                            variant_string = unpadded_aa_sequences[seqname][unpadded_aa_position]
                        variant_string += str(unpadded_aa_position + 1) + variant.variant_value
                    else:
                        if unpadded_sequences[seqname][unpadded_nt_position] in {variant.wild_value, variant.variant_value}:
                            variant_string = variant.wild_value
                        else:
                            variant_string = unpadded_sequences[seqname][unpadded_nt_position]
                        variant_string += str(unpadded_nt_position + 1) + variant.variant_value

                    lines.append('\t'.join([seqname, is_gene, is_var_only, variant_string, variant.identifier, description]))

        return lines


    @classmethod
    def _make_cluster_file(cls, sequences, filename):
        names = sorted(sequences.keys())
        with open(filename, 'w') as f:
            print(*names, sep='\t', file=f)


    def run(self, outprefix):
        original_code = pyfastaq.sequences.genetic_code
        pyfastaq.sequences.genetic_code = self.genetic_code
        unpadded_seqs = AlnToMetadata._make_unpadded_seqs(self.padded_seqs)
        insertions = AlnToMetadata._make_unpadded_insertion_coords(self.padded_seqs)
        AlnToMetadata._check_sequences(self.padded_seqs, unpadded_seqs, self.refs_are_coding, genetic_code=self.genetic_code)
        AlnToMetadata._variant_ids_are_unique(self.variants)
        AlnToMetadata._check_variants_match_sequences(unpadded_seqs, self.variants, self.refs_are_coding, genetic_code=self.genetic_code)

        tsv_lines = AlnToMetadata._variants_to_tsv_lines(self.variants, unpadded_seqs, self.padded_seqs, insertions, self.refs_are_coding, self.refs_are_variant_only)
        with open(outprefix + '.tsv', 'w') as f:
            print(*tsv_lines, sep='\n', file=f)

        with open(outprefix + '.fa', 'w') as f:
            for seqname in sorted(unpadded_seqs):
                print(unpadded_seqs[seqname], sep='\n', file=f)

        AlnToMetadata._make_cluster_file(unpadded_seqs, outprefix + '.cluster')
        pyfastaq.sequences.genetic_code = original_code
