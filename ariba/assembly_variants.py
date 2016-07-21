import operator
import pyfastaq
import pymummer
from ariba import sequence_variant
from pyfastaq import intervals


class Error (Exception): pass

class AssemblyVariants:
    def __init__(self,
      refdata,
      nucmer_snp_file,
    ):
        self.refdata = refdata
        self.nucmer_snp_file = nucmer_snp_file


    @classmethod
    def _get_codon_start(cls, gene_start, position):
        assert position >= gene_start
        while  (position - gene_start) % 3 != 0:
            position -= 1
        return position


    @classmethod
    def _get_mummer_variants(cls, snp_file):
        variants = pymummer.snp_file.get_all_variants(snp_file)
        mummer_variants = {}

        if len(variants) == 0:
            return {}

        variants.sort(key=operator.attrgetter('qry_name'))
        variants.sort(key=operator.attrgetter('ref_start'))

        for v in variants:
            if v.qry_name not in mummer_variants:
                mummer_variants[v.qry_name] = []
            mummer_variants[v.qry_name].append(v)

        for contig in mummer_variants:
            l = mummer_variants[contig]
            if len(l) > 1:
                new_l = [[l[0]]]
                previous_codon_start = AssemblyVariants._get_codon_start(0, l[0].ref_start)
                for variant in l[1:]:
                    codon_start = AssemblyVariants._get_codon_start(0, variant.ref_start)
                    if codon_start == previous_codon_start:
                        new_l[-1].append(variant)
                    else:
                        new_l.append([variant])
                        previous_codon_start = codon_start
                mummer_variants[contig] = new_l
            else:
                mummer_variants[contig] = [l]

        return mummer_variants


    @classmethod
    def _get_variant_effect(cls, variants, ref_sequence):
        '''variants = list of variants in the same codon.
           returns type of variant (cannot handle more than one indel in the same codon).'''
        assert len(variants) != 0

        var_types = [x.var_type for x in variants]
        if len(set(var_types)) != 1:
            return 'MULTIPLE', '.', '.'

        var_type = var_types[0]

        assert set([x.ref_name for x in variants]) == set([ref_sequence.id])
        codon_starts = [AssemblyVariants._get_codon_start(0, x.ref_start) for x in variants]
        assert len(set(codon_starts)) == 1
        codon_start = codon_starts[0]
        aa_start = codon_start // 3
        ref_codon = pyfastaq.sequences.Fasta('codon', ref_sequence[codon_start:codon_start+3])
        ref_aa = ref_codon.translate()

        if var_type == pymummer.variant.SNP:
            new_codon = list(ref_codon.seq)
            for v in variants:
                new_codon[v.ref_start - codon_start] = v.qry_base
            new_codon = pyfastaq.sequences.Fasta('new', ''.join(new_codon))
            qry_aa = new_codon.translate()

            if ref_aa.seq == qry_aa.seq:
                return ('SYN', '.', aa_start)
            elif qry_aa.seq == '*':
                return ('TRUNC', ref_aa.seq + str(aa_start + 1) + 'trunc', aa_start)
            else:
                return ('NONSYN', ref_aa.seq + str(aa_start + 1) + qry_aa.seq, aa_start)
        elif var_type in [pymummer.variant.INS, pymummer.variant.DEL]:
            if len(variants) > 1:
                return 'INDELS', '.', aa_start

            var = variants[0]

            if var_type == pymummer.variant.INS:
                new_seq = pyfastaq.sequences.Fasta('seq', var.qry_base)
            else:
                new_seq = pyfastaq.sequences.Fasta('seq', var.ref_base)

            if len(new_seq) % 3 != 0:
                return ('FSHIFT', ref_aa.seq + str(aa_start + 1) + 'fs', aa_start)

            new_seq_aa = new_seq.translate()
            if '*' in new_seq_aa.seq:
                return ('TRUNC', ref_aa.seq + str(aa_start + 1) + 'trunc', aa_start)
            elif var_type == pymummer.variant.INS:
                ref_codon_after_ins = pyfastaq.sequences.Fasta('codon', ref_sequence[codon_start+3:codon_start+6])
                aa_after_ins = ref_codon_after_ins.translate()
                return ('INS', ref_aa.seq + str(aa_start + 1) + '_' + aa_after_ins.seq + str(aa_start + 2) + 'ins' + new_seq_aa.seq , aa_start)
            else:
                if len(new_seq) == 3:
                    return ('DEL', ref_aa.seq + str(aa_start + 1) + 'del', aa_start)
                else:
                    assert len(new_seq) % 3 == 0
                    ref_codon_after_ins = pyfastaq.sequences.Fasta('codon', ref_sequence[codon_start+3:codon_start+6])
                    aa_after_ins = ref_codon_after_ins.translate()
                    return ('DEL', ref_aa.seq + str(aa_start + 1)+ '_' + aa_after_ins.seq + str(aa_start + 2) + 'del', aa_start)

        else:
            return ('UNKNOWN', '.', aa_start)


    @staticmethod
    def _filter_mummer_variants(mummer_variants, ref_sequence):
        if len(mummer_variants) == 0:
            return

        for contig in mummer_variants:
            variants = mummer_variants[contig]
            for i in range(len(variants)):
                t = AssemblyVariants._get_variant_effect(variants[i], ref_sequence)
                if t is not None and t[0] in ['TRUNC', 'FSHIFT']:
                    break
            mummer_variants[contig] = variants[:i+1]


    @staticmethod
    def _get_one_variant_for_one_contig_non_coding(refdata_var_dict, mummer_variant):
        var_tuple = None
        used_known_variants = set()

        # if the variant is at the same position as a known variant in the reference
        if refdata_var_dict is not None and mummer_variant.ref_start in refdata_var_dict['n']:
            if mummer_variant.var_type == pymummer.variant.SNP:
                variants_at_this_position = {x for x in refdata_var_dict['n'][mummer_variant.ref_start]}
                matching_variants = {x for x in variants_at_this_position if mummer_variant.qry_base == x.variant.variant_value}
                not_interesting_variants = {x for x in variants_at_this_position if mummer_variant.qry_base == x.variant.wild_value}
                variants_at_this_position = variants_at_this_position.difference(matching_variants)
            else:
                matching_variants = set()
                variants_at_this_position = refdata_var_dict['n'][mummer_variant.ref_start]
                not_interesting_variants = set()

            if len(not_interesting_variants) == 0:
                var_tuple = (
                    mummer_variant.ref_start,
                    'n',
                    mummer_variant.ref_base + str(mummer_variant.ref_start + 1) + mummer_variant.qry_base,
                    pymummer.variant.var_types[mummer_variant.var_type],
                    [mummer_variant],
                    matching_variants,
                    variants_at_this_position
                )

            used_known_variants.update(matching_variants, variants_at_this_position)
        else: # not at a known variant position in the reference
            var_tuple = (
                mummer_variant.ref_start,
                'n',
                mummer_variant.ref_base + str(mummer_variant.ref_start + 1) + mummer_variant.qry_base,
                pymummer.variant.var_types[mummer_variant.var_type],
                [mummer_variant],
                set(),
                set()
            )

        return var_tuple, used_known_variants


    @staticmethod
    def _get_one_variant_for_one_contig_coding(ref_sequence, refdata_var_dict, mummer_variants_list):
        aa_var_effect, aa_var_string, aa_var_position = AssemblyVariants._get_variant_effect(mummer_variants_list, ref_sequence)
        var_tuple = None
        used_known_variants = set()

        # if this variant is at the same position as a known variant in the reference
        if refdata_var_dict is not None and aa_var_position in refdata_var_dict['p']:
            if aa_var_effect == 'NONSYN':
                aa_variant = sequence_variant.Variant('p', aa_var_string, '.')
                variants_at_this_position = {x for x in refdata_var_dict['p'][aa_variant.position]}
                matching_variants = {x for x in variants_at_this_position if aa_variant.variant_value == x.variant.variant_value}
                not_interesting_variants = {x for x in variants_at_this_position if aa_variant.variant_value == x.variant.wild_value}
                variants_at_this_position = variants_at_this_position.difference(matching_variants)
            else:
                matching_variants = set()
                variants_at_this_position = refdata_var_dict['p'][aa_var_position]
                not_interesting_variants = set()

            if len(not_interesting_variants) == 0:
                var_tuple = (
                    aa_var_position,
                    'p',
                    aa_var_string,
                    aa_var_effect,
                    mummer_variants_list,
                    matching_variants,
                    variants_at_this_position
                )

            used_known_variants.update(matching_variants, variants_at_this_position)
        else: # this variant is not at a known position in the reference
            var_tuple = (
                aa_var_position,
                'p',
                aa_var_string,
                aa_var_effect,
                mummer_variants_list,
                set(),
                set()
            )

        return var_tuple, used_known_variants


    @staticmethod
    def _get_remaining_known_ref_variants(known_ref_variants, used_ref_variants, nucmer_coords):
        '''Finds variants where ref has the variant and so does the contig. Which means
           that there was no mummer call to flag it up so need to look through the known
           ref variants. Also need to check that the variant is in a nucmer match to an
           assembly contig.'''
        variants = []

        for ref_variant_pos, ref_variants_set in sorted(known_ref_variants.items()):
            for known_ref_variant in ref_variants_set:
                if known_ref_variant not in used_ref_variants:
                    variant_pos_matches_contig = False
                    pos = known_ref_variant.variant.position

                    if known_ref_variant.seq_type == 'n':
                        ref_interval = intervals.Interval(pos, pos)
                    elif known_ref_variant.seq_type == 'p':
                        ref_interval = intervals.Interval(3 * pos, 3 * pos + 2)
                    else:
                        raise Error('Unexpected variant type "' + known_ref_variant.variant_type + '" in _get_remaining_known_ref_variants. Cannot continue')

                    for interval in nucmer_coords:
                        if ref_interval.intersects(interval):
                            variant_pos_matches_contig = True
                            break

                    if variant_pos_matches_contig:
                        variants.append((None, known_ref_variant.seq_type, None, None, None, {known_ref_variant}, set()))

        return variants


    def get_variants(self, ref_sequence_name, nucmer_coords):
        '''Nucmr coords = dict. Key=contig name. Value = list of intervals of ref coords that match the contig.
           Made by assembly_compare.AssemblyCompare.nucmer_hits_to_ref_coords
           Returns dictionary. Key=contig name. Value = list of variants. Each variant
           is a tuple: (
               0 = position,
               1 = type in ['n', 'p']
               2 = Variant string, eg 'D2E',
               3 = variant effect (as returned by _get_variant_effect)
               4 = list of pymummer.variant.Variant that made up this variant (could be more than one because of
                   variants in the same codon)
               5 = set {matching known variants from metadata (=sequence_metadata.SequenceMetadata)}
               6 = set {known ref metadata (=sequence_metadata.SequenceMetadata)  at same position as SNP}, excluding those from 4
           )
        '''
        mummer_variants = self._get_mummer_variants(self.nucmer_snp_file)
        variants = {}
        seq_type, is_variant_only = self.refdata.sequence_type(ref_sequence_name)
        ref_sequence = self.refdata.sequence(ref_sequence_name)

        if ref_sequence_name in self.refdata.metadata:
            refdata_var_dict = self.refdata.metadata[ref_sequence_name]
        else:
            refdata_var_dict = None

        known_non_wild_variants_in_ref = self.refdata.all_non_wild_type_variants(ref_sequence_name)

        for contig in nucmer_coords:
            used_known_variants = set()
            variants[contig] = []

            if contig in mummer_variants:
                for mummer_variant_list in mummer_variants[contig]:
                    if seq_type == 'p':
                        new_variant, used_variants = self._get_one_variant_for_one_contig_coding(ref_sequence, refdata_var_dict, mummer_variant_list)
                    else:
                        for mummer_variant in mummer_variant_list:
                            new_variant, used_variants = self._get_one_variant_for_one_contig_non_coding(refdata_var_dict, mummer_variant)

                    if new_variant is not None:
                            variants[contig].append(new_variant)
                    used_known_variants.update(used_variants)

            # for this contig, need to know all the ref sequence and coords it maps to.
            # Then report just the unused known variants, as the contig also has these variants
            if seq_type == 'p':
                new_variants = self._get_remaining_known_ref_variants(known_non_wild_variants_in_ref['p'], used_known_variants, nucmer_coords[contig])

            else:
                new_variants = self._get_remaining_known_ref_variants(known_non_wild_variants_in_ref['n'], used_known_variants, nucmer_coords[contig])

            if is_variant_only:
                new_variants = [x for x in new_variants if len(x[5]) > 0]

            variants[contig].extend(new_variants)
            if len(variants[contig]) == 0:
                del variants[contig]

        return variants

