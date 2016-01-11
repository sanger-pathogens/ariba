import operator
import pyfastaq
import pymummer

class Error (Exception): pass

class AssemblyVariants:
    def __init__(self,
      refdata,
      ref_sequence,
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
        if len(variants) == 0:
            return None

        var_types = [x.var_type for x in variants]
        if len(set(var_types)) != 1:
            return None

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
                return ('SYN', '.')
            elif qry_aa.seq == '*':
                return ('TRUNC', ref_aa.seq + str(aa_start + 1) + 'trunc')
            else:
                return ('NONSYN', ref_aa.seq + str(aa_start + 1) + qry_aa.seq)
        elif var_type in [pymummer.variant.INS, pymummer.variant.DEL]:
            if len(variants) > 1:
                print('More than one indel in same codon not yet implemented!', ref_sequence.id, file=sys.stderr)
                return None

            var = variants[0]

            if var_type == pymummer.variant.INS:
                new_seq = pyfastaq.sequences.Fasta('seq', var.qry_base)
            else:
                new_seq = pyfastaq.sequences.Fasta('seq', var.ref_base)

            if len(new_seq) % 3 != 0:
                return ('FSHIFT', ref_aa.seq + str(aa_start + 1) + 'fs')

            new_seq_aa = new_seq.translate()
            if '*' in new_seq_aa.seq:
                return ('TRUNC', ref_aa.seq + str(aa_start + 1) + 'trunc')
            elif var_type == pymummer.variant.INS:
                ref_codon_after_ins = pyfastaq.sequences.Fasta('codon', ref_sequence[codon_start+3:codon_start+6])
                aa_after_ins = ref_codon_after_ins.translate()
                return ('INS', ref_aa.seq + str(aa_start + 1) + '_' + aa_after_ins.seq + str(aa_start + 2) + 'ins' + new_seq_aa.seq )
            else:
                if len(new_seq) == 3:
                    return ('DEL', ref_aa.seq + str(aa_start + 1) + 'del')
                else:
                    assert len(new_seq) % 3 == 0
                    new_aa = new_seq.translate()
                    ref_codon_after_ins = pyfastaq.sequences.Fasta('codon', ref_sequence[codon_start+3:codon_start+6])
                    aa_after_ins = ref_codon_after_ins.translate()
                    return ('DEL', ref_aa.seq + str(aa_start + 1)+ '_' + aa_after_ins.seq + str(aa_start + 2) + 'del')

        else:
            return ('UNKNOWN', '.')


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


    def gather_variants(self):
        pass
