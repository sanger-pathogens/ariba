import os
import copy
import pyfastaq
import pymummer

class Error (Exception): pass

class AssemblyCompare:
    def __init__(self,
      assembly_fa,
      assembly_sequences,
      ref_fa,
      ref_sequence,
      outprefix,
      refdata,
      nucmer_min_id=90,
      nucmer_min_len=20,
      nucmer_breaklen=200,
      assembled_threshold=0.95,
      unique_threshold=0.03,
      max_gene_nt_extend=30,
    ):
        self.assembly_fa = os.path.abspath(assembly_fa)
        self.assembly_sequences = assembly_sequences
        self.ref_fa = os.path.abspath(ref_fa)
        self.ref_sequence = ref_sequence
        self.outprefix = os.path.abspath(outprefix)
        self.refdata = refdata

        self.nucmer_min_id = nucmer_min_id
        self.nucmer_min_len = nucmer_min_len
        self.nucmer_breaklen = nucmer_breaklen
        self.assembled_threshold = assembled_threshold
        self.unique_threshold = unique_threshold
        self.max_gene_nt_extend = max_gene_nt_extend
        self.scaff_name_matching_ref = None
        self.gene_matching_ref = None
        self.gene_matching_ref_type = None
        self.gene_start_bases_added = None
        self.gene_end_bases_added = None

        self.nucmer_coords_file = self.outprefix + '.nucmer.coords'
        self.nucmer_snps_file = self.nucmer_coords_file + '.snps'


    def _run_nucmer(self):
        pymummer.nucmer.Runner(
            self.ref_fa,
            self.assembly_fa,
            self.nucmer_coords_file,
            min_id=self.nucmer_min_id,
            min_length=self.nucmer_min_len,
            breaklen=self.nucmer_breaklen,
            maxmatch=True,
            show_snps=True
        ).run()


    @staticmethod
    def _parse_nucmer_coords_file(coords_file, ref_name):
        '''Input is coords file made by self._run_nucmer. Reference should have one sequence only.
           ref_name is name fo the reference sequence, to sanity check the coords file.
           Returns dictionary. Key = assembly contig name. Value = list of nucmer hits to that contig'''
        file_reader = pymummer.coords_file.reader(coords_file)
        nucmer_hits = {}
        for hit in file_reader:
            assert hit.ref_name == ref_name
            contig = hit.qry_name
            if contig not in nucmer_hits:
                nucmer_hits[contig] = []
            nucmer_hits[contig].append(copy.copy(hit))

        return nucmer_hits


    @staticmethod
    def _nucmer_hits_to_percent_identity(nucmer_hits):
        '''Input is hits made by self._parse_nucmer_coords_file.
           Returns dictionary. key = contig name. Value = percent identity of hits to that contig'''
        percent_identities = {}

        for contig in nucmer_hits:
            product_sum = 0
            length_sum = 0
            for hit in nucmer_hits[contig]:
                product_sum += hit.hit_length_qry * hit.percent_identity
                length_sum += hit.hit_length_qry
            assert length_sum > 0
            percent_identities[contig] = round(product_sum / length_sum, 2)

        return percent_identities


    @staticmethod
    def _nucmer_hits_to_assembly_coords(nucmer_hits):
        '''Input is hits made by self._parse_nucmer_coords_file.
           Returns dictionary. key = contig name. Value = list of coords that match
           to the reference gene'''
        coords = {}
        for l in nucmer_hits.values():
            for hit in l:
                if hit.qry_name not in coords:
                    coords[hit.qry_name] = []
                coords[hit.qry_name].append(hit.qry_coords())

        for scaff in coords:
            pyfastaq.intervals.merge_overlapping_in_list(coords[scaff])

        return coords


    def assembly_match_coords(self):
        return self._nucmer_hits_to_assembly_coords(self.nucmer_hits)


    @classmethod
    def nucmer_hits_to_ref_coords(cls, nucmer_hits, contig=None):
        '''Input is hits made by self._parse_nucmer_coords_file.
           Returns dictionary. Key = contig name. Value = list of coords in the
           reference sequence for that contig.
           if contig=contig_name, then just gets the ref coords from that contig,
           instead of using all the contigs'''
        coords = []
        if contig is None:
            coords = {key: [] for key in nucmer_hits.keys()}
        else:
            coords = {contig: []}

        for key in coords:
            coords[key] = [hit.ref_coords() for hit in nucmer_hits[key]]
            pyfastaq.intervals.merge_overlapping_in_list(coords[key])

        return coords


    @staticmethod
    def ref_cov_per_contig(nucmer_hits):
        '''Input is hits made by self._parse_nucmer_coords_file.
           Returns dictionary. key = contig name. Value = number of bases that
           match to the reference sequence.'''
        coords = AssemblyCompare.nucmer_hits_to_ref_coords(nucmer_hits)
        return {x: pyfastaq.intervals.length_sum_from_list(coords[x]) for x in coords}


    @staticmethod
    def _get_assembled_reference_sequences(nucmer_hits, ref_sequence, assembly):
        '''nucmer_hits =  hits made by self._parse_nucmer_coords_file.
           ref_gene = reference sequence (pyfastaq.sequences.Fasta object)
           assembly = dictionary of contig name -> contig.
           Makes a set of Fasta objects of each piece of assembly that
           corresponds to the reference sequeunce.'''
        sequences = {}

        for contig in sorted(nucmer_hits):
            for hit in nucmer_hits[contig]:
                qry_coords = hit.qry_coords()
                fa = assembly[hit.qry_name].subseq(qry_coords.start, qry_coords.end + 1)
                if hit.on_same_strand():
                    strand = '+'
                else:
                    fa.revcomp()
                    strand = '-'
                ref_coords = hit.ref_coords()
                fa.id = '.'.join([
                    ref_sequence.id,
                    str(ref_coords.start + 1),
                    str(ref_coords.end + 1),
                    contig,
                    str(qry_coords.start + 1),
                    str(qry_coords.end + 1),
                    strand
                ])

                if hit.hit_length_ref == hit.ref_length:
                    fa.id += '.complete'

                sequences[fa.id] = fa

        return sequences


    @staticmethod
    def _whole_gene_covered_by_nucmer_hits(nucmer_hits, ref_seq, percent_threshold, max_nt_extend):
        '''Returns true iff the reference sequence is covered by nucmer hits.
           nucmer_hits = hits made by self._parse_nucmer_coords_file.
           Counts as covered if (total ref bases covered) / len(ref_seq) >= threshold'''
        coords = AssemblyCompare.nucmer_hits_to_ref_coords(nucmer_hits)
        covered = []
        for coords_list in coords.values():
            covered.extend(coords_list)
        pyfastaq.intervals.merge_overlapping_in_list(covered)
        return (2 * max_nt_extend + pyfastaq.intervals.length_sum_from_list(covered)) / len(ref_seq) >= percent_threshold


    @staticmethod
    def _ref_has_region_assembled_twice(nucmer_hits, ref_seq, threshold):
        '''Returns true iff there is a part of the reference that is assembled
           more than once (ie covered by >1 nucmer hit).
           Needs a minimum proportin of the ref to be assembled more than once,
           determined by threshold.
           nucmer_hits = hits made by self._parse_nucmer_coords_file.'''
        coords = AssemblyCompare.nucmer_hits_to_ref_coords(nucmer_hits)
        covered = []
        for coords_list in coords.values():
            covered.extend(coords_list)
        covered.sort()

        if len(covered) <= 1:
            return False

        coverage = {}
        for i in covered:
            for j in range(i.start, i.end + 1):
                coverage[j] = coverage.get(j, 0) + 1

        bases_depth_at_least_two = len([1 for x in coverage.values() if x > 1])
        return bases_depth_at_least_two / len(ref_seq) >= threshold


    @classmethod
    def _longest_nucmer_hit_in_ref(cls, nucmer_hits):
        max_length = None
        max_hit = None

        for l in nucmer_hits.values():
            for hit in l:
                if max_length is None or hit.hit_length_ref > max_length:
                    max_length = hit.hit_length_ref
                    max_hit = hit

        return max_hit


    @classmethod
    def _find_previous_start_codon(cls, sequence, start_coord, min_start_coord):
        for i in range(start_coord, min_start_coord - 1, -3):
            codon = pyfastaq.sequences.Fasta('x', sequence[i:i+3])
            aa = codon.translate()
            if aa.seq == '*':
                return None
            elif codon.seq in pyfastaq.genetic_codes.starts[pyfastaq.sequences.genetic_code]:
                return i

        return None


    @classmethod
    def _find_next_stop_codon(cls, sequence, end_coord, max_end_coord):
        final_i = min(len(sequence) - 3, max_end_coord - 2)
        for i in range(end_coord, final_i + 1, 3):
            codon = pyfastaq.sequences.Fasta('x', sequence[i:i+3])
            aa = codon.translate()
            if aa.seq == '*':
                return i

        return None


    @classmethod
    def _gene_from_nucmer_match(cls, nucmer_match, contig, max_end_nt_extend):
        if nucmer_match.on_same_strand():
            revcomp = False
        else:
            revcomp = True
            nucmer_match = copy.copy(nucmer_match)
            nucmer_match.reverse_query()
            contig = copy.copy(contig)
            contig.revcomp()

        ref_hit_start = min(nucmer_match.ref_start, nucmer_match.ref_end)
        contig_hit_start = min(nucmer_match.qry_start, nucmer_match.qry_end)
        min_allowed_start = max(0, contig_hit_start - max_end_nt_extend)
        if ref_hit_start % 3 != 0:
            contig_hit_start += 3 - (ref_hit_start % 3)
        contig_hit_end = max(nucmer_match.qry_start, nucmer_match.qry_end)
        max_allowed_end = min(len(contig) - 1, contig_hit_end + max_end_nt_extend)
        contig_hit_end -= (contig_hit_end - contig_hit_start + 1) % 3
        assert contig_hit_start < contig_hit_end
        gene_nt_name = nucmer_match.qry_name + '.' + str(contig_hit_start + 1) + '-' + str(contig_hit_end + 1)

        gene_nt = pyfastaq.sequences.Fasta(gene_nt_name, contig[contig_hit_start:contig_hit_end+1])
        assert len(gene_nt) % 3 == 0
        gene_aa = gene_nt.translate()
        if '*' in gene_aa[:-1]:
            return gene_nt, 'HAS_STOP', None, None

        extended_start_position = AssemblyCompare._find_previous_start_codon(contig, contig_hit_start, min_allowed_start)
        extended_end_position = AssemblyCompare._find_next_stop_codon(contig, contig_hit_end - 2, max_allowed_end)
        start = extended_start_position if extended_start_position is not None else contig_hit_start
        end = extended_end_position + 2 if extended_end_position is not None else contig_hit_end

        if revcomp:
            gene_nt_name = nucmer_match.qry_name + '.' + str(end + 1) + '-' + str(start + 1)
        else:
            gene_nt_name = nucmer_match.qry_name + '.' + str(start + 1) + '-' + str(end + 1)

        gene_nt = pyfastaq.sequences.Fasta(gene_nt_name, contig[start:end+1])
        start_nt_added = None if extended_start_position is None else min(nucmer_match.qry_start, nucmer_match.qry_end) - start
        end_nt_added = None if extended_end_position is None else end - max(nucmer_match.qry_start, nucmer_match.qry_end)

        if None in [extended_start_position, extended_end_position]:
            return gene_nt, 'START_OR_END_FAIL', start_nt_added, end_nt_added
        else:
            return gene_nt, 'GENE_FOUND', start_nt_added, end_nt_added


    @staticmethod
    def _get_gene_matching_ref(nucmer_hits, contigs, max_end_nt_extend):
        longest_match = AssemblyCompare._longest_nucmer_hit_in_ref(nucmer_hits)
        if longest_match is None:
            return None, 'NO_MATCH', None, None
        else:
            return (longest_match.qry_name,) + AssemblyCompare._gene_from_nucmer_match(longest_match, contigs[longest_match.qry_name], max_end_nt_extend)


    @staticmethod
    def _ref_covered_by_at_least_one_full_length_contig(nucmer_hits, percent_threshold, max_nt_extend):
        '''Returns true iff there exists a contig that completely
           covers the reference sequence
           nucmer_hits = hits made by self._parse_nucmer_coords_file.'''
        for l in nucmer_hits.values():
            for hit in l:
                if ( (2 * max_nt_extend) + len(hit.ref_coords()) ) / hit.ref_length >= percent_threshold:
                    return True
        return False


    def update_flag(self, flag):
        if self._whole_gene_covered_by_nucmer_hits(self.nucmer_hits, self.ref_sequence, self.assembled_threshold, self.max_gene_nt_extend):
            flag.add('assembled')

        if self.assembled_into_one_contig:
            flag.add('assembled_into_one_contig')

        if self._ref_has_region_assembled_twice(self.nucmer_hits, self.ref_sequence, self.unique_threshold):
            flag.add('region_assembled_twice')

        ref_seq_type, is_variant_only = self.refdata.sequence_type(self.ref_sequence.id)
        if ref_seq_type == 'p' and self.gene_matching_ref_type == 'GENE_FOUND':
            flag.add('complete_gene')

        if len(self.nucmer_hits) == 1:
            flag.add('unique_contig')

        return flag


    @staticmethod
    def nucmer_hit_containing_reference_position(nucmer_hits, ref_name, ref_position):
        '''Returns the first nucmer match found that contains the given
           reference location. nucmer_hits = hits made by self._parse_nucmer_coords_file.
           Returns None if no matching hit found'''
        for contig_name in nucmer_hits:
            for hit in nucmer_hits[contig_name]:
                if hit.ref_name == ref_name and hit.ref_coords().distance_to_point(ref_position) == 0:
                    return hit

        return None


    def run(self):
        self._run_nucmer()
        self.nucmer_hits = self._parse_nucmer_coords_file(self.nucmer_coords_file, self.ref_sequence.id)
        self.percent_identities = self._nucmer_hits_to_percent_identity(self.nucmer_hits)
        self.assembled_reference_sequences = self._get_assembled_reference_sequences(self.nucmer_hits, self.ref_sequence, self.assembly_sequences)
        ref_seq_type, is_variant_only = self.refdata.sequence_type(self.ref_sequence.id)
        if self._ref_covered_by_at_least_one_full_length_contig(self.nucmer_hits, self.assembled_threshold, self.max_gene_nt_extend):
            self.assembled_into_one_contig = True
            if ref_seq_type == 'p':
                self.scaff_name_matching_ref, self.gene_matching_ref, self.gene_matching_ref_type, self.gene_start_bases_added, self.gene_end_bases_added = self._get_gene_matching_ref(self.nucmer_hits, self.assembly_sequences, self.max_gene_nt_extend)
        else:
            self.assembled_into_one_contig = False
