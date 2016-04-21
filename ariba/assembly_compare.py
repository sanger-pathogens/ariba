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

        self.nucmer_coords_file = self.outprefix + '.nucmer.coords'
        self.nucmer_snps_file = self.nucmer_coords_file + '.snps'
        self.assembled_ref_seqs_file = self.outprefix + '.assembled_refs.fasta'


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
    def _write_assembled_reference_sequences(nucmer_hits, ref_sequence, assembly, outfile):
        '''nucmer_hits =  hits made by self._parse_nucmer_coords_file.
           ref_gene = reference sequence (pyfastaq.sequences.Fasta object)
           assembly = dictionary of contig name -> contig.
           Writes each piece of assembly that corresponds to the reference sequence
           to a fasta file.'''
        f = pyfastaq.utils.open_file_write(outfile)

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

                print(fa, file=f)

        pyfastaq.utils.close(f)


    @staticmethod
    def _whole_gene_covered_by_nucmer_hits(nucmer_hits, ref_seq, threshold):
        '''Returns true iff the reference sequence is covered by nucmer hits.
           nucmer_hits = hits made by self._parse_nucmer_coords_file.
           Counts as covered if (total ref bases covered) / len(ref_seq) >= threshold'''
        coords = AssemblyCompare.nucmer_hits_to_ref_coords(nucmer_hits)
        covered = []
        for coords_list in coords.values():
            covered.extend(coords_list)
        pyfastaq.intervals.merge_overlapping_in_list(covered)
        return pyfastaq.intervals.length_sum_from_list(covered) / len(ref_seq) >= threshold


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


    @staticmethod
    def _ref_covered_by_complete_contig_with_orf(nucmer_hits, contigs):
        '''Returns true iff there is a contig that covers the entire reference,
           and that contig has a complete open reading frame.
           nucmer_hits = hits made by self._parse_nucmer_coords_file.'''
        for l in nucmer_hits.values():
            for hit in l:
                if hit.hit_length_ref == hit.ref_length:
                    start = min(hit.qry_start, hit.qry_end)
                    end = max(hit.qry_start, hit.qry_end)
                    assembled_gene = pyfastaq.sequences.Fasta('x', contigs[hit.qry_name][start:end+1])
                    if (hit.ref_start < hit.ref_end) != (hit.qry_start < hit.qry_end):
                        assembled_gene.revcomp()
                    orfs = assembled_gene.orfs()
                    if len(orfs) == 0:
                        continue

                    max_orf = orfs[0]
                    for o in orfs:
                        if len(o) > len(max_orf):
                            max_orf = o

                    if len(max_orf) == len(assembled_gene):
                        return True
        return False


    @staticmethod
    def _ref_covered_by_at_least_one_full_length_contig(nucmer_hits):
        '''Returns true iff there exists a contig that completely
           covers the reference sequence
           nucmer_hits = hits made by self._parse_nucmer_coords_file.'''
        for l in nucmer_hits.values():
            for hit in l:
                if len(hit.ref_coords()) == hit.ref_length:
                    return True
        return False


    def update_flag(self, flag):
        if self._whole_gene_covered_by_nucmer_hits(self.nucmer_hits, self.ref_sequence, self.assembled_threshold):
            flag.add('assembled')

        if self._ref_covered_by_at_least_one_full_length_contig(self.nucmer_hits):
            flag.add('assembled_into_one_contig')

        if self._ref_has_region_assembled_twice(self.nucmer_hits, self.ref_sequence, self.unique_threshold):
            flag.add('region_assembled_twice')

        ref_seq_type = self.refdata.sequence_type(self.ref_sequence.id)
        if ref_seq_type != 'non_coding' and self._ref_covered_by_complete_contig_with_orf(self.nucmer_hits, self.assembly_sequences):
            flag.add('complete_orf')

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
        self._write_assembled_reference_sequences(self.nucmer_hits, self.ref_sequence, self.assembly_sequences, self.assembled_ref_seqs_file)
