import pysam
import pyfastaq
import os
from ariba import scaffold_graph

class Error (Exception): pass

class Parser:
    def __init__(self, bam, ref_seqs):
        '''Construct a Parser.
        bam: name of BAM file
        ref_seqs: dictionary of sequence name => Fasta object'''
        self.bam = os.path.abspath(bam)
        self.soft_clipped = {}
        self.ref_lengths = {seq: len(ref_seqs[seq]) for seq in ref_seqs}
        self.scaff_graph = scaffold_graph.Graph(self.ref_lengths)
        self.unmapped_mates = {}
        self.sam_reader = pysam.Samfile(self.bam, "rb")


    def _sam_to_soft_clipped(self, sam):
        '''Returns tuple of whether or not the left and right end of the mapped read in the sam record is soft-clipped'''
        if sam.is_unmapped:
            raise Error('Cannot get soft clip info from an unmapped read')

        if sam.cigar is None or len(sam.cigar) == 0:
            return False, False
        return (sam.cigar[0][0] == 4, sam.cigar[-1][0] == 4)


    def _update_soft_clipped_from_sam(self, sam):
        ref_name = self.sam_reader.getrname(sam.tid)
        left_clip, right_clip = self._sam_to_soft_clipped(sam)

        if left_clip or right_clip:
            if ref_name not in self.soft_clipped:
                self.soft_clipped[ref_name] = {}

            if left_clip:
                p = sam.pos
                if p not in self.soft_clipped[ref_name]:
                    self.soft_clipped[ref_name][p] = [0, 0]
                self.soft_clipped[ref_name][p][0] += 1

            if right_clip:
                p = sam.aend
                if p not in self.soft_clipped[ref_name]:
                    self.soft_clipped[ref_name][p] = [0, 0]
                self.soft_clipped[ref_name][p][1] += 1


    def _update_unmapped_mates_from_sam(self, sam):
        if (not sam.is_unmapped) and sam.mate_is_unmapped:
            ref_name = self.sam_reader.getrname(sam.tid)
            if ref_name not in self.unmapped_mates:
                self.unmapped_mates[ref_name] = {}
            pos = sam.reference_start
            self.unmapped_mates[ref_name][pos] = self.unmapped_mates[ref_name].get(pos, 0) + 1


    def _write_soft_clipped_to_file(self, filename):
        f = pyfastaq.utils.open_file_write(filename)
        for seq in sorted(self.soft_clipped):
            for position in sorted(self.soft_clipped[seq]):
                print(seq, position + 1, self.soft_clipped[seq][position][0], self.soft_clipped[seq][position][1], sep='\t', file=f)
        pyfastaq.utils.close(f)


    def _write_unmapped_mates_to_file(self, filename):
        f = pyfastaq.utils.open_file_write(filename)
        for seq in sorted(self.unmapped_mates):
            for position in sorted(self.unmapped_mates[seq]):
                print(seq, position + 1, self.unmapped_mates[seq][position], sep='\t', file=f)
        pyfastaq.utils.close(f)


    def parse(self):
        for sam in self.sam_reader.fetch(until_eof=True):
            if sam.is_unmapped:
                continue

            self._update_soft_clipped_from_sam(sam)
            self.scaff_graph.update_from_sam(sam, self.sam_reader)
            self._update_unmapped_mates_from_sam(sam)


    def scaff_graph_is_consistent(self, min_coverage, max_insert):
        return self.scaff_graph.is_consistent(min_coverage, max_insert)


    def write_files(self, prefix):
        self._write_soft_clipped_to_file(prefix + '.soft_clipped')
        self._write_unmapped_mates_to_file(prefix + '.unmapped_mates')
        self.scaff_graph.write_all_links_to_file(prefix + '.scaff')

