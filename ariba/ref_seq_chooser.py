import tempfile
import copy
import os
import pymummer
import pyfastaq
from ariba import common

class Error (Exception): pass

class RefSeqChooser:
    def __init__(self,
        cluster_fasta,
        all_refs_fasta,
        assembly_fasta_in,
        assembly_fasta_out,
        log_fh,
        contig_name_prefix='',
        nucmer_min_id=90,
        nucmer_min_len=20,
        nucmer_breaklen=200,
    ):
        self.cluster_fasta = os.path.abspath(cluster_fasta)
        self.all_refs_fasta = os.path.abspath(all_refs_fasta)
        self.assembly_fasta_in = os.path.abspath(assembly_fasta_in)
        self.assembly_fasta_out = os.path.abspath(assembly_fasta_out)
        self.log_fh = log_fh
        self.nucmer_min_id = nucmer_min_id
        self.nucmer_min_len = nucmer_min_len
        self.nucmer_breaklen = nucmer_breaklen
        self.closest_ref_within_cluster = None
        self.closest_ref_from_all_refs = None
        self.closest_ref_is_in_cluster = False


    @classmethod
    def _make_matching_contig_pieces_fasta(cls, contig_fasta, coords_list, outfile):
        matches = {}
        for hit in coords_list:
            if hit.qry_name not in matches:
                matches[hit.qry_name] = []
            matches[hit.qry_name].append(hit.qry_coords())

        if len(matches) == 0:
            return False

        contig_seqs = {}
        pyfastaq.tasks.file_to_dict(contig_fasta, contig_seqs)

        f_out = pyfastaq.utils.open_file_write(outfile)

        for contig_name in sorted(matches):
            pyfastaq.intervals.merge_overlapping_in_list(matches[contig_name])
            for interval in matches[contig_name]:
                subseq = contig_seqs[contig_name].subseq(interval.start, interval.end+1)
                subseq.id = contig_name + ':' + str(interval.start+1) + '-' + str(interval.end+1)
                print(subseq, file=f_out)

        pyfastaq.utils.close(f_out)
        return True


    @classmethod
    def _sequence_is_in_fasta_file(cls, seq_name, fasta_file):
        for seq in pyfastaq.sequences.file_reader(fasta_file):
            if seq.id == seq_name:
                return True
        return False


    @classmethod
    def _load_nucmer_coords_file(cls, coords_file, log_fh=None):
        matches = {}
        file_reader = pymummer.coords_file.reader(coords_file)
        for hit in file_reader:
            if log_fh is not None:
                print('[choose ref nucmer]', hit, sep='\t', file=log_fh)
            if hit.ref_name not in matches:
                matches[hit.ref_name] = []
            matches[hit.ref_name].append(copy.copy(hit))
        return matches


    @classmethod
    def _l_and_c_from_contig_name(cls, name):
        try:
            prefix, l_string, c_string, ctg, counter = name.split('.')
            l = int(l_string[1:])
            c = int(c_string[1:])
        except:
            raise Error('Error parsing contig name "' + name + '"')

        return l, c


    @classmethod
    def _best_of_two_hits(cls, hit1, hit2, use_qry_length=False, check_flanking=False):
        if use_qry_length:
            qry_length_percent1 = hit1.hit_length_qry / hit1.qry_length
            qry_length_percent2 = hit2.hit_length_qry / hit2.qry_length
            if qry_length_percent1 != qry_length_percent2:
                return hit1 if qry_length_percent1 > qry_length_percent2 else hit2

        ref_length_percent1 = hit1.hit_length_ref / hit1.ref_length
        ref_length_percent2 = hit2.hit_length_ref / hit2.ref_length
        if ref_length_percent1 != ref_length_percent2:
            return hit1 if ref_length_percent1 > ref_length_percent2 else hit2
        elif hit1.percent_identity != hit2.percent_identity:
            return hit1 if hit1.percent_identity > hit2.percent_identity else hit2
        else:
            if check_flanking:
                flank1 = min(min(hit1.qry_start, hit1.qry_end), hit1.qry_length - 1 - max(hit1.qry_start, hit1.qry_end))
                flank2 = min(min(hit2.qry_start, hit2.qry_end), hit2.qry_length - 1 - max(hit2.qry_start, hit2.qry_end))
                if flank1 != flank2:
                    return hit1 if flank1 > flank2 else hit2
            l1, c1 = RefSeqChooser._l_and_c_from_contig_name(hit1.qry_name)
            l2, c2 = RefSeqChooser._l_and_c_from_contig_name(hit2.qry_name)
            if l1 != l2:
                return hit1 if l1 > l2 else hit2
            else:
                return hit1 if c1 > c2 else hit2


    @classmethod
    def _choose_best_nucmer_match(cls, matches, use_qry_length=False, check_flanking=False):
        best_match = None
        for ref_name in matches:
            for hit in matches[ref_name]:
                if best_match is None:
                    best_match = hit
                else:
                    best_match = RefSeqChooser._best_of_two_hits(best_match, hit, use_qry_length=use_qry_length, check_flanking=check_flanking)

        return best_match


    @classmethod
    def _closest_nucmer_match_between_fastas(cls, ref_fasta, qry_fasta, log_fh, min_id, min_length, breaklen, use_qry_length, check_flanking):
        tmpdir = tempfile.mkdtemp(prefix='tmp.closest_nucmer_match.', dir=os.getcwd())
        coords_file = os.path.join(tmpdir, 'nucmer_vs_cluster_refs.coords')
        pymummer.nucmer.Runner(
            ref_fasta,
            qry_fasta,
            coords_file,
            min_id=min_id,
            min_length=min_length,
            breaklen=breaklen,
            maxmatch=True,
        ).run()
        nucmer_matches = RefSeqChooser._load_nucmer_coords_file(coords_file, log_fh=log_fh)
        common.rmtree(tmpdir)

        if len(nucmer_matches) == 0:
            return None, {}
        else:
            best_hit = RefSeqChooser._choose_best_nucmer_match(nucmer_matches, use_qry_length=use_qry_length, check_flanking=check_flanking)
            return best_hit, nucmer_matches


    def run(self):
        print('Looking for closest match from sequences within cluster', file=self.log_fh)
        best_hit_from_cluster, nucmer_matches = RefSeqChooser._closest_nucmer_match_between_fastas(self.cluster_fasta, self.assembly_fasta_in, self.log_fh, self.nucmer_min_id, self.nucmer_min_len, self.nucmer_breaklen, False, True)
        if best_hit_from_cluster is None:
            return

        self.closest_ref_within_cluster = best_hit_from_cluster.ref_name
        contigs_prefix = best_hit_from_cluster.qry_name.rsplit('.', maxsplit=2)[0]
        print('Closest cluster ref sequence is', self.closest_ref_within_cluster, 'to assembly', contigs_prefix, file=self.log_fh)
        tmpdir = tempfile.mkdtemp(prefix='tmp.choose_ref.', dir=os.getcwd())
        pieces_fasta_file = os.path.join(tmpdir, 'nucmer_vs_cluster_refs.pieces.fa')
        pieces_coords = [hit for hit in nucmer_matches[best_hit_from_cluster.ref_name] if hit.qry_name.startswith(contigs_prefix + '.')]
        RefSeqChooser._make_matching_contig_pieces_fasta(self.assembly_fasta_in, pieces_coords, pieces_fasta_file)

        print('Checking for a better match to a ref sequence outside the cluster', file=self.log_fh)
        best_hit_from_all_seqs, not_needed = RefSeqChooser._closest_nucmer_match_between_fastas(self.all_refs_fasta, pieces_fasta_file, self.log_fh, self.nucmer_min_id, self.nucmer_min_len, self.nucmer_breaklen, True, False)
        common.rmtree(tmpdir)
        self.closest_ref_from_all_refs = best_hit_from_all_seqs.ref_name
        if self.closest_ref_from_all_refs is None:
            return

        self.closest_ref_is_in_cluster = RefSeqChooser._sequence_is_in_fasta_file(self.closest_ref_from_all_refs, self.cluster_fasta)

        if self.closest_ref_is_in_cluster:
            print('Closest reference is', self.closest_ref_from_all_refs, 'and belongs to cluster', file=self.log_fh)
            l_and_c_bits = contigs_prefix.rsplit('.', maxsplit=2)[1:]
            name_regex = r'''\.''' + l_and_c_bits[0] + r'''\.''' + l_and_c_bits[1] + r'''\.ctg\.[0-9]+$'''
            print('Writing best assembly (contigs matching ', name_regex + ') to fasta file', self.assembly_fasta_out, file=self.log_fh)
            pyfastaq.tasks.filter(self.assembly_fasta_in, self.assembly_fasta_out, regex=name_regex)
        else:
            print('Closest reference is', self.closest_ref_from_all_refs, 'and does not belong to cluster', file=self.log_fh)

