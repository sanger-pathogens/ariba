import tempfile
import os
import pymummer
import pyfastaq
import shutil
from ariba import external_progs, mash

class RefSeqChooser:
    def __init__(self,
        cluster_fasta,
        all_refs_fasta,
        assembly_fasta,
        assembly_sequences,
        log_fh,
        extern_progs=None,
        nucmer_min_id=90,
        nucmer_min_len=20,
        nucmer_breaklen=200,
    ):
        self.cluster_fasta = os.path.abspath(cluster_fasta)
        self.all_refs_fasta = os.path.abspath(all_refs_fasta)
        self.assembly_fasta = os.path.abspath(assembly_fasta)
        self.assembly_sequences = assembly_sequences
        self.log_fh = log_fh
        self.nucmer_min_id = nucmer_min_id
        self.nucmer_min_len = nucmer_min_len
        self.nucmer_breaklen = nucmer_breaklen

        if extern_progs is None:
            self.extern_progs = external_progs.ExternalProgs()
        else:
            self.extern_progs = extern_progs

        self.closest_ref_within_cluster = None
        self.closest_ref_from_all_refs = None
        self.closest_ref_is_in_cluster = False


    @classmethod
    def _make_matching_contig_pieces_fasta(cls, contig_fasta, contig_seqs, coords_file, outfile, min_id=90, min_length=20, breaklen=200):
        matches = {}
        file_reader = pymummer.coords_file.reader(coords_file)
        for hit in file_reader:
            if hit.qry_name not in matches:
                matches[hit.qry_name] = []
            matches[hit.qry_name].append(hit.qry_coords())

        if len(matches) == 0:
            return False

        f_out = pyfastaq.utils.open_file_write(outfile)

        for contig_name in sorted(matches):
            pyfastaq.intervals.merge_overlapping_in_list(matches[contig_name])
            for interval in matches[contig_name]:
                subseq = contig_seqs[contig_name].subseq(interval.start, interval.end+1)
                subseq.id = contig_name + ':' + str(interval.start+1) + '-' + str(interval.end+1)
                print(subseq, file=f_out)

        pyfastaq.utils.close(f_out)


    @classmethod
    def _sequence_is_in_fasta_file(cls, seq_name, fasta_file):
        for seq in pyfastaq.sequences.file_reader(fasta_file):
            if seq.id == seq_name:
                return True
        return False


    def run(self):
        tmpdir = tempfile.mkdtemp(prefix='tmp.choose_ref.', dir=os.getcwd())

        mash_file_within_cluster = os.path.join(tmpdir, 'mash.closest_within_cluster.dist')
        masher = mash.Masher(self.cluster_fasta, self.assembly_fasta, self.log_fh, self.extern_progs)
        self.closest_ref_within_cluster = masher.run(mash_file_within_cluster)
        if self.closest_ref_within_cluster is None:
            shutil.rmtree(tmpdir)
            return

        coords_file = os.path.join(tmpdir, 'matching_pieces_to_ref.coords')

        pymummer.nucmer.Runner(
            self.cluster_fasta,
            self.assembly_fasta,
            coords_file,
            min_id=self.nucmer_min_id,
            min_length=self.nucmer_min_len,
            breaklen=self.nucmer_breaklen,
            maxmatch=True,
        ).run()

        pieces_matching_ref_fasta = os.path.join(tmpdir, 'pieces_matching_ref.fasta')
        RefSeqChooser._make_matching_contig_pieces_fasta(self.assembly_fasta, self.assembly_sequences, coords_file, pieces_matching_ref_fasta)

        mash_file_all_refs = os.path.join(tmpdir, 'mash.closest_ref_to_matching_pieces.dist')
        masher = mash.Masher(self.all_refs_fasta, pieces_matching_ref_fasta, self.log_fh, self.extern_progs)
        self.closest_ref_from_all_refs = masher.run(mash_file_all_refs)
        self.closest_ref_is_in_cluster = RefSeqChooser._sequence_is_in_fasta_file(self.closest_ref_from_all_refs, self.cluster_fasta)
        shutil.rmtree(tmpdir)
