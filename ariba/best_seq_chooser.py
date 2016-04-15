import shutil
import tempfile
import os
import pyfastaq
from ariba import mapping, faidx

class Error (Exception): pass

class BestSeqChooser:
    def __init__(self,
        reads1,
        reads2,
        references_fa,
        log_fh,
        samtools_exe='samtools',
        bowtie2_exe='bowtie2',
        bowtie2_preset='very-sensitive-local',
        threads=1,
    ):
        self.reads1 = reads1
        self.reads2 = reads2
        self.references_fa = references_fa
        self.log_fh = log_fh
        self.samtools_exe = samtools_exe
        self.bowtie2_exe = bowtie2_exe
        self.bowtie2_preset = bowtie2_preset
        self.threads = threads


    def _total_alignment_score(self, seq_name):
        tmpdir = tempfile.mkdtemp(prefix='tmp.get_total_aln_score.', dir=os.getcwd())
        tmp_bam = os.path.join(tmpdir, 'tmp.get_total_alignment_score.bam')
        tmp_fa = os.path.join(tmpdir, 'tmp.get_total_alignment_score.ref.fa')

        faidx.write_fa_subset(
            [seq_name],
            self.references_fa,
            tmp_fa,
            samtools_exe=self.samtools_exe,
            verbose=True,
            verbose_filehandle=self.log_fh
        )

        mapping.run_bowtie2(
            self.reads1,
            self.reads2,
            tmp_fa,
            tmp_bam[:-4],
            threads=self.threads,
            samtools=self.samtools_exe,
            bowtie2=self.bowtie2_exe,
            bowtie2_preset=self.bowtie2_preset,
            verbose=True,
            verbose_filehandle=self.log_fh
        )

        score = mapping.get_total_alignment_score(tmp_bam)
        shutil.rmtree(tmpdir)
        return score


    def _get_best_seq_by_alignment_score(self):
        total_sequences = pyfastaq.tasks.count_sequences(self.references_fa)
        if total_sequences == 1:
            seqs = {}
            pyfastaq.tasks.file_to_dict(self.references_fa, seqs)
            assert len(seqs) == 1
            seq_name = list(seqs.values())[0].id
            print('No need to choose sequence for this cluster because only has one sequence:', seq_name, file=self.log_fh)
            return seq_name

        print('\nChoosing best sequence from cluster of', total_sequences, 'sequences...', file=self.log_fh)
        file_reader = pyfastaq.sequences.file_reader(self.references_fa)
        best_score = 0
        best_seq_name = None
        for seq in file_reader:
            score = self._total_alignment_score(seq.id)
            print('Total alignment score for sequence', seq.id, 'is', score, file=self.log_fh)
            if score > best_score:
                best_score = score
                best_seq_name = seq.id

        print('\nBest sequence is', best_seq_name, 'with total alignment score of', best_score, file=self.log_fh)
        print(file=self.log_fh)
        return best_seq_name


    def best_seq(self, outfile):
        '''Finds the closest matchng sequence, writes it to a FASTA file, and returns it as a pyfastaq.sequences.Fasta object'''
        seq_name = self._get_best_seq_by_alignment_score()
        if seq_name is None:
            return None
        faidx.write_fa_subset([seq_name], self.references_fa, outfile, samtools_exe=self.samtools_exe, verbose=True, verbose_filehandle=self.log_fh)
        seqs = {}
        pyfastaq.tasks.file_to_dict(outfile, seqs)
        assert len(seqs) == 1
        return list(seqs.values())[0]

