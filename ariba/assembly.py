import os
import sys
import shutil
import pyfastaq
import pymummer
import fermilite_ariba
from ariba import common, faidx, mapping, bam_parse, external_progs, ref_seq_chooser

class Error (Exception): pass

class Assembly:
    def __init__(self,
      reads1,
      reads2,
      ref_fasta,
      ref_fastas,
      working_dir,
      final_assembly_fa,
      final_assembly_bam,
      log_fh,
      all_reference_fasta,
      contig_name_prefix='ctg',
      assembler='fermilite',
      max_insert=1000,
      min_scaff_depth=10,
      min_scaff_length=50,
      nucmer_min_id=90,
      nucmer_min_len=20,
      nucmer_breaklen=200,
      extern_progs=None,
      clean=True,
    ):
        self.reads1 = os.path.abspath(reads1)
        self.reads2 = os.path.abspath(reads2)
        self.ref_fasta = os.path.abspath(ref_fasta)
        self.ref_fastas = os.path.abspath(ref_fastas)
        self.working_dir = os.path.abspath(working_dir)
        self.final_assembly_fa = os.path.abspath(final_assembly_fa)
        self.final_assembly_bam = os.path.abspath(final_assembly_bam)
        self.log_fh = log_fh
        self.all_reference_fasta = os.path.abspath(all_reference_fasta)
        self.contig_name_prefix = contig_name_prefix

        self.ref_seq_name = None
        self.assembler = assembler
        self.max_insert = max_insert
        self.min_scaff_depth = min_scaff_depth
        self.min_scaff_length = min_scaff_length
        self.nucmer_min_id = nucmer_min_id
        self.nucmer_min_len = nucmer_min_len
        self.nucmer_breaklen = nucmer_breaklen
        self.clean = clean

        if extern_progs is None:
            self.extern_progs = external_progs.ExternalProgs()
        else:
            self.extern_progs = extern_progs

        try:
            os.mkdir(self.working_dir)
        except:
            raise Error('Error mkdir ' + self.working_dir)

        self.assembler_dir = os.path.join(self.working_dir, 'Assemble')
        self.all_assembly_contigs_fa = os.path.join(self.working_dir, 'debug_all_contigs.fa')
        self.best_assembly_fa = os.path.join(self.working_dir, 'debug_best_assembly.fa')
        self.final_assembly_fa = os.path.abspath(final_assembly_fa)


    @staticmethod
    def _run_fermilite(reads_in, fasta_out, log_out, name_prefix):
        return fermilite_ariba.fermilite_ariba(reads_in, fasta_out, log_out, name_prefix)


    def _assemble_with_fermilite(self):
        cwd = os.getcwd()
        try:
            os.chdir(self.working_dir)
        except:
            raise Error('Error chdir ' + self.working_dir)

        interleaved_reads = 'reads.fq'
        pyfastaq.tasks.interleave(self.reads1, self.reads2, interleaved_reads)
        fermilite_log = self.all_assembly_contigs_fa + '.log'
        got_from_fermilite = self._run_fermilite(interleaved_reads, self.all_assembly_contigs_fa, fermilite_log, self.contig_name_prefix)
        os.unlink(interleaved_reads)
        if os.path.exists(fermilite_log):
            with open(fermilite_log) as f:
                for line in f:
                    print(line, end='', file=self.log_fh)

            os.unlink(fermilite_log)

        self.assembled_ok = (got_from_fermilite == 0)
        os.chdir(cwd)


    @staticmethod
    def _fix_contig_orientation(contigs_fa, ref_fa, outfile, min_id=90, min_length=20, breaklen=200):
        '''Changes orientation of each contig to match the reference, when possible.
           Returns a set of names of contigs that had hits in both orientations to the reference'''
        if not os.path.exists(contigs_fa):
            raise Error('Cannot fix orientation of assembly contigs because file not found: ' + contigs_fa)

        tmp_coords = os.path.join(outfile + '.tmp.rename.coords')
        pymummer.nucmer.Runner(
            ref_fa,
            contigs_fa,
            tmp_coords,
            min_id=min_id,
            min_length=min_length,
            breaklen=breaklen,
            maxmatch=True,
        ).run()

        to_revcomp = set()
        not_revcomp = set()
        file_reader = pymummer.coords_file.reader(tmp_coords)
        for hit in file_reader:
            if hit.on_same_strand():
                not_revcomp.add(hit.qry_name)
            else:
                to_revcomp.add(hit.qry_name)

        os.unlink(tmp_coords)
        in_both = to_revcomp.intersection(not_revcomp)

        f = pyfastaq.utils.open_file_write(outfile)
        seq_reader = pyfastaq.sequences.file_reader(contigs_fa)
        for seq in seq_reader:
            if seq.id in to_revcomp and seq.id not in in_both:
                seq.revcomp()
            print(seq, file=f)
        pyfastaq.utils.close(f)

        return in_both


    @staticmethod
    def _parse_bam(sequences, bam, min_scaff_depth, max_insert):
        if not os.path.exists(bam):
            raise Error('File not found: ' + bam)

        bam_parser = bam_parse.Parser(bam, sequences)
        bam_parser.parse()
        bam_parser.write_files(bam)
        return bam_parser.scaff_graph_is_consistent(min_scaff_depth, max_insert)


    def run(self):
        self._assemble_with_fermilite()
        print('Finished running assemblies', flush=True, file=self.log_fh)
        self.sequences = {}

        # double-check we got some contigs
        number_of_contigs = pyfastaq.tasks.count_sequences(self.all_assembly_contigs_fa) if os.path.exists(self.all_assembly_contigs_fa) else 0
        if number_of_contigs == 0:
            self.assembled_ok = False
            # This is to make this object picklable, to keep multithreading happy
            self.log_fh = None
            return
        else:
            self.assembled_ok = True

        if self.assembled_ok:
            ref_chooser = ref_seq_chooser.RefSeqChooser(
                self.ref_fastas,
                self.all_reference_fasta,
                self.all_assembly_contigs_fa,
                self.best_assembly_fa,
                self.log_fh,
                nucmer_min_id=self.nucmer_min_id,
                nucmer_min_len=self.nucmer_min_len,
                nucmer_breaklen=self.nucmer_breaklen,
            )
            ref_chooser.run()

            if ref_chooser.closest_ref_from_all_refs is None:
                print('Could not find match to reference sequences', file=self.log_fh)
                self.ref_seq_name = None
                self.log_fh = None
                return
            elif not ref_chooser.closest_ref_is_in_cluster:
                print('Closest reference', ref_chooser.closest_ref_from_all_refs, 'was not in cluster', file=self.log_fh)
                self.ref_seq_name = None
                self.log_fh = None
                return
            else:
                assert ref_chooser.closest_ref_from_all_refs is not None
                self.ref_seq_name = ref_chooser.closest_ref_from_all_refs

            print('Closest reference sequence:', self.ref_seq_name, file=self.log_fh)

            file_reader = pyfastaq.sequences.file_reader(self.ref_fastas)
            for ref_seq in file_reader:
                if self.ref_seq_name == ref_seq.id:
                    f_out = pyfastaq.utils.open_file_write(self.ref_fasta)
                    print(ref_seq, file=f_out)
                    pyfastaq.utils.close(f_out)
                    break

            contigs_both_strands = self._fix_contig_orientation(self.best_assembly_fa, self.ref_fasta, self.final_assembly_fa, min_id=self.nucmer_min_id, min_length=self.nucmer_min_len, breaklen=self.nucmer_breaklen)
            self.has_contigs_on_both_strands = len(contigs_both_strands) > 0
            pyfastaq.tasks.file_to_dict(self.final_assembly_fa, self.sequences)

            mapping.run_bowtie2(
                self.reads1,
                self.reads2,
                self.final_assembly_fa,
                self.final_assembly_bam[:-4],
                threads=1,
                sort=True,
                bowtie2=self.extern_progs.exe('bowtie2'),
                verbose=True,
                verbose_filehandle=self.log_fh
            )

            self.scaff_graph_ok = self._parse_bam(self.sequences, self.final_assembly_bam, self.min_scaff_depth, self.max_insert)
            print('Scaffolding graph is OK:', self.scaff_graph_ok, file=self.log_fh)

            if self.clean:
                for suffix in ['soft_clipped', 'unmapped_mates', 'scaff']:
                    filename = self.final_assembly_bam + '.' + suffix
                    print('Deleting file', filename, file=self.log_fh)
                    os.unlink(filename)


        # This is to make this object picklable, to keep multithreading happy
        self.log_fh = None
