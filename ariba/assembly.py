import os
import sys
import pyfastaq
import pymummer
import fermilite_ariba
from ariba import common, mapping, bam_parse, external_progs, ref_seq_chooser
import shlex

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
      spades_mode="wgs",
      spades_options=None,
      threads=1
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
        self.spades_mode = spades_mode
        self.spades_options = spades_options
        self.threads = threads

        if extern_progs is None:
            self.extern_progs = external_progs.ExternalProgs(using_spades=self.assembler == 'spades')
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
    def _check_spades_log_file(logfile):
        '''SPAdes can fail with a strange error. Stop everything if this happens'''
        f = pyfastaq.utils.open_file_read(logfile)

        for line in f:
            if line.startswith('== Error ==  system call for:') and line.rstrip().endswith('finished abnormally, err code: -7'):
                pyfastaq.utils.close(f)
                print('Error running SPAdes. Cannot continue. This is the error from the log file', logfile, '...', file=sys.stderr)
                print(line, file=sys.stderr)
                raise Error('Fatal error ("err code: -7") running spades. Cannot continue')

        pyfastaq.utils.close(f)
        return True

    def _assemble_with_spades(self):
        cwd = os.getcwd()
        self.assembled_ok = False
        try:
            try:
                os.chdir(self.working_dir)
            except:
                raise Error('Error chdir ' + self.working_dir)
            spades_exe = self.extern_progs.exe('spades')
            if not spades_exe:
                raise Error("Spades executable has not been found")
            spades_options = self.spades_options
            if spades_options is not None:
                spades_options = shlex.split(self.spades_options)
            if self.spades_mode == "rna":
                spades_options = ["--rna"] + (["-k","127"] if spades_options is None else spades_options)
                spades_out_seq_base = "transcripts.fasta"
            elif self.spades_mode == "sc":
                spades_options = ["--sc"] + (["-k", "33,55,77,99,127","--careful"] if spades_options is None else spades_options)
                spades_out_seq_base = "contigs.fasta"
            elif self.spades_mode == "wgs":
                spades_options = ["-k", "33,55,77,99,127","--careful"] if spades_options is None else spades_options
                spades_out_seq_base = "contigs.fasta"
            else:
                raise ValueError("Unknown spades_mode value: {}".format(self.spades_mode))
            asm_cmd = ['python3', spades_exe, "-t", str(self.threads), "--pe1-1", self.reads1, "--pe1-2", self.reads2, "-o", self.assembler_dir] + \
                spades_options
            asm_ok,err = common.syscall(asm_cmd, verbose=True, verbose_filehandle=self.log_fh, shell=False, allow_fail=True)
            if not asm_ok:
                print('Assembly finished with errors. These are the errors:', file=self.log_fh)
                print(err, file=self.log_fh)
                print('\nEnd of spades errors\n', file=self.log_fh)
            else:

                spades_log = os.path.join(self.assembler_dir, 'spades.log')
                if os.path.exists(spades_log):
                    self._check_spades_log_file(spades_log)

                    with open(spades_log) as f:
                        print('\n______________ SPAdes log ___________________\n', file=self.log_fh)
                        for line in f:
                            print(line.rstrip(), file=self.log_fh)
                        print('\n______________ End of SPAdes log _________________\n', file=self.log_fh)

                spades_warnings = os.path.join(self.assembler_dir, 'warnings.log')
                if os.path.exists(spades_warnings):
                    with open(spades_warnings) as f:
                        print('\n______________ SPAdes warnings ___________________\n', file=self.log_fh)
                        for line in f:
                            print(line.rstrip(), file=self.log_fh)
                        print('\n______________ End of SPAdes warnings _________________\n', file=self.log_fh)

                ## fermilight module generates contig names that look like `cluster_1.l15.c17.ctg.1` where 'cluster_1'==self.contig_name_prefix
                ## the whole structure of the contig name is expected in several places downstream where it is parsed into individual components.
                ## For example, it is parsed into to l and c parts in ref_seq_chooser (although the parts are not actually used).
                ## This is the code from fermilight module that generates the contig ID string:
                ## ofs << ">" << namePrefix << ".l" << overlap << ".c" << minCount << ".ctg." << i + 1 << '\n'
                ##
                ## We generate the same contig name structure here using dummy values for overlap and minCount, in order
                ## to avoid distrupting the downstream code.
                ## Note that the fermilight module generates multiple versions of the assembly on a grid of l and c values,
                ## and ref_seq_chooser then picks a single "best" (l,c) version based on coverage/identity of the nucmer
                ## alignment to the reference. Spades generates a single version of the assembly, so ref_seq_chooser
                ## can only pick that one version.

                spades_out_seq = os.path.join(self.assembler_dir,spades_out_seq_base)
                ## No need really to use general-purpose pyfastaq.sequences.file_reader here and pay performance cost for
                ## its multi-format line tests since we are only replacing the IDs in a pre-defined format
                if os.path.exists(spades_out_seq):
                    with open(spades_out_seq,"r") as inp, open(self.all_assembly_contigs_fa,"w") as out:
                        pref = self.contig_name_prefix
                        i_cont = 0
                        for line in inp:
                            if line.startswith(">"):
                                i_cont += 1
                                line = ">{}.l15.c17.ctg.{}\n".format(pref,i_cont)
                            out.write(line)
                        if i_cont > 0:
                            self.assembled_ok = True
            if self.clean:
                print('Deleting assembly directory', self.assembler_dir, file=self.log_fh)
                common.rmtree(self.assembler_dir)
        finally:
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
        if self.assembler == 'fermilite':
            self._assemble_with_fermilite()
        elif self.assembler == "spades":
            self._assemble_with_spades()
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
                threads=self.threads,
                sort=True,
                bowtie2=self.extern_progs.exe('bowtie2'),
                bowtie2_version=self.extern_progs.version('bowtie2'),
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
