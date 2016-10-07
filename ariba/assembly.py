import os
import sys
import shutil
import pyfastaq
import pymummer
import fermilite_ariba
from ariba import common, faidx, mapping, bam_parse, external_progs, mash

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
      mash_reference_fasta,
      scaff_name_prefix='scaffold',
      kmer=0,
      assembler='fermilite',
      max_insert=1000,
      min_scaff_depth=10,
      min_scaff_length=50,
      nucmer_min_id=90,
      nucmer_min_len=20,
      nucmer_breaklen=200,
      spades_other_options=None,
      sspace_k=20,
      sspace_sd=0.4,
      reads_insert=500,
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
        self.mash_reference_fasta = os.path.abspath(mash_reference_fasta)
        self.scaff_name_prefix = scaff_name_prefix

        self.ref_seq_name = None
        self.assembly_kmer = self._get_assembly_kmer(kmer, reads1, reads2)
        self.assembler = assembler
        self.max_insert = max_insert
        self.min_scaff_depth = min_scaff_depth
        self.min_scaff_length = min_scaff_length
        self.nucmer_min_id = nucmer_min_id
        self.nucmer_min_len = nucmer_min_len
        self.nucmer_breaklen = nucmer_breaklen
        self.spades_other_options = spades_other_options
        self.sspace_k = sspace_k
        self.sspace_sd = sspace_sd
        self.reads_insert = reads_insert
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
        self.assembly_contigs = os.path.join(self.working_dir, 'contigs.fa')
        self.scaffold_dir = os.path.join(self.working_dir, 'Scaffold')
        self.scaffolder_scaffolds = os.path.join(self.working_dir, 'scaffolds.fa')
        self.gapfill_dir = os.path.join(self.working_dir, 'Gapfill')
        self.gapfilled_scaffolds = os.path.join(self.working_dir, 'scaffolds.gapfilled.fa')
        self.gapfilled_length_filtered = os.path.join(self.working_dir, 'scaffolds.gapfilled.length_filtered.fa')
        self.mash_dist_file = os.path.join(self.working_dir, 'mash_dist_to_ref_seqs')


    @staticmethod
    def _get_assembly_kmer(k, reads1, reads2):
        '''If the kmer not given, uses 2/3 of the mean read length (using first 1000 forward and first 1000 reverse reads)'''
        if k == 0:
            read_length1 = pyfastaq.tasks.mean_length(reads1, limit=1000)
            read_length2 = pyfastaq.tasks.mean_length(reads2, limit=1000)
            assembly_kmer = round( (read_length1 + read_length2) / 3)
            if assembly_kmer % 2 == 0:
                assembly_kmer += 1
        else:
            assembly_kmer = k

        return assembly_kmer


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


    @staticmethod
    def _run_fermilite(reads_in, fasta_out, log_out):
        return fermilite_ariba.fermilite_ariba(reads_in, fasta_out, log_out)


    def _assemble_with_fermilite(self):
        cwd = os.getcwd()
        try:
            os.chdir(self.working_dir)
        except:
            raise Error('Error chdir ' + self.working_dir)

        interleaved_reads = 'reads.fq'
        pyfastaq.tasks.interleave(self.reads1, self.reads2, interleaved_reads)
        fermilite_log = self.assembly_contigs + '.log'
        got_from_fermilite = self._run_fermilite(interleaved_reads, self.assembly_contigs, fermilite_log)
        os.unlink(interleaved_reads)
        if os.path.exists(fermilite_log):
            with open(fermilite_log) as f:
                for line in f:
                    print(line, end='', file=self.log_fh)

            os.unlink(fermilite_log)

        self.assembled_ok = (got_from_fermilite == 0)
        os.chdir(cwd)


    def _assemble_with_spades(self, unittest=False):
        cmd = ' '.join([
            self.extern_progs.exe('spades'),
            '-1', self.reads1,
            '-2', self.reads2,
            '-o', self.assembler_dir,
            '-k', str(self.assembly_kmer),
            '--threads 1', # otherwise defaults to 16!
        ])
        if self.spades_other_options is not None:
            cmd += ' ' + self.spades_other_options

        cwd = os.getcwd()
        try:
            os.chdir(self.working_dir)
        except:
            raise Error('Error chdir ' + self.working_dir)
        spades_contigs = os.path.join(os.path.split(self.assembler_dir)[1], 'scaffolds.fasta')

        if unittest:
            os.mkdir(self.assembler_dir)
            open(spades_contigs, 'w').close()
            self.assembled_ok = True
        else:
            self.assembled_ok, err = common.syscall(cmd, verbose=True, allow_fail=True, verbose_filehandle=self.log_fh, print_errors=False)
        if self.assembled_ok:
            os.rename(spades_contigs, os.path.basename(self.assembly_contigs))
        else:
            print('Assembly finished with errors. These are the errors:', file=self.log_fh)
            print(err, file=self.log_fh)
            print('\nEnd of spades errors\n', file=self.log_fh)

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

        os.chdir(cwd)

        if self.clean:
            print('Deleting assembly directory', self.assembler_dir, file=self.log_fh)
            shutil.rmtree(self.assembler_dir)


    def _scaffold_with_sspace(self):
        if not os.path.exists(self.assembly_contigs):
            raise Error('Cannot scaffold because contigs file not found: ' + self.assembly_contigs)

        try:
            os.mkdir(self.scaffold_dir)
        except:
            raise Error('Error mkdir '+  self.scaffold_dir)

        cwd = os.getcwd()

        #if self.extern_progs.exe('sspace') is None:
        if True:  # no longer use sspace, but leave the option here just in case
            os.chdir(self.working_dir)
            os.symlink(self.assembly_contigs, os.path.basename(self.scaffolder_scaffolds))
            os.chdir(cwd)
            return

        os.chdir(self.scaffold_dir)
        lib_file = 'lib'
        with open(lib_file, 'w') as f:
            print('LIB', self.reads1, self.reads2, int(self.reads_insert), self.sspace_sd, 'FR', file=f)

        cmd = ' '.join([
            'perl', self.extern_progs.exe('sspace'),
            '-k', str(self.sspace_k),
            '-l', lib_file,
            '-s', self.assembly_contigs
        ])

        common.syscall(cmd, verbose=True, verbose_filehandle=self.log_fh)
        sspace_scaffolds = os.path.abspath('standard_output.final.scaffolds.fasta')
        sspace_log = os.path.abspath('standard_output.logfile.txt')
        with open(sspace_log) as f:
            print('\n_______________ SSPACE log __________________\n', file=self.log_fh)
            for line in f:
                print(line.rstrip(), file=self.log_fh)
            print('_______________ End of SSPACE log __________________\n', file=self.log_fh)

        os.rename(sspace_scaffolds, self.scaffolder_scaffolds)
        os.chdir(cwd)

        if self.clean:
            print('Deleting scaffolding directory', self.scaffold_dir, file=self.log_fh)
            shutil.rmtree(self.scaffold_dir)


    @staticmethod
    def _has_gaps_to_fill(filename):
        seq_reader = pyfastaq.sequences.file_reader(filename)
        for seq in seq_reader:
            if 'n' in seq.seq or 'N' in seq.seq:
                return True
        return False


    @staticmethod
    def _rename_scaffolds(infile, outfile, prefix):
        freader = pyfastaq.sequences.file_reader(infile)
        f_out = pyfastaq.utils.open_file_write(outfile)
        i = 1
        for scaff in freader:
            scaff.id = prefix + '.scaffold.' + str(i)
            i += 1
            print(scaff, file=f_out)
        pyfastaq.utils.close(f_out)


    def _gap_fill_with_gapfiller(self):
        if not os.path.exists(self.scaffolder_scaffolds):
            raise Error('Cannot gap fill because scaffolds file not found: ' + self.scaffolder_scaffolds)

        cwd = os.getcwd()

        #if self.extern_progs.exe('gapfiller') is None or not self._has_gaps_to_fill(self.scaffolder_scaffolds):
        if True: # no longer use gapfiller, but leave the option in, just in case
            self._rename_scaffolds(self.scaffolder_scaffolds, self.gapfilled_scaffolds, self.scaff_name_prefix)
            return

        try:
            os.mkdir(self.gapfill_dir)
        except:
            raise Error('Error mkdir '+  self.gapfill_dir)

        os.chdir(self.gapfill_dir)
        lib_file = 'lib'
        with open(lib_file, 'w') as f:
            print('LIB', 'bwa', self.reads1, self.reads2, self.reads_insert, self.sspace_sd, 'FR', file=f)

        cmd = ' '.join([
            'perl', self.extern_progs.exe('gapfiller'),
            '-l', lib_file,
            '-s', self.scaffolder_scaffolds
        ])

        gapfilled_scaffolds = os.path.join(self.gapfill_dir, 'standard_output', 'standard_output.gapfilled.final.fa')
        common.syscall(cmd, verbose=True, verbose_filehandle=self.log_fh)
        self._rename_scaffolds(gapfilled_scaffolds, self.gapfilled_scaffolds, self.scaff_name_prefix)
        os.chdir(cwd)
        if self.clean:
            print('Deleting GapFiller directory', self.gapfill_dir, file=self.log_fh)
            shutil.rmtree(self.gapfill_dir)


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
        self.sequences = {}

        # double-check we got some contigs
        number_of_contigs = pyfastaq.tasks.count_sequences(self.assembly_contigs) if os.path.exists(self.assembly_contigs) else 0
        if number_of_contigs == 0:
            self.assembled_ok = False
            # This is to make this object picklable, to keep multithreading happy
            self.log_fh = None
            return
        else:
            self.assembled_ok = True

        if self.assembled_ok:
            self._scaffold_with_sspace()
            self._gap_fill_with_gapfiller()

            pyfastaq.tasks.filter(self.gapfilled_scaffolds, self.gapfilled_length_filtered, minlength=self.min_scaff_length)
            if pyfastaq.tasks.count_sequences(self.gapfilled_length_filtered) == 0:
                self.assembled_ok = False
                # This is to make this object picklable, to keep multithreading happy
                self.log_fh = None
                return

            masher = mash.Masher(self.mash_reference_fasta, self.gapfilled_length_filtered, self.log_fh, self.extern_progs)
            self.ref_seq_name = masher.run(self.mash_dist_file)
            if self.ref_seq_name is None:
                print('Could not determine closest reference sequence', file=self.log_fh)
                self.log_fh = None
                return

            file_reader = pyfastaq.sequences.file_reader(self.ref_fastas)
            for ref_seq in file_reader:
                if self.ref_seq_name == ref_seq.id:
                    f_out = pyfastaq.utils.open_file_write(self.ref_fasta)
                    print(ref_seq, file=f_out)
                    pyfastaq.utils.close(f_out)
                    break
            else:
                print('Closest reference sequence ', self.ref_seq_name, ' does not belong to this cluster', file=self.log_fh)
                self.ref_seq_name = None
                self.log_fh = None
                return

            print('Closest reference sequence according to mash: ', self.ref_seq_name, file=self.log_fh)

            contigs_both_strands = self._fix_contig_orientation(self.gapfilled_length_filtered, self.ref_fasta, self.final_assembly_fa, min_id=self.nucmer_min_id, min_length=self.nucmer_min_len, breaklen=self.nucmer_breaklen)
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
