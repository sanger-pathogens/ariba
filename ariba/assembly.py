import os
import shutil
import pyfastaq
import pymummer
from ariba import common, mapping, bam_parse

class Error (Exception): pass

class Assembly:
    def __init__(self,
      reads1,
      reads2,
      ref_fasta,
      working_dir,
      final_assembly_fa,
      final_assembly_bam,
      log_fh,
      scaff_name_prefix='scaffold',
      kmer=0,
      assembler='spades',
      bowtie2_exe='bowtie2',
      bowtie2_preset='very-sensitive-local',
      max_insert=1000,
      min_scaff_depth=10,
      nucmer_min_id=90,
      nucmer_min_len=50,
      nucmer_breaklen=50,
      spades_other_options=None,
      sspace_k=20,
      sspace_sd=0.4,
      reads_insert=500,
      samtools_exe='samtools',
      spades_exe='spades.py',
      sspace_exe='SSPACE_Basic_v2.0.pl',
      gapfiller_exe='GapFiller.pl',
    ):
        self.reads1 = os.path.abspath(reads1)
        self.reads2 = os.path.abspath(reads2)
        self.ref_fasta = os.path.abspath(ref_fasta)
        self.working_dir = os.path.abspath(working_dir)
        self.final_assembly_fa = os.path.abspath(final_assembly_fa)
        self.log_fh = log_fh
        self.scaff_name_prefix = scaff_name_prefix

        self.assembly_kmer = self._get_assembly_kmer(kmer, reads1, reads2)
        self.assembler = assembler
        self.bowtie2_exe = bowtie2_exe
        self.bowtie2_preset = bowtie2_preset
        self.max_insert = max_insert
        self.min_scaff_depth = min_scaff_depth
        self.nucmer_min_id = nucmer_min_id
        self.nucmer_min_len = nucmer_min_len
        self.nucmer_breaklen = nucmer_breaklen
        self.spades_other_options = spades_other_options
        self.sspace_k = sspace_k
        self.sspace_sd = sspace_sd
        self.reads_insert = reads_insert
        self.samtools_exe = samtools_exe
        self.spades_exe = spades_exe

        try:
            os.mkdir(self.working_dir)
        except:
            raise Error('Error mkdir ' + self.working_dir)

        self.sspace_exe = shutil.which(sspace_exe)
        if self.sspace_exe is None:
            self.gapfiller_exe = None
        else:
            self.sspace_exe = os.path.realpath(self.sspace_exe) # otherwise sspace dies loading packages
            self.gapfiller_exe = shutil.which(gapfiller_exe)
            if self.gapfiller_exe is not None:
                self.gapfiller_exe = os.path.realpath(self.gapfiller_exe) # otherwise gapfiller dies loading packages

        #self.assembly_dir = os.path.join(self.working_dir, 'Assembly')
        self.assembler_dir = os.path.join(self.working_dir, 'Assemble')
        self.assembly_contigs = os.path.join(self.working_dir, 'contigs.fa')
        self.scaffold_dir = os.path.join(self.working_dir, 'Scaffold')
        self.scaffolder_scaffolds = os.path.join(self.working_dir, 'scaffolds.fa')
        self.gapfill_dir = os.path.join(self.working_dir, 'Gapfill')
        self.gapfilled_scaffolds = os.path.join(self.working_dir, 'scaffolds.gapfilled.fa')
        #self.final_assembly_fa = os.path.join(self.assembly_dir,


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


    def _assemble_with_spades(self, unittest=False):
        cmd = ' '.join([
            self.spades_exe,
            '-1', self.reads1,
            '-2', self.reads2,
            '-o', self.assembler_dir,
            '-k', str(self.assembly_kmer),
            '--untrusted-contigs', self.ref_fasta,
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
            os.symlink(spades_contigs, os.path.basename(self.assembly_contigs))
        else:
            spades_errors_file = os.path.join(self.working_dir, 'spades_errors')
            with open(spades_errors_file, 'w') as f:
                print(err, file=f)
            f.close()

        os.chdir(cwd)


    def _scaffold_with_sspace(self):
        if not os.path.exists(self.assembly_contigs):
            raise Error('Cannot scaffold because contigs file not found: ' + self.assembly_contigs)

        try:
            os.mkdir(self.scaffold_dir)
        except:
            raise Error('Error mkdir '+  self.scaffold_dir)

        cwd = os.getcwd()

        if self.sspace_exe is None:
            os.chdir(self.assembly_dir)
            os.symlink(os.path.basename(self.assembly_contigs), os.path.basename(self.scaffolder_scaffolds))
            os.chdir(cwd)
            return

        os.chdir(self.scaffold_dir)
        lib_file = 'lib'
        with open(lib_file, 'w') as f:
            print('LIB', self.reads1, self.reads2, int(self.reads_insert), self.sspace_sd, 'FR', file=f)

        cmd = ' '.join([
            'perl', self.sspace_exe,
            '-k', str(self.sspace_k),
            '-l', lib_file,
            '-s', self.assembly_contigs
        ])

        sspace_scaffolds = os.path.abspath('standard_output.final.scaffolds.fasta')
        common.syscall(cmd, verbose=True, verbose_filehandle=self.log_fh)
        os.chdir(self.working_dir)
        os.symlink(os.path.relpath(sspace_scaffolds), os.path.basename(self.scaffolder_scaffolds))
        os.chdir(cwd)


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

        if self.gapfiller_exe is None or not self._has_gaps_to_fill(self.scaffolder_scaffolds):
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
            'perl', self.gapfiller_exe,
            '-l', lib_file,
            '-s', self.scaffolder_scaffolds
        ])

        gapfilled_scaffolds = os.path.join(self.gapfill_dir, 'standard_output', 'standard_output.gapfilled.final.fa')
        common.syscall(cmd, verbose=True, verbose_filehandle=self.log_fh)
        self._rename_scaffolds(gapfilled_scaffolds, self.gapfilled_scaffolds, self.scaff_name_prefix)
        os.chdir(cwd)


    @staticmethod
    def _fix_contig_orientation(contigs_fa, ref_fa, outfile, min_id=90, min_length=50, breaklen=50):
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
        self._assemble_with_spades()
        self.sequences = {}

        # double-check we got some contigs
        number_of_contigs = pyfastaq.tasks.count_sequences(self.assembly_contigs) if os.path.exists(self.assembly_contigs) else 0
        if number_of_contigs == 0:
            self.assembled_ok = False
            return
        else:
            self.assembled_ok = True

        if self.assembled_ok:
            self._scaffold_with_sspace()
            self._gap_fill_with_gapfiller()
            contigs_both_strands = self._fix_contig_orientation(self.gapfilled_scaffolds, self.final_assembly_fa)
            self.has_contigs_on_both_strands = len(contigs_both_strands) > 0
            pyfastaq.tasks.file_to_dict(self.final_assembly_fa, self.sequences)

            mapping.run_bowtie2(
                self.reads1,
                self.reads2,
                self.final_assembly_fa,
                self.final_assembly_bam[:-4],
                threads=1,
                sort=True,
                samtools=self.samtools_exe,
                bowtie2=self.bowtie2_exe,
                bowtie2_preset=self.bowtie2_preset,
                verbose=True,
                verbose_filehandle=self.log_fh
            )

            self.scaff_graph_ok = self._parse_bam(self.sequences, self.final_assembly_bam, self.min_scaff_depth, self.max_insert)
