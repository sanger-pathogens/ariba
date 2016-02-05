import os
import sys
import pyfastaq
from ariba import assembly, assembly_compare, assembly_variants, bam_parse, best_seq_chooser, flag, mapping, report, samtools_variants

class Error (Exception): pass

unittest=False

class Cluster:
    def __init__(self,
      root_dir,
      name,
      refdata,
      assembly_kmer=21,
      assembler='spades',
      max_insert=1000,
      min_scaff_depth=10,
      nucmer_min_id=90,
      nucmer_min_len=50,
      nucmer_breaklen=50,
      reads_insert=500,
      sspace_k=20,
      sspace_sd=0.4,
      threads=1,
      bcf_min_dp=10,
      bcf_min_dv=5,
      bcf_min_dv_over_dp=0.3,
      bcf_min_qual=20,
      assembled_threshold=0.95,
      unique_threshold=0.03,
      bcftools_exe='bcftools',
      gapfiller_exe='GapFiller.pl',
      samtools_exe='samtools',
      bowtie2_exe='bowtie2',
      bowtie2_preset='very-sensitive-local',
      spades_exe='spades.py',
      sspace_exe='SSPACE_Basic_v2.0.pl',
      spades_other_options=None,
      clean=1,
    ):

        self.root_dir = os.path.abspath(root_dir)
        if not os.path.exists(self.root_dir):
            raise Error('Directory ' + self.root_dir + ' not found. Cannot continue')

        self.name = name
        self.refdata = refdata
        self.assembly_kmer = assembly_kmer
        self.assembler = assembler
        self.sspace_k = sspace_k
        self.sspace_sd = sspace_sd
        self.reads_insert = reads_insert
        self.spades_exe = spades_exe
        self.sspace_exe = sspace_exe
        self.spades_other_options = spades_other_options
        self.gapfiller_exe = gapfiller_exe

        self.reads1 = os.path.join(self.root_dir, 'reads_1.fq')
        self.reads2 = os.path.join(self.root_dir, 'reads_2.fq')
        self.reference_fa = os.path.join(self.root_dir, 'reference.fa')
        self.references_fa = os.path.join(self.root_dir, 'references.fa')

        for fname in [self.reads1, self.reads2, self.references_fa]:
            if not os.path.exists(fname):
                raise Error('File ' + fname + ' not found. Cannot continue')


        self.ref_sequence = None

        self.max_insert = max_insert
        self.min_scaff_depth = min_scaff_depth

        self.nucmer_min_id = nucmer_min_id
        self.nucmer_min_len = nucmer_min_len
        self.nucmer_breaklen = nucmer_breaklen
        self.assembly_vs_gene_coords = os.path.join(self.root_dir, 'assembly_vs_gene.coords')

        self.bcf_min_dp = bcf_min_dp
        self.bcf_min_dv = bcf_min_dv
        self.bcf_min_dv_over_dp = bcf_min_dv_over_dp
        self.bcf_min_qual = bcf_min_qual

        self.bcftools_exe = bcftools_exe
        self.samtools_exe = samtools_exe
        self.bowtie2_exe = bowtie2_exe
        self.bowtie2_preset = bowtie2_preset

        self.threads = threads
        self.assembled_threshold = assembled_threshold
        self.unique_threshold = unique_threshold
        self.status_flag = flag.Flag()
        self.flag_file = os.path.join(self.root_dir, 'flag')
        self.clean = clean

        self.assembly_dir = os.path.join(self.root_dir, 'Assembly')
        self.final_assembly_fa = os.path.join(self.root_dir, 'assembly.fa')
        self.final_assembly_bam = os.path.join(self.root_dir, 'assembly.reads_mapped.bam')
        self.final_assembly_read_depths = os.path.join(self.root_dir, 'assembly.reads_mapped.bam.read_depths.gz')
        self.final_assembly_vcf = os.path.join(self.root_dir, 'assembly.reads_mapped.bam.vcf')
        self.samtools_vars_prefix = self.final_assembly_bam
        self.assembly_compare_prefix = os.path.join(self.root_dir, 'assembly_compare')

        self.mummer_variants = {}
        self.variant_depths = {}
        self.percent_identities = {}
        self.total_reads = self._count_reads(self.reads1, self.reads2)

        # The log filehandle self.log_fh is set at the start of the run() method.
        # Lots of other methods use self.log_fh. But for unit testing, run() isn't
        # run. So we need to set this to something for unit testing.
        # On the other hand, setting it here breaks a real run of ARIBA because
        # multiprocessing complains with the error:
        # TypeError: cannot serialize '_io.TextIOWrapper' object.
        # Hence the following two lines...
        if unittest:
            self.log_fh = sys.stdout


    @staticmethod
    def _count_reads(reads1, reads2):
        count1 = pyfastaq.tasks.count_sequences(reads1)
        count2 = pyfastaq.tasks.count_sequences(reads2)
        assert(count1 == count2)
        return count1 + count2


    def run(self):
        self.logfile = os.path.join(self.root_dir, 'log.txt')
        self.log_fh = pyfastaq.utils.open_file_write(self.logfile)

        seq_chooser = best_seq_chooser.BestSeqChooser(
            self.reads1,
            self.reads2,
            self.references_fa,
            self.log_fh,
            samtools_exe=self.samtools_exe,
            bowtie2_exe=self.bowtie2_exe,
            bowtie2_preset=self.bowtie2_preset,
            threads=1,
        )
        self.ref_sequence = seq_chooser.best_seq(self.reference_fa)

        if self.ref_sequence is None:
            self.status_flag.add('ref_seq_choose_fail')
            self.assembled_ok = False
        else:
            self.ref_sequence_type = self.refdata.sequence_type(self.ref_sequence.id)
            assert self.ref_sequence_type is not None
            self.assembly = assembly.Assembly(
              self.reads1,
              self.reads2,
              self.reference_fa,
              self.assembly_dir,
              self.final_assembly_fa,
              self.final_assembly_bam,
              self.log_fh,
              scaff_name_prefix=self.ref_sequence.id,
              kmer=self.assembly_kmer,
              assembler=self.assembler,
              spades_other_options=self.spades_other_options,
              sspace_k=self.sspace_k,
              sspace_sd=self.sspace_sd,
              reads_insert=self.reads_insert,
              spades_exe=self.spades_exe,
              sspace_exe=self.sspace_exe,
              gapfiller_exe=self.gapfiller_exe,
            )

            self.assembly.run()
            self.assembled_ok = self.assembly.assembled_ok

        if self.assembled_ok:
            print('\nAssembly was successful\n', file=self.log_fh)

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
            bam_parser = bam_parse.Parser(self.final_assembly_bam, self.assembly.sequences)
            bam_parser.parse()
            if not bam_parser.scaff_graph_is_consistent(self.min_scaff_depth, self.max_insert):
                self.status_flag.add('scaffold_graph_bad')

            self.assembly_compare = assembly_compare.AssemblyCompare(
              self.final_assembly_fa,
              self.assembly.sequences,
              self.reference_fa,
              self.ref_sequence,
              self.assembly_compare_prefix,
              self.refdata,
              nucmer_min_id=self.nucmer_min_id,
              nucmer_min_len=self.nucmer_min_len,
              nucmer_breaklen=self.nucmer_breaklen,
              assembled_threshold=self.assembled_threshold,
              unique_threshold=self.unique_threshold,
            )
            self.assembly_compare.run()
            self.status_flag = self.assembly_compare.update_flag(self.status_flag)

            nucmer_hits_to_ref = assembly_compare.AssemblyCompare.nucmer_hits_to_ref_coords(self.assembly_compare.nucmer_hits)
            assembly_variants_obj = assembly_variants.AssemblyVariants(self.refdata, self.assembly_compare.nucmer_snps_file)
            self.assembly_variants = assembly_variants_obj.get_variants(self.ref_sequence.id, nucmer_hits_to_ref)

            self.samtools_vars = samtools_variants.SamtoolsVariants(
                self.final_assembly_fa,
                self.final_assembly_bam,
                self.samtools_vars_prefix,
                log_fh=self.log_fh,
                samtools_exe=self.samtools_exe,
                bcftools_exe=self.bcftools_exe,
                bcf_min_dp=self.bcf_min_dp,
                bcf_min_dv=self.bcf_min_dv,
                bcf_min_dv_over_dp=self.bcf_min_dv_over_dp,
                bcf_min_qual=self.bcf_min_qual,
            )
            self.samtools_vars.run()
            if self.samtools_vars.variants_in_coords(self.assembly_compare.assembly_match_coords(), self.samtools_vars.vcf_file):
                self.status_flag.add('variants_suggest_collapsed_repeat')
        else:
            print('\nAssembly failed\n', file=self.log_fh)
            self.status_flag.add('assembly_fail')

        self.report_lines = report.report_lines(self)
        #self._clean()
        pyfastaq.utils.close(self.log_fh)

        # This stops multiprocessing complaining with the error:
        # multiprocessing.pool.MaybeEncodingError: Error sending result: '[<ariba.cluster.Cluster object at 0x7ffa50f8bcd0>]'. Reason: 'TypeError("cannot serialize '_io.TextIOWrapper' object",)'
        self.log_fh = None

