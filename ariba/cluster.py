import os
import copy
from operator import itemgetter
import sys
import shutil
import pysam
import operator
import pyfastaq
import pymummer
from ariba import assembly, assembly_compare, common, mapping, bam_parse, flag, faidx, reference_data, report_line

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
      smalt_exe='smalt',
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
        self.references_fa = os.path.join(self.root_dir, 'referencees.fa')

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
        self.smalt_exe = smalt_exe
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


    def _initial_make_report_lines(self):
        '''Makes report lines. While they are being made, we discover if there were
        and non-synonymous variants. This affects the flag, which also gets updated
        by the function. To then fix the report lines, must run _update_flag_in_report_lines()'''
        self.report_lines = []
        if self.ref_sequence is not None:
            this_sequence_type = self.refdata.sequence_type(self.name)
        else:
            this_sequence_type = '.'
        print('self.name:', self.name)
        assert this_sequence_type is not None

        if not self.assembled_ok:
            #gene_name = 'NA' if self.gene is None else self.gene.id
            #gene_length = '.' if self.gene is None else len(self.gene)
            new_report_line = report_line.ReportLine(
                self.ref_sequence.id,
                self.refdata,
                self.status_flag,
                self.total_reads,
                self.name,
            )
            self.report_lines.append(new_report_line)
            return

        cov_per_contig = self._nucmer_hits_to_gene_cov_per_contig()
        samtools_variants = self._get_samtools_variants()
        ref_non_wild_variants = self.refdata.all_non_wild_type_variants(self.name)


        for contig in self.mummer_variants:
            # 1. get all mummer variants = self.mummer_variants[contig]
            # 2. get all known variants that this sequence has (ie not wild type)
            #       - if not found by mummer, then we want all variants where
            #         ref gene has the non-wild type

            # 3. get all samtools variants
            # 4. for each mummer variant:
            #        - check if same as a known variant
            #        - make report line
            # 5. check each known variant not already covered by 4.
            #        - make report line if we have any variants not wild type
            # 6. report any left over samtools variants
            # 7. If no report lines made yet, make non-variant report line


            for variants in self.mummer_variants[contig]:
                t = self._get_variant_effect(variants)
                if t is not None:
                    effect, new_bases = t
                    if effect != 'SYN':
                        self.status_flag.add('has_nonsynonymous_variants')

                    for v in variants:
                        depths = self._get_assembly_read_depths(contig, v.qry_start)
                        if depths is None:
                            # this happens with low coverage contigs. It can get assembled, but
                            # there are some bases that do not have reads mapped to them.
                            # If mummer called a variant at one of these, then we're looking
                            # for read dpeth where there is none.
                            print('Warning: could not get read depth info on contig "' + contig + '" at position ', str(v.qry_start + 1), 'from file', self.final_assembly_read_depths, file=sys.stderr)
                            print(' - a variant was called at this position using nucmer, but there is no read depth (probably a mapping artifact)', file=sys.stderr)
                            depths = ['.'] * 4

                        ref_base, alt_base, ref_counts, alt_counts = depths

                        self.report_lines.append([
                            self.ref_sequence.id,
                            self.status_flag.to_number(),
                            self.total_reads,
                            self.name,
                            len(self.ref_sequence),
                            cov_per_contig[contig],
                            self.percent_identities[contig],
                            pymummer.variant.var_types[v.var_type],
                            effect,
                            new_bases,
                            v.ref_start + 1,
                            v.ref_end + 1,
                            v.ref_base,
                            v.qry_name,
                            v.qry_length,
                            v.qry_start + 1,
                            v.qry_end + 1,
                            v.qry_base,
                            ref_counts,
                            alt_base,
                            alt_counts,
                        ])

                        if contig in samtools_variants and v.qry_start in samtools_variants[contig]:
                            del samtools_variants[contig][v.qry_start]
                            if len(samtools_variants[contig]) == 0:
                                del samtools_variants[contig]

            if contig in samtools_variants:
                for pos in samtools_variants[contig]:
                    ref_base, alt_base, ref_counts, alt_counts = samtools_variants[contig][pos]
                    self.report_lines.append(
                      [
                        self.ref_sequence.id,
                        self.status_flag.to_number(),
                        self.total_reads,
                        self.name,
                        len(self.ref_sequence),
                        cov_per_contig[contig],
                        self.percent_identities[contig],
                      ] + \
                      ['.'] * 6 + \
                      [
                        contig,
                        len(self.final_assembly[contig]),
                        pos + 1,
                        pos + 1,
                        ref_base,
                        ref_counts,
                        alt_base,
                        alt_counts
                      ]
                    )

        if len(self.report_lines) == 0:
            for contig in self.percent_identities:
                self.report_lines.append([
                    self.ref_sequence.id,
                    self.status_flag.to_number(),
                    self.total_reads,
                    self.name,
                    len(self.ref_sequence),
                    cov_per_contig[contig],
                    self.percent_identities[contig],
                  ] + \
                  ['.'] * 6 + [contig, len(self.final_assembly[contig])] + ['.'] * 6
                )

        self.report_lines.sort(key=itemgetter(0, 14, 15))


    def _update_flag_in_report_lines(self):
        '''This corrects the flag in all the report lines made by _initial_make_report_lines()'''
        flag_column = 1
        if self.status_flag.has('has_nonsynonymous_variants'):
            for line in self.report_lines:
                line[flag_column] = self.status_flag.to_number()


    def _make_report_lines(self):
        self._initial_make_report_lines()
        self._update_flag_in_report_lines()


    def _clean(self):
        print('Cleaning', self.root_dir, file=self.log_fh)

        if self.clean > 0:
            print('  rm -r', self.assembly_dir, file=self.log_fh)
            shutil.rmtree(self.assembly_dir)

        to_clean = [
            [
                'assembly.reads_mapped.unsorted.bam',
            ],
            [
                'assembly.fa.fai',
                'assembly.reads_mapped.bam.scaff',
                'assembly.reads_mapped.bam.soft_clipped',
                'assembly.reads_mapped.bam.unmapped_mates',
                'assembly_vs_gene.coords',
                'assembly_vs_gene.coords.snps',
                'genes.fa',
                'genes.fa.fai',
                'reads_1.fq',
                'reads_2.fq',
            ],
            [
                'assembly.fa.fai',
                'assembly.reads_mapped.bam',
                'assembly.reads_mapped.bam.vcf',
                'assembly_vs_gene.coords',
                'assembly_vs_gene.coords.snps',
            ]
        ]

        for i in range(self.clean + 1):
            for fname in to_clean[i]:
                fullname = os.path.join(self.root_dir, fname)
                if os.path.exists(fullname):
                    print('  rm', fname, file=self.log_fh)
                    os.unlink(fullname)


    def run(self):
        self.logfile = os.path.join(self.root_dir, 'log.txt')
        self.log_fh = pyfastaq.utils.open_file_write(self.logfile)
        self.assembly_compare_prefix = os.path.join(self.root_dir, 'assembly_compare')
        self.ref_sequence = self._choose_best_gene()

        if self.ref_sequence is None:
            self.assembled_ok = False
        else:
            self.assembly = assembly.Assembly(
              self.reads1,
              self.reads2,
              self.reference_fa,
              self.assembly_dir,
              self.final_assembly_fa,
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

            assembly.run()
            self.assembled_ok = assembly.assembled_ok

        if self.assembled_ok:
            self.assembly_compare = assembly_compare.AssemblyCompare(
              self.final_assembly_fa,
              self.assembly.sequences,
              self.final_assembly_bam,
              self.reference_fa,
              self.ref_sequence,
              self.assembly_compare_prefix,
              self.refdata,
              nucmer_min_id=self.nucmer_min_id,
              nucmer_min_len=self.nucmer_min_len,
              nucmer_breaklen=self.nucmer_breaklen,
              assembled_threshold=self.assembled_threshold,
              unique_threhsold=self.unique_threhsold,
            )

            self.assembly_compare.run()
            self.flag = self.assembly_compare.update_flag(self.flag)
            number_of_assembly_variants_in_ref_matches = samtools_vars.variants_in_coords(self.assembly_compare.assembly_match_coords())
            if number_of_assembly_variants_in_ref_matches > 1:
                self.status_flag.add('variants_suggest_collapsed_repeat')




            self._get_mummer_variants()
            self._filter_mummer_variants()
            self._update_flag_from_nucmer_file()
            self._make_assembly_vcf()
            self._get_vcf_variant_counts()
        else:
            self.status_flag.add('assembly_fail')

        self._make_report_lines()
        self._clean()
        pyfastaq.utils.close(self.log_fh)

        # This stops multiprocessing complaining with the error:
        # multiprocessing.pool.MaybeEncodingError: Error sending result: '[<ariba.cluster.Cluster object at 0x7ffa50f8bcd0>]'. Reason: 'TypeError("cannot serialize '_io.TextIOWrapper' object",)'
        self.log_fh = None
