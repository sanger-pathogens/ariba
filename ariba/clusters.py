import os
import copy
import itertools
import sys
import shutil
import openpyxl
import multiprocessing
import pysam
import pyfastaq
from ariba import cdhit, cluster, common, mapping, histogram, faidx, report, report_filter

class Error (Exception): pass


def _run_cluster(obj, verbose):
    if verbose:
        print('Start running cluster', obj.name, 'in directory', obj.root_dir, flush=True)
    obj.run()
    if verbose:
        print('Finished running cluster', obj.name, 'in directory', obj.root_dir, flush=True)
    return obj


class Clusters:
    def __init__(self,
      refdata,
      reads_1,
      reads_2,
      outdir,
      extern_progs,
      assembly_kmer=21,
      threads=1,
      verbose=False,
      assembler='spades',
      spades_other=None,
      max_insert=1000,
      min_scaff_depth=10,
      nucmer_min_id=90,
      nucmer_min_len=50,
      nucmer_breaklen=50,
      assembled_threshold=0.95,
      unique_threshold=0.03,
      bowtie2_preset='very-sensitive-local',
      cdhit_seq_identity_threshold=0.9,
      cdhit_length_diff_cutoff=0.9,
      run_cd_hit=True,
      clean=1,
    ):
        self.refdata = refdata
        self.reads_1 = os.path.abspath(reads_1)
        self.reads_2 = os.path.abspath(reads_2)
        self.outdir = os.path.abspath(outdir)
        self.extern_progs = extern_progs
        self.clusters_outdir = os.path.join(self.outdir, 'Clusters')
        self.clean = clean

        self.assembler = assembler
        assert self.assembler in ['spades']
        self.assembly_kmer = assembly_kmer
        self.spades_other = spades_other

        self.refdata_files_prefix = os.path.join(self.outdir, 'refdata')
        self.cdhit_files_prefix = os.path.join(self.outdir, 'cdhit')
        self.cdhit_cluster_representatives_fa = self.cdhit_files_prefix + '.cluster_representatives.fa'
        self.cluster_ids = {}
        self.bam_prefix = self.cdhit_cluster_representatives_fa + '.map_reads'
        self.bam = self.bam_prefix + '.bam'
        self.report_file_all_tsv = os.path.join(self.outdir, 'report.all.tsv')
        self.report_file_all_xls = os.path.join(self.outdir, 'report.all.xls')
        self.report_file_filtered_prefix = os.path.join(self.outdir, 'report')
        self.catted_assembled_seqs_fasta = os.path.join(self.outdir, 'assembled_seqs.fa')
        self.threads = threads
        self.verbose = verbose

        self.max_insert = max_insert
        self.bowtie2_preset = bowtie2_preset

        self.insert_hist_bin = 10
        self.insert_hist = histogram.Histogram(self.insert_hist_bin)
        self.insert_size = None
        self.insert_sspace_sd = None
        self.insert_proper_pair_max = None

        self.min_scaff_depth = min_scaff_depth
        self.nucmer_min_id = nucmer_min_id
        self.nucmer_min_len = nucmer_min_len
        self.nucmer_breaklen = nucmer_breaklen

        self.assembled_threshold = assembled_threshold
        self.unique_threshold = unique_threshold

        self.cluster_to_dir = {}  # gene name -> abs path of cluster directory
        self.clusters = {}        # gene name -> Cluster object

        self.cdhit_seq_identity_threshold = cdhit_seq_identity_threshold
        self.cdhit_length_diff_cutoff = cdhit_length_diff_cutoff
        self.run_cd_hit = run_cd_hit

        for d in [self.outdir, self.clusters_outdir]:
            try:
                os.mkdir(d)
            except:
                raise Error('Error mkdir ' + d)


    def _run_cdhit(self):
        self.cluster_ids = self.refdata.cluster_with_cdhit(
            self.refdata_files_prefix + '.01.check_variants',
            self.cdhit_files_prefix,
            seq_identity_threshold=self.cdhit_seq_identity_threshold,
            threads=self.threads,
            length_diff_cutoff=self.cdhit_length_diff_cutoff,
            nocluster=not self.run_cd_hit,
            verbose=self.verbose,
        )


    def _map_reads_to_clustered_genes(self):
        mapping.run_bowtie2(
            self.reads_1,
            self.reads_2,
            self.cdhit_cluster_representatives_fa,
            self.bam_prefix,
            threads=self.threads,
            samtools=self.extern_progs.exe('samtools'),
            bowtie2=self.extern_progs.exe('bowtie2'),
            bowtie2_preset=self.bowtie2_preset,
            verbose=self.verbose,
            remove_both_unmapped=True,
        )


    def _sam_to_fastq(self, s):
        name = s.qname
        if s.is_read1:
            name += '/1'
        elif s.is_read2:
            name += '/2'
        else:
            raise Error('Read ' + name + ' must be first or second of pair according to flag. Cannot continue')

        seq = pyfastaq.sequences.Fastq(name, common.decode(s.seq), common.decode(s.qual))
        if s.is_reverse:
            seq.revcomp()

        return seq


    def _sam_pair_to_insert(self, s1, s2):
        if s1.is_unmapped or s2.is_unmapped or (s1.tid != s2.tid) or (s1.is_reverse == s2.is_reverse):
            return None

        # If here, reads are both mapped to the same ref, and in opposite orientations
        if s1.is_reverse:
            end = s1.reference_end - 1
            start = s2.reference_start
        else:
            end = s2.reference_end - 1
            start = s1.reference_start

        if start < end:
            return end - start + 1
        else:
            return None


    def _bam_to_clusters_reads(self):
        '''Sets up Cluster directories (one for each gene that has reads that mapped to it), writes reads fwd and rev files. Also gathers histogram data of insert size'''
        filehandles_1 = {} # gene name -> filehandle of fwd_reads
        filehandles_2 = {} # gene name -> filehandle of rev_reads
        sam_reader = pysam.Samfile(self.bam, "rb")
        sam1 = None

        for s in sam_reader.fetch(until_eof=True):
            if sam1 is None:
                sam1 = s
                continue

            ref_seqs = set()
            if not s.is_unmapped:
                ref_seqs.add(sam_reader.getrname(s.tid))
            if not sam1.is_unmapped:
                ref_seqs.add(sam_reader.getrname(sam1.tid))

            read1 = self._sam_to_fastq(sam1)
            read2 = self._sam_to_fastq(s)
            if read1.id.endswith('/2'):
                read1, read2 = read2, read1

            insert = self._sam_pair_to_insert(s, sam1)
            if insert is not None:
                self.insert_hist.add(insert)

            for ref in ref_seqs:
                if ref not in self.cluster_to_dir:
                    assert ref not in filehandles_1
                    assert ref not in filehandles_2

                    new_dir = os.path.join(self.clusters_outdir, ref)
                    try:
                        os.mkdir(new_dir)
                    except:
                        raise Error('Error mkdir ' + new_dir)

                    self.cluster_to_dir[ref] = new_dir
                    filehandles_1[ref] = pyfastaq.utils.open_file_write(os.path.join(new_dir, 'reads_1.fq'))
                    filehandles_2[ref] = pyfastaq.utils.open_file_write(os.path.join(new_dir, 'reads_2.fq'))
                    if self.verbose:
                        print('New cluster with reads that hit:', ref, flush=True)

                print(read1, file=filehandles_1[ref])
                print(read2, file=filehandles_2[ref])

            sam1 = None

        for ref in filehandles_1:
            pyfastaq.utils.close(filehandles_1[ref])
            pyfastaq.utils.close(filehandles_2[ref])

        if self.verbose:
            print('Total clusters to perform local assemblies:', len(self.cluster_to_dir), flush=True)


    def _set_insert_size_data(self):
        assert len(self.insert_hist) > 0
        (x, self.insert_size, pc95, self.insert_sspace_sd) = self.insert_hist.stats()
        self.insert_sspace_sd = min(1, self.insert_sspace_sd)
        self.insert_proper_pair_max = 1.1 * pc95
        if self.verbose:
            print('\nInsert size information from reads mapped to reference genes:')
            print('Insert size:', self.insert_size, sep='\t')
            print('Insert sspace sd:', self.insert_sspace_sd, sep='\t')
            print('Max insert:', self.insert_proper_pair_max, sep='\t')
            print()


    def _init_and_run_clusters(self):
        if len(self.cluster_to_dir) == 0:
            raise Error('Did not get any reads mapped to genes. Cannot continue')

        counter = 0
        cluster_list = []

        for seq_type in sorted(self.cluster_ids):
            if self.cluster_ids[seq_type] is None:
                continue

            for seq_name in sorted(self.cluster_ids[seq_type]):
                if seq_name not in self.cluster_to_dir:
                    continue
                counter += 1
                if self.verbose:
                    print('Constructing cluster', seq_name + '.', counter, 'of', str(len(self.cluster_to_dir)))
                new_dir = self.cluster_to_dir[seq_name]
                self.refdata.write_seqs_to_fasta(os.path.join(new_dir, 'references.fa'), self.cluster_ids[seq_type][seq_name])

                cluster_list.append(cluster.Cluster(
                    new_dir,
                    seq_name,
                    refdata=self.refdata,
                    assembly_kmer=self.assembly_kmer,
                    assembler=self.assembler,
                    max_insert=self.insert_proper_pair_max,
                    min_scaff_depth=self.min_scaff_depth,
                    nucmer_min_id=self.nucmer_min_id,
                    nucmer_min_len=self.nucmer_min_len,
                    nucmer_breaklen=self.nucmer_breaklen,
                    reads_insert=self.insert_size,
                    sspace_k=self.min_scaff_depth,
                    sspace_sd=self.insert_sspace_sd,
                    threads=1, # clusters now run in parallel, so this should always be 1!
                    bcf_min_dp=10,            # let the user change this in a future version?
                    bcf_min_dv=5,             # let the user change this in a future version?
                    bcf_min_dv_over_dp=0.3,   # let the user change this in a future version?
                    bcf_min_qual=20,          # let the user change this in a future version?
                    assembled_threshold=self.assembled_threshold,
                    unique_threshold=self.unique_threshold,
                    bowtie2_preset=self.bowtie2_preset,
                    spades_other_options=self.spades_other,
                    clean=self.clean,
                    extern_progs=self.extern_progs,
                ))


        pool = multiprocessing.Pool(self.threads)
        cluster_list = pool.starmap(_run_cluster, zip(cluster_list, itertools.repeat(self.verbose)))
        self.clusters = {c.name: c for c in cluster_list}


    @staticmethod
    def _write_reports(clusters_in, tsv_out, xls_out):
        columns = copy.copy(report.columns)
        columns[0] = '#' + columns[0]

        f = pyfastaq.utils.open_file_write(tsv_out)
        print('\t'.join(columns), file=f)

        columns[0] = columns[0][1:]
        workbook = openpyxl.Workbook()
        worksheet = workbook.worksheets[0]
        worksheet.title = 'ARIBA_report'
        worksheet.append(columns)

        for seq_name in sorted(clusters_in):
            if clusters_in[seq_name].report_lines is None:
                continue

            for line in clusters_in[seq_name].report_lines:
                print(line, file=f)
                worksheet.append(line.split('\t'))

        pyfastaq.utils.close(f)
        workbook.save(xls_out)


    def _write_catted_assembled_seqs_fasta(self, outfile):
        f = pyfastaq.utils.open_file_write(outfile)

        for gene in sorted(self.clusters):
            try:
                cluster_fasta = self.clusters[gene].assembly_compare.assembled_ref_seqs_file
            except:
                continue
            if os.path.exists(cluster_fasta):
                file_reader = pyfastaq.sequences.file_reader(cluster_fasta)
                for seq in file_reader:
                    print(seq, file=f)

        pyfastaq.utils.close(f)


    def _clean(self):
        if self.clean == 0:
            if self.verbose:
                print('   ... not deleting anything because --clean 0 used')
            return

        to_delete= [
            self.bam,
            self.cdhit_cluster_representatives_fa,
            self.cdhit_cluster_representatives_fa + '.fai',
            self.cdhit_files_prefix + '.non_coding.cdhit',
            self.cdhit_files_prefix + '.presence_absence.cdhit',
            self.cdhit_files_prefix + '.variants_only.cdhit',
        ]

        if self.clean == 2:
            if self.verbose:
                print('    rm -r', self.clusters_outdir)
                shutil.rmtree(self.clusters_outdir)

            to_delete.extend([
                self.cdhit_files_prefix + '.clusters.tsv',
                self.refdata_files_prefix + '.00.check_fasta_presence_absence.log',
                self.refdata_files_prefix + '.00.check_fasta_variants_only.log',
                self.refdata_files_prefix + '.01.check_variants.log',
                self.refdata_files_prefix + '.01.check_variants.non_coding.fa',
                self.refdata_files_prefix + '.01.check_variants.presence_absence.fa',
                self.refdata_files_prefix + '.01.check_variants.tsv',
                self.refdata_files_prefix + '.01.check_variants.variants_only.fa',
            ])

        for filename in to_delete:
            if os.path.exists(filename):
                if self.verbose:
                    print('    rm', filename)
                try:
                    os.unlink(filename)
                except:
                    raise Error('Error deleting file', filename)


    def write_versions_file(self, original_dir):
        with open('version_info.txt', 'w') as f:
            print('ARIBA version', common.version, 'run with this command:', file=f)
            print(' '.join(sys.argv), file=f)
            print('from this directory:', original_dir, file=f)
            print(file=f)
            print('Versions of dependencies:\n', file=f)
            print(*self.extern_progs.version_report, sep='\n', file=f)


    def run(self):
        cwd = os.getcwd()
        os.chdir(self.outdir)
        self.write_versions_file(cwd)

        if self.verbose:
            print('{:_^79}'.format(' Checking reference data '), flush=True)
        self.refdata.sanity_check(self.refdata_files_prefix)

        if self.verbose:
            print()
            print('{:_^79}'.format(' Running cd-hit '), flush=True)
        self._run_cdhit()

        if self.verbose:
            print('Finished cd-hit\n')
            print('{:_^79}'.format(' Mapping reads to clustered genes '), flush=True)
        self._map_reads_to_clustered_genes()

        if self.verbose:
            print('Finished mapping\n')
            print('{:_^79}'.format(' Generating clusters '), flush=True)
        self._bam_to_clusters_reads()

        if len(self.cluster_to_dir) > 0:
            self._set_insert_size_data()
            if self.verbose:
                print('{:_^79}'.format(' Assembling each cluster '))
                print('Will run', self.threads, 'cluster(s) in parallel', flush=True)
            self._init_and_run_clusters()
            if self.verbose:
                print('Finished assembling clusters\n')
        else:
            if self.verbose:
                print('No reads mapped. Skipping all assemblies', flush=True)
            print('WARNING: no reads mapped to reference genes. Therefore no local assemblies will be run', file=sys.stderr)

        if self.verbose:
            print('{:_^79}'.format(' Writing report files '), flush=True)
            print(self.report_file_all_tsv)
            print(self.report_file_all_xls)
        self._write_reports(self.clusters, self.report_file_all_tsv, self.report_file_all_xls)

        if self.verbose:
            print(self.report_file_filtered_prefix + '.tsv')
            print(self.report_file_filtered_prefix + '.xls')
        rf = report_filter.ReportFilter(infile=self.report_file_all_tsv)
        rf.run(self.report_file_filtered_prefix)

        if self.verbose:
            print('{:_^79}'.format(' Writing fasta of assembled sequences '), flush=True)
            print(self.catted_assembled_seqs_fasta)
        self._write_catted_assembled_seqs_fasta(self.catted_assembled_seqs_fasta)

        if self.verbose:
            print('\n\nCleaning files:', flush=True)
        self._clean()

        if self.verbose:
            print('\nAll done!\n')

        os.chdir(cwd)
