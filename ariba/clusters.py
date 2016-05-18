import signal
import time
import atexit
import os
import copy
import tempfile
import pickle
import itertools
import sys
import shutil
import openpyxl
import multiprocessing
import pysam
import pyfastaq
from ariba import cluster, common, mapping, histogram, read_store, report, report_filter, reference_data

class Error (Exception): pass


def _run_cluster(obj, verbose, clean, fails_dir):
    failed_clusters = os.listdir(fails_dir)

    if len(failed_clusters) > 0:
        print('Other clusters failed. Will not start cluster', obj.name, file=sys.stderr)
        return obj

    if verbose:
        print('Start running cluster', obj.name, 'in directory', obj.root_dir, flush=True)
    try:
        obj.run()
    except:
        print('Failed cluster:', obj.name, file=sys.stderr)
        with open(os.path.join(fails_dir, obj.name), 'w'):
            pass

    if verbose:
        print('Finished running cluster', obj.name, 'in directory', obj.root_dir, flush=True)

    if clean:
        if verbose:
            print('Deleting cluster dir', obj.root_dir, flush=True)
        if os.path.exists(obj.root_dir):
            shutil.rmtree(obj.root_dir)

    return obj


class Clusters:
    def __init__(self,
      refdata_dir,
      reads_1,
      reads_2,
      outdir,
      extern_progs,
      version_report_lines=None,
      assembly_kmer=21,
      assembly_coverage=100,
      threads=1,
      verbose=False,
      assembler='spades',
      spades_other=None,
      max_insert=1000,
      min_scaff_depth=10,
      nucmer_min_id=90,
      nucmer_min_len=20,
      nucmer_breaklen=200,
      assembled_threshold=0.95,
      unique_threshold=0.03,
      max_gene_nt_extend=30,
      bowtie2_preset='very-sensitive-local',
      clean=True,
      tmp_dir=None,
    ):
        self.refdata_dir = os.path.abspath(refdata_dir)
        self.refdata, self.cluster_ids = self._load_reference_data_from_dir(refdata_dir)
        self.reads_1 = os.path.abspath(reads_1)
        self.reads_2 = os.path.abspath(reads_2)
        self.outdir = os.path.abspath(outdir)
        self.extern_progs = extern_progs

        if version_report_lines is None:
            self.version_report_lines = []
        else:
            self.version_report_lines = version_report_lines

        self.clean = clean
        self.logs_dir = os.path.join(self.outdir, 'Logs')

        self.assembler = assembler
        assert self.assembler in ['spades']
        self.assembly_kmer = assembly_kmer
        self.assembly_coverage = assembly_coverage
        self.spades_other = spades_other

        self.cdhit_files_prefix = os.path.join(self.refdata_dir, 'cdhit')
        self.cdhit_cluster_representatives_fa = self.cdhit_files_prefix + '.cluster_representatives.fa'
        self.bam_prefix = os.path.join(self.outdir, 'map_reads_to_cluster_reps')
        self.bam = self.bam_prefix + '.bam'
        self.report_file_all_tsv = os.path.join(self.outdir, 'report.all.tsv')
        self.report_file_all_xls = os.path.join(self.outdir, 'report.all.xls')
        self.report_file_filtered_prefix = os.path.join(self.outdir, 'report')
        self.catted_assembled_seqs_fasta = os.path.join(self.outdir, 'assembled_seqs.fa.gz')
        self.catted_genes_matching_refs_fasta = os.path.join(self.outdir, 'assembled_genes.fa.gz')
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
        self.max_gene_nt_extend = max_gene_nt_extend

        self.cluster_to_dir = {}  # gene name -> abs path of cluster directory
        self.clusters = {}        # gene name -> Cluster object
        self.cluster_read_counts = {} # gene name -> number of reads
        self.cluster_base_counts = {} # gene name -> number of bases
        self.pool = None
        self.fails_dir = os.path.join(self.outdir ,'.fails')
        self.clusters_all_ran_ok = True

        for d in [self.outdir, self.logs_dir, self.fails_dir]:
            try:
                os.mkdir(d)
            except:
                raise Error('Error mkdir ' + d)

        if tmp_dir is None:
            if 'ARIBA_TMPDIR' in os.environ:
                tmp_dir = os.path.abspath(os.environ['ARIBA_TMPDIR'])
            else:
                tmp_dir = self.outdir

        if not os.path.exists(tmp_dir):
            raise Error('Temporary directory ' + tmp_dir + ' not found. Cannot continue')

        if self.clean:
            self.tmp_dir_obj = tempfile.TemporaryDirectory(prefix='ariba.tmp.', dir=os.path.abspath(tmp_dir))
            self.tmp_dir = self.tmp_dir_obj.name
        else:
            self.tmp_dir_obj = None
            self.tmp_dir = os.path.join(self.outdir, 'clusters')
            try:
                os.mkdir(self.tmp_dir)
            except:
                raise Error('Error making directory ' + self.tmp_dir)

        if self.verbose:
            print('Temporary directory:', self.tmp_dir)

        for i in [x for x in dir(signal) if x.startswith("SIG") and x not in {'SIGCHLD', 'SIGCLD'}]:
            try:
                signum = getattr(signal, i)
                signal.signal(signum, self._receive_signal)
            except:
                pass


    def _stop_pool(self):
        if self.pool is None:
            return
        self.pool.close()
        self.pool.terminate()
        while len(multiprocessing.active_children()) > 0:
            time.sleep(1)


    def _emergency_stop(self):
        self._stop_pool()
        if self.clean:
            try:
                self.tmp_dir_obj.cleanup()
            except:
                pass


    def _receive_signal(self, signum, stack):
        print('Stopping! Signal received:', signum, file=sys.stderr, flush=True)
        self._emergency_stop()
        sys.exit(1)


    @classmethod
    def _load_reference_data_info_file(cls, filename):
        data = {
            'genetic_code': None
        }

        with open(filename) as f:
            for line in f:
                key, val = line.rstrip().split('\t')
                if key in data:
                    data[key] = val

        if None in data.values():
            missing_values = [x for x in data if data[x] is None]
            raise Error('Error reading reference info file ' + filename + '. These values not found: ' + ','.join(missing_values))

        data['genetic_code'] = int(data['genetic_code'])
        return data


    @staticmethod
    def _load_reference_data_from_dir(indir):
        if not os.path.exists(indir):
            raise Error('Error loading reference data. Input directory ' + indir + ' not found. Cannot continue')

        variants_only_fa = os.path.join(indir, 'refcheck.01.check_variants.variants_only.fa')
        presence_absence_fa = os.path.join(indir, 'refcheck.01.check_variants.presence_absence.fa')
        non_coding_fa = os.path.join(indir, 'refcheck.01.check_variants.non_coding.fa')
        metadata_tsv = os.path.join(indir, 'refcheck.01.check_variants.tsv')
        info_file = os.path.join(indir, 'info.txt')
        clusters_file = os.path.join(indir, 'cdhit.clusters.pickle')
        params = Clusters._load_reference_data_info_file(info_file)
        refdata = reference_data.ReferenceData(
            presence_absence_fa=presence_absence_fa if os.path.exists(presence_absence_fa) else None,
            variants_only_fa=variants_only_fa if os.path.exists(variants_only_fa) else None,
            non_coding_fa=non_coding_fa if os.path.exists(non_coding_fa) else None,
            metadata_tsv=metadata_tsv if os.path.exists(metadata_tsv) else None,
            genetic_code=params['genetic_code'],
        )

        with open(clusters_file, 'rb') as f:
            cluster_ids = pickle.load(f)

        return refdata, cluster_ids


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


    def _bam_to_clusters_reads(self):
        '''Sets up ReadStore of reads for all the clusters. Also gathers histogram data of insert size'''
        reads_file_for_read_store = os.path.join(self.outdir, 'reads')
        f_out = pyfastaq.utils.open_file_write(reads_file_for_read_store)

        sam_reader = pysam.Samfile(self.bam, "rb")
        sam1 = None
        self.proper_pairs = 0

        for s in sam_reader.fetch(until_eof=True):
            if sam1 is None:
                sam1 = s
                continue

            ref_seqs = set()
            if not s.is_unmapped:
                ref_seqs.add(sam_reader.getrname(s.tid))
            if not sam1.is_unmapped:
                ref_seqs.add(sam_reader.getrname(sam1.tid))

            read1 = mapping.sam_to_fastq(sam1)
            read2 = mapping.sam_to_fastq(s)
            if read1.id.endswith('/2'):
                read1, read2 = read2, read1

            insert = mapping.sam_pair_to_insert(s, sam1)
            if insert is not None:
                self.insert_hist.add(insert)
                self.proper_pairs += 1

            for ref in ref_seqs:
                if ref not in self.cluster_to_dir:
                    new_dir = os.path.join(self.tmp_dir, ref)
                    self.cluster_to_dir[ref] = new_dir
                    if self.verbose:
                        print('New cluster with reads that hit:', ref, flush=True)

                self.cluster_read_counts[ref] = self.cluster_read_counts.get(ref, 0) + 2
                self.cluster_base_counts[ref] = self.cluster_base_counts.get(ref, 0) + len(read1) + len(read2)
                print(ref, self.cluster_read_counts[ref] - 1, read1.seq, read1.qual, sep='\t', file=f_out)
                print(ref, self.cluster_read_counts[ref], read2.seq, read2.qual, sep='\t', file=f_out)

            sam1 = None

        pyfastaq.utils.close(f_out)

        if len(self.cluster_read_counts):
            if self.verbose:
                filehandle = sys.stdout
            else:
                filehandle = None

            self.read_store = read_store.ReadStore(
              reads_file_for_read_store,
              os.path.join(self.outdir, 'read_store'),
              log_fh=filehandle
            )

        os.unlink(reads_file_for_read_store)

        if self.verbose:
            print('Found', self.proper_pairs, 'proper read pairs')
            print('Total clusters to perform local assemblies:', len(self.cluster_to_dir), flush=True)


    def _set_insert_size_data(self):
        if len(self.insert_hist) == 0:
            return False
        else:
            (x, self.insert_size, pc95, self.insert_sspace_sd) = self.insert_hist.stats()
            self.insert_sspace_sd = min(1, self.insert_sspace_sd)
            self.insert_proper_pair_max = 1.1 * pc95
            if self.verbose:
                print('\nInsert size information from reads mapped to reference genes:')
                print('Insert size:', self.insert_size, sep='\t')
                print('Insert sspace sd:', self.insert_sspace_sd, sep='\t')
                print('Max insert:', self.insert_proper_pair_max, sep='\t')
                print()
            return True


    def _init_and_run_clusters(self):
        if len(self.cluster_to_dir) == 0:
            raise Error('Did not get any reads mapped to genes. Cannot continue')

        counter = 0
        cluster_list = []
        self.log_files = []

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
                self.log_files.append(os.path.join(self.logs_dir, seq_name + '.log'))

                cluster_list.append(cluster.Cluster(
                    new_dir,
                    seq_name,
                    self.refdata,
                    self.cluster_read_counts[seq_name],
                    self.cluster_base_counts[seq_name],
                    fail_file=os.path.join(self.fails_dir, seq_name),
                    read_store=self.read_store,
                    reference_names=self.cluster_ids[seq_type][seq_name],
                    logfile=self.log_files[-1],
                    assembly_coverage=self.assembly_coverage,
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
                    max_gene_nt_extend=self.max_gene_nt_extend,
                    bowtie2_preset=self.bowtie2_preset,
                    spades_other_options=self.spades_other,
                    clean=self.clean,
                    extern_progs=self.extern_progs,
                ))

        try:
            if self.threads > 1:
                self.pool = multiprocessing.Pool(self.threads)
                cluster_list = self.pool.starmap(_run_cluster, zip(cluster_list, itertools.repeat(self.verbose), itertools.repeat(self.clean), itertools.repeat(self.fails_dir)))
            else:
                for c in cluster_list:
                    _run_cluster(c, self.verbose, self.clean, self.fails_dir)
        except:
            self.clusters_all_ran_ok = False

        if len(os.listdir(self.fails_dir)) > 0:
            self.clusters_all_ran_ok = False

        self.clusters = {c.name: c for c in cluster_list}


    @staticmethod
    def _write_reports(clusters_in, tsv_out, xls_out=None):
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
        if xls_out is not None:
            workbook.save(xls_out)


    def _write_catted_assembled_seqs_fasta(self, outfile):
        f = pyfastaq.utils.open_file_write(outfile)

        for gene in sorted(self.clusters):
            try:
                seq_dict = self.clusters[gene].assembly_compare.assembled_reference_sequences
            except:
                continue

            for seq_name in sorted(seq_dict):
                print(seq_dict[seq_name], file=f)

        pyfastaq.utils.close(f)


    def _write_catted_genes_matching_refs_fasta(self, outfile):
        f = pyfastaq.utils.open_file_write(outfile)

        for gene in sorted(self.clusters):
            if self.clusters[gene].assembly_compare is not None and self.clusters[gene].assembly_compare.gene_matching_ref is not None:
                seq = copy.copy(self.clusters[gene].assembly_compare.gene_matching_ref)
                seq.id += '.' + '.'.join([
                    self.clusters[gene].assembly_compare.gene_matching_ref_type,
                    str(self.clusters[gene].assembly_compare.gene_start_bases_added),
                    str(self.clusters[gene].assembly_compare.gene_end_bases_added)
                ])
                print(seq, file=f)

        pyfastaq.utils.close(f)


    def _clean(self):
        if self.clean:
            shutil.rmtree(self.fails_dir)

            try:
                self.tmp_dir_obj.cleanup()
            except:
                pass

            if self.verbose:
                print('Deleting Logs directory', self.logs_dir)
            try:
                shutil.rmtree(self.logs_dir)
            except:
                pass

            if self.verbose:
                print('Deleting reads store files', self.read_store.outfile + '[.tbi]')
            try:
                self.read_store.clean()
            except:
                pass
        else:
            if self.verbose:
                print('Not deleting anything because --noclean used')


    def write_versions_file(self, original_dir):
        with open('version_info.txt', 'w') as f:
            print('ARIBA run with this command:', file=f)
            print(' '.join([sys.argv[0]] + ['run'] + sys.argv[1:]), file=f)
            print('from this directory:', original_dir, file=f)
            print(file=f)
            print(*self.version_report_lines, sep='\n', file=f)


    def run(self):
        try:
            self._run()
        except Error as err:
            self._emergency_stop()
            raise Error('Something went wrong during ariba run. Cannot continue. Error was:\n' + str(err))


    def _run(self):
        cwd = os.getcwd()
        os.chdir(self.outdir)
        self.write_versions_file(cwd)

        if self.verbose:
            print('{:_^79}'.format(' Mapping reads to clustered genes '), flush=True)
        self._map_reads_to_clustered_genes()

        if self.verbose:
            print('Finished mapping\n')
            print('{:_^79}'.format(' Generating clusters '), flush=True)
        self._bam_to_clusters_reads()
        if self.clean:
            if self.verbose:
                print('Deleting BAM', self.bam, flush=True)
            os.unlink(self.bam)

        if len(self.cluster_to_dir) > 0:
            got_insert_data_ok = self._set_insert_size_data()
            if not got_insert_data_ok:
                print('WARNING: not enough proper read pairs (found ' + str(self.proper_pairs) + ') to determine insert size.', file=sys.stderr)
                print('This probably means that very few reads were mapped at all. No local assemblies will be run', file=sys.stderr)
                if self.verbose:
                    print('Not enough proper read pairs mapped to determine insert size. Skipping all assemblies.', flush=True)
            else:
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

        if not self.clusters_all_ran_ok:
            raise Error('At least one cluster failed! Stopping...')

        if self.verbose:
            print('{:_^79}'.format(' Writing reports '), flush=True)
            print('Making', self.report_file_all_tsv)
        self._write_reports(self.clusters, self.report_file_all_tsv)

        if self.verbose:
            print('Making', self.report_file_filtered_prefix + '.tsv')
        rf = report_filter.ReportFilter(infile=self.report_file_all_tsv)
        rf.run(self.report_file_filtered_prefix)

        if self.verbose:
            print()
            print('{:_^79}'.format(' Writing fasta of assembled sequences '), flush=True)
            print(self.catted_assembled_seqs_fasta, 'and', self.catted_genes_matching_refs_fasta, flush=True)
        self._write_catted_assembled_seqs_fasta(self.catted_assembled_seqs_fasta)
        self._write_catted_genes_matching_refs_fasta(self.catted_genes_matching_refs_fasta)

        clusters_log_file = os.path.join(self.outdir, 'log.clusters.gz')
        if self.verbose:
            print()
            print('{:_^79}'.format(' Catting cluster log files '), flush=True)
            print('Writing file', clusters_log_file, flush=True)
        common.cat_files(self.log_files, clusters_log_file)

        if self.verbose:
            print()
            print('{:_^79}'.format(' Cleaning files '), flush=True)
        self._clean()

        if self.clusters_all_ran_ok and self.verbose:
            print('\nAll done!\n')

        os.chdir(cwd)
