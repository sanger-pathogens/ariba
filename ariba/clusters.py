import os
import copy
import pickle
import itertools
import sys
import shutil
import openpyxl
import multiprocessing
import pysam
import pyfastaq
from ariba import cluster, common, mapping, histogram, report, report_filter, reference_data

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
      refdata_dir,
      reads_1,
      reads_2,
      outdir,
      extern_progs,
      assembly_kmer=21,
      assembly_coverage=100,
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
      clean=1,
    ):
        self.refdata_dir = os.path.abspath(refdata_dir)
        self.refdata, self.cluster_ids = self._load_reference_data_from_dir(refdata_dir)
        self.reads_1 = os.path.abspath(reads_1)
        self.reads_2 = os.path.abspath(reads_2)
        self.outdir = os.path.abspath(outdir)
        self.extern_progs = extern_progs
        self.clusters_outdir = os.path.join(self.outdir, 'Clusters')
        self.clean = clean

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
        self.cluster_read_counts = {} # gene name -> number of reads
        self.cluster_base_counts = {} # gene name -> number of bases

        for d in [self.outdir, self.clusters_outdir]:
            try:
                os.mkdir(d)
            except:
                raise Error('Error mkdir ' + d)


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

            read1 = mapping.sam_to_fastq(sam1)
            read2 = mapping.sam_to_fastq(s)
            if read1.id.endswith('/2'):
                read1, read2 = read2, read1

            insert = mapping.sam_pair_to_insert(s, sam1)
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
                self.cluster_read_counts[ref] = self.cluster_read_counts.get(ref, 0) + 2
                self.cluster_base_counts[ref] = self.cluster_base_counts.get(ref, 0) + len(read1) + len(read2)

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
                    self.refdata,
                    self.cluster_read_counts[seq_name],
                    self.cluster_base_counts[seq_name],
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
                    bowtie2_preset=self.bowtie2_preset,
                    spades_other_options=self.spades_other,
                    clean=self.clean,
                    extern_progs=self.extern_progs,
                ))


        if self.threads > 1:
            pool = multiprocessing.Pool(self.threads)
            cluster_list = pool.starmap(_run_cluster, zip(cluster_list, itertools.repeat(self.verbose)))
        else:
            for c in cluster_list:
                _run_cluster(c, self.verbose)

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

        to_delete= [self.bam]

        if self.clean == 2:
            if self.verbose:
                print('    rm -r', self.clusters_outdir)
                shutil.rmtree(self.clusters_outdir)

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
            print(' '.join([sys.argv[0]] + ['run'] + sys.argv[1:]), file=f)
            print('from this directory:', original_dir, file=f)
            print(file=f)
            print('Versions of dependencies:\n', file=f)
            print(*self.extern_progs.version_report, sep='\n', file=f)


    def run(self):
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
