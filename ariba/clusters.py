import os
import sys
import shutil
import openpyxl
import pysam
import pyfastaq
from ariba import cluster, common, mapping, histogram

class Error (Exception): pass


class Clusters:
    def __init__(self,
      db_fasta,
      reads_1,
      reads_2,
      outdir,
      assembly_kmer=21,
      threads=1,
      verbose=False,
      assembler='velvet',
      smalt_k=13,
      smalt_s=2,
      smalt_min_id=0.9,
      spades_other=None,
      max_insert=1000,
      min_scaff_depth=10,
      nucmer_min_id=90,
      nucmer_min_len=50,
      nucmer_breaklen=50,
      assembled_threshold=0.95,
      unique_threshold=0.03,
      bcftools_exe='bcftools',
      gapfiller_exe='GapFiller.pl',
      samtools_exe='samtools',
      smalt_exe='smalt',
      spades_exe='spades.py',
      sspace_exe='SSPACE_Basic_v2.0.pl',
      velvet_exe='velvet', # prefix of velvet{g,h}
    ):
        self.db_fasta = os.path.abspath(db_fasta)
        self.reads_1 = os.path.abspath(reads_1)
        self.reads_2 = os.path.abspath(reads_2)
        self.outdir = os.path.abspath(outdir)

        self.assembler = assembler
        assert self.assembler in ['velvet', 'spades']
        self.assembly_kmer = assembly_kmer
        self.spades_other = spades_other

        self.bam_prefix = os.path.join(self.outdir, 'map_all_reads')
        self.bam = self.bam_prefix + '.bam'
        self.report_file_tsv = os.path.join(self.outdir, 'report.tsv')
        self.report_file_xls = os.path.join(self.outdir, 'report.xls')
        self.threads = threads
        self.verbose = verbose

        self.smalt_k = smalt_k
        self.smalt_s = smalt_s
        self.smalt_min_id = smalt_min_id
        self.max_insert = max_insert
        self.smalt_exe = smalt_exe

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

        self.bcftools_exe = bcftools_exe

        self.gapfiller_exe = shutil.which(gapfiller_exe)
        if self.gapfiller_exe is None:
            raise Error('Error! ' + gapfiller_exe + ' not found in path')
        self.gapfiller_exe = os.path.realpath(self.gapfiller_exe) # otherwise gapfiller dies loading packages

        self.samtools_exe = samtools_exe
        self.spades_exe = spades_exe

        self.sspace_exe = shutil.which(sspace_exe)
        if self.sspace_exe is None:
            raise Error('Error! ' + sspace_exe + ' not found in path')
        self.sspace_exe = os.path.realpath(self.sspace_exe) # otherwise sspace dies loading packages

        self.velvet = velvet_exe

        try:
            os.mkdir(self.outdir)
        except:
            raise Error('Error mkdir ' + self.outdir)


    def _map_reads(self):
        mapping.run_smalt(
            self.reads_1,
            self.reads_2,
            self.db_fasta,
            self.bam_prefix,
            index_k=self.smalt_k,
            index_s=self.smalt_s,
            threads=self.threads,
            samtools=self.samtools_exe,
            smalt=self.smalt_exe,
            minid=self.smalt_min_id,
            verbose=self.verbose,
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

                    new_dir = os.path.join(self.outdir, ref)
                    try:
                        os.mkdir(new_dir)
                    except:
                        raise Error('Error mkdir ' + new_dir)

                    self.cluster_to_dir[ref] = new_dir
                    filehandles_1[ref] = pyfastaq.utils.open_file_write(os.path.join(new_dir, 'reads_1.fq'))
                    filehandles_2[ref] = pyfastaq.utils.open_file_write(os.path.join(new_dir, 'reads_2.fq'))

                print(read1, file=filehandles_1[ref])
                print(read2, file=filehandles_2[ref])

            sam1 = None

        for ref in filehandles_1:
            pyfastaq.utils.close(filehandles_1[ref])
            pyfastaq.utils.close(filehandles_2[ref])


    def _set_insert_size_data(self):
        assert len(self.insert_hist) > 0
        (x, self.insert_size, pc95, self.insert_sspace_sd) = self.insert_hist.stats()
        self.insert_proper_pair_max = 1.1 * pc95
        if self.verbose:
            print('\nInsert size information from reads mapped to reference genes:')
            print('Insert size:', self.insert_size, sep='\t')
            print('Insert sspace sd:', self.insert_sspace_sd, sep='\t')
            print('Max insert:', self.insert_proper_pair_max, sep='\t')
            print()


    def _write_gene_fa(self, gene_name, outfile):
        if not os.path.exists(self.db_fasta + '.fai'):
            common.syscall(self.samtools_exe + ' faidx ' + self.db_fasta, verbose=self.verbose)

        common.syscall(' '.join([
            self.samtools_exe + ' faidx',
            self.db_fasta,
            gene_name,
            '>', outfile
        ]))



    def _init_and_run_clusters(self):
        if len(self.cluster_to_dir) == 0:
            raise Error('Did not get any reads mapped to genes. Cannot continue')

        counter = 0

        for gene in sorted(self.cluster_to_dir):
            counter += 1
            if self.verbose:
                print('\nAssembling cluster', counter, 'of', str(len(self.cluster_to_dir)) + ':', gene)
            new_dir = self.cluster_to_dir[gene]
            self._write_gene_fa(gene, os.path.join(new_dir, 'gene.fa'))
            self.clusters[gene] = cluster.Cluster(
                new_dir,
                assembly_kmer=self.assembly_kmer,
                assembler=self.assembler,
                max_insert=self.insert_proper_pair_max,
                min_scaff_depth=self.min_scaff_depth,
                nucmer_min_id=self.nucmer_min_id,
                nucmer_min_len=self.nucmer_min_len,
                nucmer_breaklen=self.nucmer_breaklen,
                sspace_k=self.min_scaff_depth,
                reads_insert=self.insert_size,
                sspace_sd=self.insert_sspace_sd,
                threads=self.threads,
                assembled_threshold=self.assembled_threshold,
                unique_threshold=self.unique_threshold,
                verbose=self.verbose,
                bcftools_exe=self.bcftools_exe,
                gapfiller_exe=self.gapfiller_exe,
                samtools_exe=self.samtools_exe,
                spades_exe=self.spades_exe,
                sspace_exe=self.sspace_exe,
                velvet_exe=self.velvet,
                spades_other=self.spades_other
            )

            self.clusters[gene].run()


    def _write_reports(self):
        columns = [
            '#gene',
            'flag',
            'gene_len',
            'var_type',
            'var_effect',
            'new_aa',
            'gene_start',
            'gene_end',
            'gene_nt',
            'scaffold',
            'scaff_len',
            'scaff_start',
            'scaff_end',
            'scaff_nt',
        ]

        f = pyfastaq.utils.open_file_write(self.report_file_tsv)
        print('\t'.join(columns), file=f)

        columns[0] = 'gene'
        workbook = openpyxl.Workbook()
        worksheet = workbook.worksheets[0] 
        worksheet.title = 'ARIBA_report'
        worksheet.append(columns)

        for gene in sorted(self.clusters):
            for line in self.clusters[gene].report_lines:
                print('\t'.join([str(x) for x in line]), file=f)
                worksheet.append(line)
        pyfastaq.utils.close(f)
        workbook.save(self.report_file_xls)




    def run(self):
        if self.verbose:
            print('{:_^79}'.format(' Mapping reads to reference genes '))
        self._map_reads()
        if self.verbose:
            print('Finished mapping\n')
            print('{:_^79}'.format(' Generating clusters '))
        self._bam_to_clusters_reads()
        self._set_insert_size_data()
        if self.verbose:
            print('{:_^79}'.format(' Assembling each cluster '))
        self._init_and_run_clusters()
        if self.verbose:
            print('Finished assembling clusters\n')
            print('{:_^79}'.format(' Writing report files '))
        self._write_reports()
        if self.verbose:
            print('Finished writing report files. All done!')
