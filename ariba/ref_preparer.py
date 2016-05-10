import sys
import os
import pickle
from ariba import common, mapping, reference_data

class Error (Exception): pass


class RefPreparer:
    def __init__(self,
        extern_progs,
        version_report_lines=None,
        ref_prefix=None,
        presabs=None,
        varonly=None,
        noncoding=None,
        metadata=None,
        min_gene_length=6,
        max_gene_length=10000,
        genetic_code=11,
        cdhit_min_id=0.9,
        cdhit_min_length=0.9,
        run_cdhit=True,
        clusters_file=None,
        threads=1,
        verbose=False,
    ):
        self.extern_progs = extern_progs

        if version_report_lines is None:
            self.version_report_lines = []
        else:
            self.version_report_lines = version_report_lines

        self.ref_prefix = ref_prefix
        self.presabs = presabs
        self.varonly = varonly
        self.noncoding = noncoding
        self.metadata = metadata
        self.min_gene_length = min_gene_length
        self.max_gene_length = max_gene_length
        self.genetic_code = genetic_code
        self.cdhit_min_id = cdhit_min_id
        self.cdhit_min_length = cdhit_min_length
        self.run_cdhit = run_cdhit
        self.clusters_file = clusters_file
        self.threads = threads
        self.verbose = verbose


    @staticmethod
    def _get_ref_files(ref_prefix, presabs, varonly, noncoding, metadata, verbose=False):
        if {None} == {ref_prefix, presabs, varonly, noncoding}:
            raise Error('Error in RefPreparer._get_ref_files. All input files and ref_prefix were None. Cannot continue')

        filenames = {
            'presabs': presabs,
            'varonly': varonly,
            'noncoding': noncoding,
            'metadata': metadata,
        }

        file_suffixes = {
            'presabs': 'presence_absence.fa',
            'varonly': 'variants_only.fa',
            'noncoding': 'noncoding.fa',
            'metadata': 'metadata.tsv',
        }

        if verbose:
            print('\nLooking for input files ...')

        for key in file_suffixes:
            if ref_prefix is not None:
                filename = os.path.abspath(ref_prefix + '.' + file_suffixes[key])

                if os.path.exists(filename):
                    if verbose:
                        print('Found: ', filename, '.\n    ...treating it as if this was used: --', key, ' ', filename, sep='')
                    filenames[key] = filename
                else:
                    if verbose:
                        print('Not found:', filename)
                    filenames[key] = None
            elif filenames[key] is not None:
                if os.path.exists(filenames[key]):
                    filenames[key] = os.path.abspath(filenames[key])
                    if verbose:
                        print('Found: ', filenames[key], ' from option --', key, sep='')
                else:
                    raise Error('File not found! Cannot continue. Looked for: ' + filenames[key])

        if {None} == {filenames['presabs'], filenames['varonly'], filenames['noncoding']}:
            raise Error('Error in RefPreparer._get_ref_files. No FASTA files given! Cannot continue')

        return filenames


    def _write_info_file(self, outfile):
        with open(outfile, 'w') as fout:
            for key in ('presabs', 'varonly', 'noncoding', 'metadata'):
                print('input_' + key, self.filenames[key], sep='\t', file=fout)

            print('genetic_code', self.genetic_code, sep='\t', file=fout)


    def run(self, outdir):
        original_dir = os.getcwd()

        if os.path.exists(outdir):
            raise Error('Error! Output directory ' + outdir + ' already exists. Cannot continue')

        try:
            os.mkdir(outdir)
        except:
            raise Error('Error making output directory ' + outdir + '. Cannot continue')

        with open(os.path.join(outdir, 'version_info.txt'), 'w') as f:
            print('ARIBA run with this command:', file=f)
            print(' '.join([sys.argv[0]] + ['prepareref'] + sys.argv[1:]), file=f)
            print('from this directory:', original_dir, file=f)
            print(file=f)
            print(*self.version_report_lines, sep='\n', file=f)

        self.filenames = self._get_ref_files(self.ref_prefix, self.presabs, self.varonly, self.noncoding, self.metadata, self.verbose)
        self._write_info_file(os.path.join(outdir, 'info.txt'))

        self.refdata = reference_data.ReferenceData(
            presence_absence_fa=self.filenames['presabs'],
            variants_only_fa=self.filenames['varonly'],
            non_coding_fa=self.filenames['noncoding'],
            metadata_tsv=self.filenames['metadata'],
            min_gene_length=self.min_gene_length,
            max_gene_length=self.max_gene_length,
            genetic_code=self.genetic_code,
        )

        if self.verbose:
            print('\nLoading and checking input data', flush=True)

        refdata_outprefix = os.path.join(outdir, 'refcheck')
        self.refdata.rename_sequences(refdata_outprefix + '.rename_info')
        self.refdata.sanity_check(refdata_outprefix)

        if self.verbose:
            print('\nRunning cdhit', flush=True)
        cdhit_outprefix = os.path.join(outdir, 'cdhit')

        clusters = self.refdata.cluster_with_cdhit(
            refdata_outprefix + '.01.check_variants',
            cdhit_outprefix,
            seq_identity_threshold=self.cdhit_min_id,
            threads=self.threads,
            length_diff_cutoff=self.cdhit_min_length,
            nocluster=not self.run_cdhit,
            verbose=self.verbose,
            cd_hit_est=self.extern_progs.exe('cdhit'),
            clusters_file=self.clusters_file,
        )

        if self.verbose:
            print('\nWriting clusters to file.', len(clusters), 'in total', flush=True)

        clusters_pickle_file = cdhit_outprefix + '.clusters.pickle'
        with open(clusters_pickle_file, 'wb') as f:
            pickle.dump(clusters, f)

        cluster_representatives_fa = cdhit_outprefix + '.cluster_representatives.fa'

        if self.verbose:
            print('\nRunning bowtie2-build on FASTA of cluster representatives')

        mapping.bowtie2_index(
            cluster_representatives_fa,
            cluster_representatives_fa,
            bowtie2=self.extern_progs.exe('bowtie2'),
            verbose=self.verbose,
        )

        if self.verbose:
            print('\nRunning samtools faidx on FASTA of cluster representatives')

        cmd = ' '.join([
            self.extern_progs.exe('samtools'),
            'faidx',
            cluster_representatives_fa
        ])

        common.syscall(cmd, verbose=self.verbose)
