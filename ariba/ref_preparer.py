import sys
import os
import shutil
import pickle
import pyfastaq
from ariba import reference_data, mash

class Error (Exception): pass


class RefPreparer:
    def __init__(self,
        fasta_files,
        extern_progs,
        metadata_tsv_files=None,
        all_coding=None,
        version_report_lines=None,
        min_gene_length=6,
        max_gene_length=10000,
        genetic_code=11,
        cdhit_min_id=0.9,
        cdhit_min_length=0.9,
        run_cdhit=True,
        clusters_file=None,
        threads=1,
        verbose=False,
        force=False,
    ):
        self.extern_progs = extern_progs

        if version_report_lines is None:
            self.version_report_lines = []
        else:
            self.version_report_lines = version_report_lines

        self.fasta_files = fasta_files
        self.metadata_tsv_files = [] if metadata_tsv_files is None else metadata_tsv_files
        self.all_coding = all_coding
        self.min_gene_length = min_gene_length
        self.max_gene_length = max_gene_length
        self.genetic_code = genetic_code
        self.cdhit_min_id = cdhit_min_id
        self.cdhit_min_length = cdhit_min_length
        self.run_cdhit = run_cdhit
        self.clusters_file = clusters_file
        self.threads = threads
        self.verbose = verbose
        self.force = force


    @classmethod
    def _fasta_to_metadata(cls, infile, out_fh, all_coding):
       seq_reader = pyfastaq.sequences.file_reader(infile)
       coding = '1' if all_coding else '0'

       for seq in seq_reader:
           fields = seq.id.split(maxsplit=1)
           if len(fields) > 1:
               info_column = 'Original name: ' + seq.id
               seq.id = fields[0]
           else:
               info_column = '.'
           print(seq.id, coding, 0, '.', '.', info_column, sep='\t', file=out_fh)


    def _write_info_file(self, outfile):
        with open(outfile, 'w') as fout:
            for filename in self.fasta_files:
                print('input fasta file:', filename, sep='\t', file=fout)

            for filename in self.metadata_tsv_files:
                print('input tsv file:', filename, sep='\t', file=fout)

            print('genetic_code', self.genetic_code, sep='\t', file=fout)


    @staticmethod
    def _rename_clusters(clusters_in):
        new_clusters = {}
        key_count = {}
        min_prefix_length = 3

        for old_name, name_set in sorted(clusters_in.items()):
            names_before_dots = {}
            for name in name_set:
                if '.' in name:
                    prefix = name.split('.')[0]
                    if len(prefix) >= min_prefix_length:
                        names_before_dots[prefix] = names_before_dots.get(prefix, 0) + 1

            if len(names_before_dots) == 0:
                new_key = 'cluster'
            elif len(names_before_dots) == 1:
                new_key = list(names_before_dots.keys())[0]
                if sum(list(names_before_dots.values())) < len(name_set):
                    new_key += '+'
            else:
                common_prefix = os.path.commonprefix(list(names_before_dots.keys()))
                if common_prefix == '' or len(common_prefix) < min_prefix_length:
                    max_value = max(list(names_before_dots.values()))
                    possible_keys = [x for x in names_before_dots if names_before_dots[x] == max_value]
                    possible_keys.sort()
                    new_key = possible_keys[0] + '+'
                else:
                    new_key = common_prefix + '-'

            i = 1
            new_new_key = new_key
            while new_new_key in new_clusters:
                new_new_key = new_key + '_' + str(i)
                i += 1
            new_key = new_new_key

            if new_key in key_count:
                if new_key in new_clusters:
                    assert key_count[new_key] == 1

                    i = 1
                    while new_key + '_' + str(i) in new_clusters:
                        i += 1

                    assert new_key + '_' + str(i) not in new_clusters
                    new_clusters[new_key + '_' + str(i)] = new_clusters[new_key]
                    del new_clusters[new_key]
                    key_count[new_key] = i

                key_count[new_key] += 1
                new_name = new_key + '_' + str(key_count[new_key])
            else:
                key_count[new_key] = 1
                new_name = new_key

            assert new_name not in new_clusters
            new_clusters[new_name] = name_set

        return new_clusters


    def run(self, outdir):
        original_dir = os.getcwd()

        if self.force and os.path.exists(outdir):
            shutil.rmtree(outdir)

        if os.path.exists(outdir):
            raise Error('Error! Output directory ' + outdir + ' already exists. Cannot continue')

        try:
            os.mkdir(outdir)
        except:
            raise Error('Error making output directory ' + outdir + '. Cannot continue')

        with open(os.path.join(outdir, '00.version_info.txt'), 'w') as f:
            print('ARIBA run with this command:', file=f)
            print(' '.join([sys.argv[0]] + ['prepareref'] + sys.argv[1:]), file=f)
            print('from this directory:', original_dir, file=f)
            print(file=f)
            print(*self.version_report_lines, sep='\n', file=f)

        if self.all_coding is not None:
            assert len(self.metadata_tsv_files) == 0
            assert self.all_coding in {'yes', 'no'}
            self.metadata_tsv_files = [os.path.join(outdir, '00.auto_metadata.tsv')]
            f_out = pyfastaq.utils.open_file_write(self.metadata_tsv_files[0])
            for fasta_file in self.fasta_files:
                RefPreparer._fasta_to_metadata(fasta_file, f_out, self.all_coding=='yes')
            pyfastaq.utils.close(f_out)
        else:
            assert self.all_coding is None
            assert len(self.metadata_tsv_files) > 0

        self._write_info_file(os.path.join(outdir, '00.info.txt'))

        self.refdata = reference_data.ReferenceData(
            self.fasta_files,
            self.metadata_tsv_files,
            min_gene_length=self.min_gene_length,
            max_gene_length=self.max_gene_length,
            genetic_code=self.genetic_code,
        )

        if self.verbose:
            print('\nLoading and checking input data', flush=True)

        self.refdata.rename_sequences(os.path.join(outdir, '00.rename_info'))
        self.refdata.sanity_check(os.path.join(outdir, '01.filter'))

        if self.verbose:
            print('\nRunning cdhit', flush=True)
        cdhit_outprefix = os.path.join(outdir, '02.cdhit')

        clusters = self.refdata.cluster_with_cdhit(
            cdhit_outprefix,
            seq_identity_threshold=self.cdhit_min_id,
            threads=self.threads,
            length_diff_cutoff=self.cdhit_min_length,
            nocluster=not self.run_cdhit,
            verbose=self.verbose,
            clusters_file=self.clusters_file,
        )

        clusters = self._rename_clusters(clusters)
        reference_data.ReferenceData.write_cluster_allocation_file(clusters, cdhit_outprefix + '.clusters.tsv')

        if self.verbose:
            print('\nWriting clusters to file.', len(clusters), 'in total', flush=True)

        clusters_pickle_file = cdhit_outprefix + '.clusters.pickle'

        with open(clusters_pickle_file, 'wb') as f:
            pickle.dump(clusters, f)

        if self.verbose:
            print('\nMash-sketching all reference sequences', flush=True)

        mash.Masher.sketch(os.path.join(outdir, '02.cdhit.all.fa'), True, self.extern_progs, self.verbose, sys.stdout)

