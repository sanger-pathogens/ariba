import shutil
import sys
import os
import pyfastaq
from ariba import mlst_profile, pubmlst_getter, ref_preparer, versions

class Error (Exception): pass


class PubmlstRefPreparer:
    def __init__(self, species, outdir, debug=False, verbose=False):
        self.species = species
        self.outdir = outdir
        self.debug = debug
        self.verbose = verbose
        self.clusters_file = os.path.join(outdir, 'clusters.tsv')
        self.mlst_download_dir = os.path.join(self.outdir, 'pubmlst_download')
        self.prepareref_dir = os.path.join(self.outdir, 'ref_db')
        self.extern_progs, version_report_lines = versions.get_all_versions()


    @classmethod
    def _filter_seq_dict(cls, d):
        lengths = [len(x) for x in d.values()]
        lengths.sort()
        median_length = lengths[int(len(d) / 2)]
        min_length = 0.9 * median_length
        max_length = 1.1 * median_length
        to_delete = []
        for seq in d:
            if not min_length <= len(d[seq]) <= max_length:
                to_delete.append(seq)
                print('WARNING: Median sequence length is', median_length, 'but', seq, 'has length', len(d[seq]), 'which is too long or short. Removing.', file=sys.stderr)

        for seq in to_delete:
            del d[seq]


    def _load_fasta_files_and_write_clusters_file(self, indir):
        self.sequences = {}
        self.fasta_files = []
        clusters_fh = pyfastaq.utils.open_file_write(self.clusters_file)

        for gene_name in self.profile.genes_list:
            infile = os.path.join(indir, gene_name + '.tfa')

            if not os.path.exists(infile):
                pyfastaq.utils.close(clusters_fh)
                raise Error('Cannot find file "' + infile + '" for gene ' + gene_name)

            self.sequences[gene_name] = {}
            pyfastaq.tasks.file_to_dict(infile, self.sequences[gene_name])
            PubmlstRefPreparer._filter_seq_dict(self.sequences[gene_name])
            seq_names = sorted(list(self.sequences[gene_name].keys()))
            print(*seq_names, sep='\t', file=clusters_fh)

            if self.verbose:
                print('Loaded fasta file for gene', gene_name)

            with open(infile, 'w') as f:
                for seq in sorted(self.sequences[gene_name]):
                    print(self.sequences[gene_name][seq], file=f)

            self.fasta_files.append(infile)

        pyfastaq.utils.close(clusters_fh)


    def run(self):
        try:
            os.mkdir(self.outdir)
        except:
            raise Error('Error making output directory ' + self.outdir)

        pubmlst = pubmlst_getter.PubmlstGetter(debug=self.debug, verbose=self.verbose)
        pubmlst.get_species_files(self.species, self.mlst_download_dir)
        if self.verbose:
            print('Downloaded data from pubmlst')

        profile_file = os.path.join(self.mlst_download_dir, 'profile.txt')
        self.profile = mlst_profile.MlstProfile(profile_file, duplicate_warnings=True)
        if self.verbose:
            print('Loaded mlst profile file', profile_file)

        self._load_fasta_files_and_write_clusters_file(self.mlst_download_dir)
        if self.verbose:
            print('Loaded fasta files and wrote clusters file')
            print('Putting data in ariba db directory', self.prepareref_dir)

        refprep = ref_preparer.RefPreparer(
            self.fasta_files,
            self.extern_progs,
            all_coding='no',
            clusters_file=self.clusters_file,
            verbose=self.verbose,
        )
        refprep.run(self.prepareref_dir)
        shutil.copy(profile_file, os.path.join(self.prepareref_dir, 'pubmlst.profile.txt'))

        print('ariba db directory prepared. You can use it like this:')
        print('ariba run', self.prepareref_dir, 'reads_1.fq reads_2.fq output_directory')
