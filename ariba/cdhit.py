import tempfile
import shutil
import sys
import os
import pyfastaq
from ariba import common, external_progs

class Error (Exception): pass

class Runner:
    def __init__(
      self,
      infile,
      seq_identity_threshold=0.9,
      threads=1,
      length_diff_cutoff=0.9,
      verbose=False,
      min_cluster_number=0
    ):

        if not os.path.exists(infile):
            raise Error('File not found: "' + infile + '". Cannot continue')

        self.infile = os.path.abspath(infile)
        self.seq_identity_threshold = seq_identity_threshold
        self.threads = threads
        self.length_diff_cutoff = length_diff_cutoff
        self.verbose = verbose
        self.min_cluster_number = min_cluster_number
        extern_progs = external_progs.ExternalProgs(fail_on_error=True)
        self.cd_hit_est = extern_progs.exe('cdhit')


    def fake_run(self):
        '''Doesn't actually run cd-hit. Instead, puts each input sequence into its own cluster. So it's as if cdhit was run, but didn't cluster anything'''
        clusters = {}
        used_names = set()
        seq_reader = pyfastaq.sequences.file_reader(self.infile)

        for seq in seq_reader:
            if seq.id in used_names:
                raise Error('Sequence name "' + seq.id + '" not unique. Cannot continue')

            clusters[str(len(clusters) + self.min_cluster_number)] = {seq.id}
            used_names.add(seq.id)

        return clusters


    @staticmethod
    def _load_user_clusters_file(filename, all_ref_seqs, rename_dict=None):
        if rename_dict is None:
            rename_dict = {}

        f = pyfastaq.utils.open_file_read(filename)
        clusters = {}
        used_names = set()

        for line in f:
            names_list = line.rstrip().split()
            to_remove = set()

            for name in names_list:
                new_name = rename_dict.get(name, name)
                if new_name not in all_ref_seqs:
                    to_remove.add(name)
                    print('WARNING: ignoring sequence', name, 'from clusters file because not in fasta file. This probably means it failed sanity checks - see the log files 01.filter.check_genes.log, 01.filter.check_metadata.log.', file=sys.stderr)

            names_list = [x for x in names_list if x not in to_remove]
            new_names = set([rename_dict.get(name, name) for name in names_list])
            if len(names_list) != len(new_names) or not new_names.isdisjoint(used_names):
                pyfastaq.utils.close(f)
                raise Error('Error in user-provided clusters file ' + filename + '. Non unique name found at this line:\n' + line)

            clusters[str(len(clusters))] = new_names
            used_names.update(new_names)

        pyfastaq.utils.close(f)
        return clusters


    def run_get_clusters_from_file(self, clusters_infile, all_ref_seqs, rename_dict=None):
        '''Instead of running cdhit, gets the clusters info from the input file.'''
        if rename_dict is None:
            rename_dict = {}

        # check that every sequence in the clusters file can be
        # found in the fasta file
        seq_reader = pyfastaq.sequences.file_reader(self.infile)
        names_list_from_fasta_file = [seq.id for seq in seq_reader]
        names_set_from_fasta_file = set(names_list_from_fasta_file)

        clusters = self._load_user_clusters_file(clusters_infile, all_ref_seqs, rename_dict=rename_dict)

        if len(names_set_from_fasta_file) != len(names_list_from_fasta_file):
            raise Error('At least one duplicate name in fasta file ' + self.infile + '. Cannot continue')

        names_from_clusters_file = set()
        for new_names in clusters.values():
            names_from_clusters_file.update(new_names)

        if not names_set_from_fasta_file.issubset(names_from_clusters_file):
            raise Error('Some names in fasta file "' + self.infile + '" not given in cluster file. Cannot continue')

        return clusters


    @staticmethod
    def _get_clusters_from_bak_file(filename, min_cluster_number=0):
        f = pyfastaq.utils.open_file_read(filename)
        clusters = {}

        for line in f:
            try:
                cluster_number, length, name, *spam = line.rstrip().split()
                cluster_number = int(cluster_number) + min_cluster_number
            except:
                pyfastaq.utils.close(f)
                raise Error('Error parsing cdhit output at this line:\n' + line)

            # keep cluster names as strings in case we want to change them
            # at a later date.
            cluster = str(cluster_number)

            if not (name.startswith('>') and name.endswith('...')):
                pyfastaq.utils.close(f)
                raise Error('Error getting sequence name from cdhit output at this line:\n' + line)

            if cluster not in clusters:
                clusters[cluster] = set()
            clusters[str(cluster)].add(name[1:-3])

        pyfastaq.utils.close(f)
        return clusters


    def run(self):
        tmpdir = tempfile.mkdtemp(prefix='tmp.run_cd-hit.', dir=os.getcwd())
        cdhit_fasta = os.path.join(tmpdir, 'cdhit')
        cluster_info_outfile = cdhit_fasta + '.bak.clstr'

        cmd = ' '.join([
            self.cd_hit_est,
            '-i', self.infile,
            '-o', cdhit_fasta,
            '-c', str(self.seq_identity_threshold),
            '-T', str(self.threads),
            '-s', str(self.length_diff_cutoff),
            '-d 0',
            '-bak 1',
        ])

        common.syscall(cmd, verbose=self.verbose)
        clusters = self._get_clusters_from_bak_file(cluster_info_outfile, self.min_cluster_number)
        shutil.rmtree(tmpdir)
        return clusters

