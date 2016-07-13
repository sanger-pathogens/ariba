import tempfile
import shutil
import os
import pyfastaq
from ariba import common, external_progs

class Error (Exception): pass

class Runner:
    def __init__(
      self,
      infile,
      outfile,
      seq_identity_threshold=0.9,
      threads=1,
      length_diff_cutoff=0.9,
      verbose=False,
      rename_suffix='x',
      min_cluster_number=0
    ):

        if not os.path.exists(infile):
            raise Error('File not found: "' + infile + '". Cannot continue')

        self.infile = os.path.abspath(infile)
        self.outfile = os.path.abspath(outfile)
        self.seq_identity_threshold = seq_identity_threshold
        self.threads = threads
        self.length_diff_cutoff = length_diff_cutoff
        self.verbose = verbose
        extern_progs = external_progs.ExternalProgs(fail_on_error=True)
        self.cd_hit_est = extern_progs.exe('cdhit')
        self.rename_suffix = rename_suffix
        self.min_cluster_number = min_cluster_number


    def fake_run(self):
        '''Doesn't actually run cd-hit. Instead, puts each input sequence into its own cluster. So it's as if cdhit was run, but didn't cluster anything'''
        clusters = {}
        used_names = set()
        seq_reader = pyfastaq.sequences.file_reader(self.infile)

        for seq in seq_reader:
            if seq.id in used_names:
                raise Error('Sequence name "' + seq.id + '" not unique. Cannot continue')

            clusters[str(len(clusters))] = {seq.id}
            used_names.add(seq.id)

        return clusters


    @staticmethod
    def _load_user_clusters_file(filename):
        f = pyfastaq.utils.open_file_read(filename)
        seq_to_cluster = {}
        for line in f:
            data = line.rstrip().split()

            for seq_name in data:
                if seq_name in seq_to_cluster:
                    pyfastaq.utils.close(f)
                    raise Error('Error reading clusters file. The sequence "' + seq_name + '" was found more than once in the file ' + filename)
                seq_to_cluster[seq_name] = data[0]

        pyfastaq.utils.close(f)
        return seq_to_cluster


    def run_get_clusters_from_file(self, infile):
        '''Instead of running cdhit, gets the clusters info from the input dict.
           Dict expected to be key=sequence name, value=name of cluster'''
        seq_to_cluster = self._load_user_clusters_file(infile)
        cluster_names = set(seq_to_cluster.values())
        tmpdir = tempfile.mkdtemp(prefix='tmp.run_cd-hit.', dir=os.getcwd())
        tmp_fa = os.path.join(tmpdir, 'cdhit.fa')
        clusters = {}
        seq_reader = pyfastaq.sequences.file_reader(self.infile)
        f = pyfastaq.utils.open_file_write(tmp_fa)

        for seq in seq_reader:
            if seq.id in clusters and seq.id in clusters[seq.id]:
                pyfastaq.utils.close(f)
                shutil.rmtree(tmpdir)
                raise Error('Sequence name "' + seq.id + '" not unique. Cannot continue')

            if seq.id not in seq_to_cluster:
                raise Error('Error forcing cdhit clustering. Found sequence ' + seq.id + ' in FASTA file, but not in provided clusters info from file ' + infile)

            cluster = seq_to_cluster[seq.id]
            if cluster not in clusters:
                clusters[cluster] = set()

            clusters[cluster].add(seq.id)
            if seq.id in cluster_names:
                print(seq, file=f)

        pyfastaq.utils.close(f)
        clusters = self._rename_clusters(clusters, tmp_fa, self.outfile, rename_suffix=self.rename_suffix)
        shutil.rmtree(tmpdir)
        return clusters


    def _get_ids(self, infile):
        seq_reader = pyfastaq.sequences.file_reader(infile)
        return set([seq.id for seq in seq_reader])


    @staticmethod
    def _rename_clusters(clusters_dict, infile, outfile, rename_suffix='x'):
        new_clusters_dict = {}
        freader = pyfastaq.sequences.file_reader(infile)
        f_out = pyfastaq.utils.open_file_write(outfile)

        for seq in freader:
            original_name = seq.id
            assert original_name in clusters_dict
            new_name = original_name.split('.')[0] + '.' + rename_suffix

            if new_name in new_clusters_dict:
                suffix = 2
                while new_name + '.' + str(suffix) in new_clusters_dict:
                    suffix += 1
                new_name += '.' + str(suffix)

            new_clusters_dict[new_name] = clusters_dict[original_name]
            seq.id = new_name
            print(seq, file=f_out)

        pyfastaq.utils.close(f_out)

        return new_clusters_dict


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

