import tempfile
import shutil
import os
import pyfastaq
from ariba import common

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
      cd_hit_est='cd-hit-est',
      rename_suffix='x',
    ):

        if not os.path.exists(infile):
            raise Error('File not found: "' + infile + '". Cannot continue')

        self.infile = os.path.abspath(infile)
        self.outfile = os.path.abspath(outfile)
        self.seq_identity_threshold = seq_identity_threshold
        self.threads = threads
        self.length_diff_cutoff = length_diff_cutoff
        self.verbose = verbose
        self.cd_hit_est = cd_hit_est
        self.rename_suffix = rename_suffix


    def fake_run(self):
        '''Doesn't actually run cd-hit. Instead, puts each input sequence into its own cluster. So it's as if cdhit was run, but didn't cluster anything'''
        tmpdir = tempfile.mkdtemp(prefix='tmp.run_cd-hit.', dir=os.getcwd())
        tmp_fa = os.path.join(tmpdir, 'cdhit.fa')
        clusters = {}
        seq_reader = pyfastaq.sequences.file_reader(self.infile)
        f = pyfastaq.utils.open_file_write(tmp_fa)

        for seq in seq_reader:
            if seq.id in clusters:
                pyfastaq.utils.close(f)
                shutil.rmtree(tmpdir)
                raise Error('Sequence name "' + seq.id + '" not unique. Cannot continue')

            clusters[seq.id] = {seq.id}
            print(seq, file=f)

        pyfastaq.utils.close(f)
        clusters = self._rename_clusters(clusters, tmp_fa, self.outfile, rename_suffix=self.rename_suffix)
        shutil.rmtree(tmpdir)
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
    def _parse_cluster_info_file(infile, cluster_representatives):
        f = pyfastaq.utils.open_file_read(infile)
        cluster_sets = {}
        found_representatives = {}  # store cluster number -> representative name

        for line in f:
            data = line.rstrip().split()
            seqname = data[2]
            if not (seqname.startswith('>') and seqname.endswith('...')):
                raise Error('Unexpected format of line from cdhit output file "' + infile + '". Line is:\n' + line)
            seqname = seqname[1:-3]

            cluster_number = int(data[0]) # this is the cluster number used by cdhit
            if cluster_number not in cluster_sets:
                cluster_sets[cluster_number] = set()

            cluster_sets[cluster_number].add(seqname)

            if data[3] == '*':
                found_representatives[cluster_number] = seqname

        pyfastaq.utils.close(f)

        if set(found_representatives.values()) != cluster_representatives:
            raise Error('Mismatch in cdhit output sequence names between fasta file and clusters file. Cannot continue')

        clusters = {}
        for cluster_number, cluster_name in found_representatives.items():
            clusters[cluster_name] = cluster_sets[cluster_number]

        return clusters


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
        cluster_representatives = self._get_ids(cdhit_fasta)
        clusters = self._parse_cluster_info_file(cluster_info_outfile, cluster_representatives)
        clusters = self._rename_clusters(clusters, cdhit_fasta, self.outfile, rename_suffix=self.rename_suffix)

        shutil.rmtree(tmpdir)
        return clusters

