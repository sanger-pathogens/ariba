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
    ):

        if not os.path.exists(infile):
            raise Error('File not found: "' + infile + '". Cannot continue')

        self.infile = os.path.abspath(infile)
        self.outfile = os.path.abspath(outfile)
        self.seq_identity_threshold = seq_identity_threshold
        self.threads = threads
        self.length_diff_cutoff = length_diff_cutoff
        self.verbose = verbose


    def fake_run(self):
        '''Doesn't actually run cd-hit. Instead, puts each input sequence into its own cluster. So it's as if cdhit was run, but didn't cluster anything'''
        clusters = {}
        seq_reader = pyfastaq.sequences.file_reader(self.infile)
        f = pyfastaq.utils.open_file_write(self.outfile)

        for seq in seq_reader:
            if seq.id in clusters:
                pyfastaq.utils.close(f)
                raise Error('Sequence name "' + seq.id + '" not unique. Cannot continue')

            clusters[seq.id] = {seq.id}
            print(seq, file=f)

        pyfastaq.utils.close(f)
        return clusters


    def _get_ids(self, infile):
        seq_reader = pyfastaq.sequences.file_reader(infile)
        return set([seq.id for seq in seq_reader])


    @staticmethod
    def _parse_cluster_info_file(infile, cluster_representatives):
        f = pyfastaq.utils.open_file_read(infile)
        clusters = {}
        current_cluster = None

        for line in f:
            data = line.rstrip().split()
            seqname = data[2]
            if not (seqname.startswith('>') and seqname.endswith('...')):
                raise Error('Unexpected format of line from cdhit output file "' + infile + '". Line is:\n' + line)
            seqname = seqname[1:-3]

            if data[3] == '*':
                current_cluster = seqname
                assert current_cluster not in clusters
                clusters[current_cluster] = {current_cluster}
            else:
                assert current_cluster in clusters
                if seqname in clusters[current_cluster]:
                    raise Error('Duplicate name "' + seqname + '" found in cluster ' + cluster)

                clusters[current_cluster].add(seqname)

        pyfastaq.utils.close(f)
        if set(clusters.keys()) != cluster_representatives:
            raise Error('Mismatch in cdhit output sequence names between fasta file and clusters file. Cannot continue')

        return clusters


    def run(self):
        tmpdir = tempfile.mkdtemp(prefix='tmp.run_cd-hit.', dir=os.getcwd())
        cdhit_fasta = os.path.join(tmpdir, 'cdhit')
        cluster_info_outfile = cdhit_fasta + '.bak.clstr'

        cmd = ' '.join([
            'cd-hit-est',
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

        try:
            os.rename(cdhit_fasta, self.outfile)
        except:
            raise Error('Error rname ' + cdhit_fasta + ' ' + self.outfile + '. Cannot continue')

        shutil.rmtree(tmpdir)
        return clusters

