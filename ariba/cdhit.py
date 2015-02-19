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


    def run(self):
        tmpdir = tempfile.mkdtemp(prefix='tmp.run_cd-hit.', dir=os.getcwd())
        cdhit_fasta = os.path.join(tmpdir, 'cdhit')
        cluster_info_outfile = cdhit_fasta + '.bak.clstr'
        infile_renamed = os.path.join(tmpdir, 'input.renamed.fa')

        # cd-hit truncates all names to 19 bases in its report of which
        # sequences belong to which clusters. So need to temporarily
        # rename all sequences to have short enough names. Grrr.
        new_to_old_name = self._enumerate_fasta(self.infile, infile_renamed)

        cmd = ' '.join([
            'cd-hit',
            '-i', infile_renamed,
            '-o', cdhit_fasta,
            '-c', str(self.seq_identity_threshold),
            '-T', str(self.threads),
            '-s', str(self.length_diff_cutoff),
            '-bak 1',
        ])

        common.syscall(cmd, verbose=self.verbose) 

        cluster_representatives = self._get_ids(cdhit_fasta)
        clusters, cluster_rep_to_cluster = self._parse_cluster_info_file(cluster_info_outfile, new_to_old_name, cluster_representatives)
        self._rename_fasta(cdhit_fasta, self.outfile, cluster_rep_to_cluster)
        shutil.rmtree(tmpdir)
        return clusters


    def _enumerate_fasta(self, infile, outfile):
        rename_file = outfile + '.tmp.rename_info'
        assert not os.path.exists(rename_file)
        pyfastaq.tasks.enumerate_names(infile, outfile, rename_file=rename_file)
        
        with open(rename_file) as f:
            lines = [x.rstrip().split('\t') for x in f.readlines() if x != '#old\tnew\n']
            new_to_old_name = {x[1]: x[0] for x in lines}
            if len(lines) != len(new_to_old_name):
                raise Error('Sequence names in input file not unique! Cannot continue')

        os.unlink(rename_file)
        return new_to_old_name


    def _rename_fasta(self, infile, outfile, names_dict):
        seq_reader = pyfastaq.sequences.file_reader(infile)
        f = pyfastaq.utils.open_file_write(outfile)
        for seq in seq_reader:
            seq.id = names_dict[seq.id]
            print(seq, file=f)

        pyfastaq.utils.close(f)


    def _parse_cluster_info_file(self, infile, names_dict, cluster_representatives):
        f = pyfastaq.utils.open_file_read(infile)
        clusters = {}
        cluster_representative_to_cluster_number = {}
        for line in f:
            data = line.rstrip().split()
            cluster = data[0]
            seqname = data[2]
            if not (seqname.startswith('>') and seqname.endswith('...')):
                raise Error('Unexpected format of sequence name in line:\n' + line)
            seqname = seqname[1:-3]

            if seqname in cluster_representatives:
                cluster_representative_to_cluster_number[seqname] = cluster

            seqname = names_dict[seqname]

            if cluster not in clusters:
                clusters[cluster] = set()

            if seqname in clusters[cluster]:
                raise Error('Duplicate name "' + seqname + '" found in cluster ' + str(cluster))

            clusters[cluster].add(seqname)

        pyfastaq.utils.close(f)

        return clusters, cluster_representative_to_cluster_number


    def _get_ids(self, infile):
        seq_reader = pyfastaq.sequences.file_reader(infile)
        return set([seq.id for seq in seq_reader])

