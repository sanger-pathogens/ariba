import os
import tempfile
import shutil

from ariba import common, external_progs

class Error (Exception): pass

class ReadFilter:
    def __init__(self,
        readstore_obj,
        references_fa,
        cluster_name,
        log_fh,
        extern_progs=None,
    ):
        self.readstore = readstore_obj
        self.references_fa = references_fa
        self.cluster_name = cluster_name
        self.log_fh = log_fh

        if extern_progs is None:
            self.extern_progs = external_progs.ExternalProgs()
        else:
            self.extern_progs = extern_progs


    @staticmethod
    def _run_cdhit_est_2d(reference, reads, outfile, cdhitest2d, verbose=False, verbose_fh=None):
        cmd = ' '.join([
            cdhitest2d,
            '-i', reference,
            '-i2', reads,
            '-G 0 -M 0 -d 0 -aS 0.95',
            '-o', outfile
        ])
        common.syscall(cmd, verbose=verbose, verbose_filehandle=verbose_fh)
        os.unlink(outfile)


    @staticmethod
    def _cdhit_clstr_to_reads(infile):
        reads = set()

        with open(infile) as f:
            for line in f:
                if line.startswith('>') or line.endswith('*\n'):
                    continue

                fields = line.split()
                if not (fields[2].startswith('>') and fields[2].endswith('...')):
                    raise Error('Error parsing the following line of cd-hit-est-2d output from file "' + infile + '":\n' + line)

                reads.add(int(fields[2][1:-5]))

        return reads


    def run(self, reads_out1, reads_out2):
        tmpdir = tempfile.mkdtemp(prefix='tmp.filter_reads.', dir=os.getcwd())
        all_reads_fasta = os.path.join(tmpdir, 'all_reads_for_cdhit.fa')
        self.readstore.get_reads(self.cluster_name, all_reads_fasta, fasta=True, log_fh=self.log_fh)
        cdhit_out = os.path.join(tmpdir, 'cdhit')
        ReadFilter._run_cdhit_est_2d(
            self.references_fa,
            all_reads_fasta,
            cdhit_out,
            self.extern_progs.exe('cdhit2d'),
            verbose=True,
            verbose_fh=self.log_fh
        )

        wanted_read_ids = ReadFilter._cdhit_clstr_to_reads(cdhit_out + '.clstr')
        total_reads, total_bases = self.readstore.get_reads(
            self.cluster_name,
            reads_out1,
            out2=reads_out2,
            log_fh=self.log_fh,
            wanted_ids=wanted_read_ids
        )

        shutil.rmtree(tmpdir)
        return total_reads, total_bases
