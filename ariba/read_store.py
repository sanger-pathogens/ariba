import pyfastaq
import os
import pysam
from ariba import common

class Error (Exception): pass

class ReadStore:
    def __init__(self, infile, outprefix, log_fh=None):
        assert infile != outprefix
        self.infile = os.path.abspath(infile)
        self.outprefix = os.path.abspath(outprefix)
        self.outfile = os.path.abspath(outprefix) + '.gz'

        if not os.path.exists(self.infile):
            raise Error('File not found ' + self.infile + '. Cannot continue')

        self._sort_file(self.infile, self.outprefix, log_fh)
        self._compress_and_index_file(self.outprefix, log_fh)
        os.unlink(self.outprefix)


    @staticmethod
    def _sort_file(infile, outfile, log_fh=None):
        cmd = 'sort -k1,1 -k 2,2n ' + infile + ' > ' + outfile
        verbose = log_fh is not None
        common.syscall(cmd, verbose=verbose, verbose_filehandle=log_fh)


    @staticmethod
    def _compress_and_index_file(infile, log_fh=None):
        if log_fh is not None:
            print('Compressing file', infile, file=log_fh, flush=True)
        pysam.tabix_compress(infile, infile + '.gz')
        pysam.tabix_index(infile + '.gz', seq_col=0, start_col=1, end_col=1)


    def get_reads(self, cluster_name, out1, out2=None, fasta=False, log_fh=None, wanted_ids=None):
        total_reads = 0
        total_bases = 0

        if log_fh is not None:
            print('Getting reads for', cluster_name, 'from', self.outfile, file=log_fh)
        tabix_file = pysam.TabixFile(self.outfile)
        f_out1 = pyfastaq.utils.open_file_write(out1)
        if out2 is None:
            f_out2 = f_out1
        else:
            f_out2 = pyfastaq.utils.open_file_write(out2)

        for line in tabix_file.fetch(reference=cluster_name):
            cluster, number, seq, qual = line.rstrip().split()
            number = int(number)
            if wanted_ids is not None:
                new_number = number if number % 2 else number - 1
                if new_number not in wanted_ids:
                    continue

            if number % 2 == 0:
                if fasta:
                    print('>' + str(number - 1) + '/2', seq, sep='\n', file=f_out2)
                else:
                    print('@' + str(number - 1) + '/2', seq, '+', qual, sep='\n', file=f_out2)
            else:
                if fasta:
                    print('>' + str(number) + '/1', seq, sep='\n', file=f_out1)
                else:
                    print('@' + str(number) + '/1', seq, '+', qual, sep='\n', file=f_out1)

            total_reads += 1
            total_bases += len(qual)

        pyfastaq.utils.close(f_out1)
        if out2 is not None:
            pyfastaq.utils.close(f_out2)
        if log_fh is not None:
            print('Finished getting reads for', cluster_name, 'from', self.outfile, file=log_fh)

        return total_reads, total_bases

    def clean(self):
        os.unlink(self.outfile)
        os.unlink(self.outfile + '.tbi')
