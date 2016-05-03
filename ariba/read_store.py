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


    def get_reads(self, cluster_name, out1, out2, log_fh=None):
        if log_fh is not None:
            print('Getting reads for', cluster_name, 'from', self.outfile, file=log_fh)
        tabix_file = pysam.TabixFile(self.outfile)
        f_out1 = pyfastaq.utils.open_file_write(out1)
        f_out2 = pyfastaq.utils.open_file_write(out2)

        for line in tabix_file.fetch(reference=cluster_name):
            cluster, number, seq, qual = line.rstrip().split()
            number = int(number)
            if number % 2 == 0:
                print('@' + str(number - 1) + '/2', seq, '+', qual, sep='\n', file=f_out2)
            else:
                print('@' + str(number) + '/1', seq, '+', qual, sep='\n', file=f_out1)

        pyfastaq.utils.close(f_out1)
        pyfastaq.utils.close(f_out2)
        if log_fh is not None:
            print('Finished getting reads for', cluster_name, 'from', self.outfile, file=log_fh)

    def clean(self):
        os.unlink(self.outfile)
        os.unlink(self.outfile + '.tbi')
