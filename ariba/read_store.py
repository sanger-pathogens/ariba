import pysam
import pyfastaq
import os
import pysam
from ariba import common

class Error (Exception): pass

class ReadStore:
    def __init__(self, infile, outfile, log_fh=None):
        self.infile = os.path.abspath(infile)
        self.outfile = os.path.abspath(outfile)
        self.log_fh = log_fh

        if not os.path.exists(self.infile):
            raise Error('File not found ' + self.infile + '. Cannot continue')


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
