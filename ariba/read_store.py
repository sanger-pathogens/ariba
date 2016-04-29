import pysam
import pyfastaq
import os
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
