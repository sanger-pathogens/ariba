import sys
import os
from ariba import common, external_progs

class Masher:
    def __init__(self,
        reference_fa,
        query_fa,
        log_fh=sys.stdout,
        extern_progs=None,
    ):
        self.reference_fa = reference_fa
        self.query_fa = query_fa
        self.log_fh = log_fh

        if extern_progs is None:
            self.extern_progs = external_progs.ExternalProgs()
        else:
            self.extern_progs = extern_progs


    @classmethod
    def sketch(cls, infile, individual, extern_progs, verbose=True, verbose_filehandle=None):
        if verbose:
            assert verbose_filehandle is not None

        cmd_list = [
            extern_progs.exe('mash'),
            'sketch',
            '-s 100000'
        ]

        if individual:
            cmd_list.append('-i')

        cmd_list.append(infile)
        common.syscall(' '.join(cmd_list), verbose=verbose, verbose_filehandle=verbose_filehandle)


    def _dist(self, outfile):
        cmd = ' '.join([
            self.extern_progs.exe('mash'),
            'dist',
            self.reference_fa + '.msh',
            self.query_fa + '.msh',
            '| sort -k3n >', outfile
        ])
        common.syscall(cmd, verbose=True, verbose_filehandle=self.log_fh)


    def run(self, outfile):
        if not os.path.exists(self.reference_fa + '.msh'):
            Masher.sketch(self.reference_fa, True, self.extern_progs, verbose=True, verbose_filehandle=self.log_fh)
        Masher.sketch(self.query_fa, False, self.extern_progs, verbose=True, verbose_filehandle=self.log_fh)
        self._dist(outfile)
        if os.path.getsize(outfile) == 0:
            return None

        try:
            with open(outfile) as f:
                line = f.readline().rstrip()
                print('best mash match:\t', line, file=self.log_fh)
                matching_hashes = int(line.split()[-1].split('/')[0])
                if matching_hashes == 0:
                    best_seq = None
                else:
                    best_seq = line.split()[0]
        except:
            best_seq = None

        return best_seq

