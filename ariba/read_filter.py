import os
import tempfile

from ariba import common


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

