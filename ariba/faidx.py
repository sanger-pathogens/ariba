import sys
import os
from ariba import common


def write_fa_subset(seq_names, infile, outfile, samtools_exe='samtools', verbose=False, verbose_filehandle=sys.stdout):
    if not os.path.exists(infile + '.fai'):
        common.syscall(samtools_exe + ' faidx ' + infile, verbose=verbose, verbose_filehandle=verbose_filehandle)

    if os.path.exists(outfile):
        os.path.unlink(outfile)

    for name in seq_names:
        common.syscall(' '.join([
            samtools_exe + ' faidx',
            infile,
            '"' + name + '"',
            '>>', outfile
        ]))

