import os
from ariba import common


def write_fa_subset(seq_names, infile, outfile, samtools_exe='samtools', verbose=False):
    if not os.path.exists(infile + '.fai'):
        common.syscall(samtools_exe + ' faidx ' + infile, verbose=verbose)

    if os.path.exists(outfile):
        os.path.unlink(outfile)

    for name in seq_names:
        common.syscall(' '.join([
            samtools_exe + ' faidx',
            infile,
            name,
            '>>', outfile
        ]))

