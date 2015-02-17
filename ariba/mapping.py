import os
from ariba import common

class Error (Exception): pass


def run_smalt(
      reads_fwd,
      reads_rev,
      ref_fa,
      out_prefix,
      index_k=9,
      index_s=2,
      threads=1,
      max_insert=1000,
      minid=0.9,
      sort=False,
      extra_smalt_map_ops='-x',
      samtools='samtools',
      smalt='smalt',
      verbose=False
    ):
    if extra_smalt_map_ops is None:
        extra_smalt_map_ops = ''
    map_index = out_prefix + '.map_index'
    clean_files = [map_index + '.' + x for x in ['smi', 'sma']]
    index_cmd = ' '.join([
        smalt, 'index',
        '-k', str(index_k),
        '-s', str(index_s),
        map_index,
        ref_fa
    ])

    map_cmd = smalt + ' map ' + extra_smalt_map_ops + ' '

    # depending on OS, -n can break smalt, so only use -n if it's > 1.
    if threads > 1:
        map_cmd += '-n ' + str(threads) + ' -O '

    if reads_rev is None:
        map_cmd += ' '.join([
            '-y', str(minid),
            map_index,
            reads_fwd,
        ])
    else:
        map_cmd += ' '.join([
            '-i', str(max_insert),
            '-y', str(minid),
            map_index,
            reads_fwd,
            reads_rev,
        ])

    map_cmd += ' | ' + samtools + ' view'

    final_bam = out_prefix + '.bam'
    if sort:
        intermediate_bam = out_prefix + '.unsorted.bam'
    else:
        intermediate_bam = final_bam

    map_cmd += ' -bS -T ' + ref_fa + '  - > ' + intermediate_bam
    common.syscall(index_cmd, verbose=verbose)
    common.syscall(map_cmd, verbose=verbose)

    if sort:
        threads = min(4, threads)
        thread_mem = int(500 / threads)
        sort_cmd = samtools + ' sort -@' + str(threads) + ' -m ' + str(thread_mem) + 'M ' + intermediate_bam + ' ' + out_prefix
        index_cmd = samtools + ' index ' + final_bam
        common.syscall(sort_cmd, verbose=verbose)
        common.syscall(index_cmd, verbose=verbose)
    for fname in clean_files:
        os.unlink(fname)
