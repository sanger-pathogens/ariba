import os
import sys
import pysam
import pyfastaq
from ariba import common

class Error (Exception): pass

bowtie2_index_extensions = [x + '.bt2' for x in ['1', '2', '3', '4', 'rev.1', 'rev.2']]

def bowtie2_index(ref_fa, outprefix, bowtie2='bowtie2', verbose=False, verbose_filehandle=sys.stdout):
    expected_files = [outprefix + '.' + x + '.bt2' for x in ['1', '2', '3', '4', 'rev.1', 'rev.2']]
    file_missing = False
    for filename in expected_files:
        if not os.path.exists(filename):
            file_missing = True
            break

    if not file_missing:
        return

    cmd = ' '.join([
        bowtie2 + '-build',
        '-q',
        ref_fa,
        outprefix
    ])

    common.syscall(cmd, verbose=verbose, verbose_filehandle=verbose_filehandle)


def run_bowtie2(
      reads_fwd,
      reads_rev,
      ref_fa,
      out_prefix,
      threads=1,
      max_insert=1000,
      sort=False,
      bowtie2='bowtie2',
      bowtie2_preset='very-sensitive-local',
      verbose=False,
      verbose_filehandle=sys.stdout,
      remove_both_unmapped=False,
      clean_index=True,
    ):

    ref_is_indexed = True
    for ext in bowtie2_index_extensions:
        if not os.path.exists(ref_fa + '.' + ext):
            ref_is_indexed = False
            break

    clean_files = []

    if ref_is_indexed:
        if verbose:
            print('Bowtie2 index files found (', ref_fa, '.*.bt2) so no need to index', sep='', file=verbose_filehandle)
        map_index = ref_fa
    else:
        map_index = out_prefix + '.map_index'
        bowtie2_index(ref_fa, map_index, bowtie2=bowtie2, verbose=verbose, verbose_filehandle=verbose_filehandle)

        if clean_index:
            clean_files = [map_index + '.' + x + '.bt2' for x in ['1', '2', '3', '4', 'rev.1', 'rev.2']]

    final_bam = out_prefix + '.bam'
    if sort:
        intermediate_bam = out_prefix + '.unsorted.bam'
    else:
        intermediate_bam = final_bam

    map_cmd = [
        bowtie2,
        '--threads', str(threads),
        '--reorder',
        '--' + bowtie2_preset,
        '-X', str(max_insert),
        '-x', map_index,
        '-1', reads_fwd,
        '-2', reads_rev,
    ]

    if remove_both_unmapped:
        map_cmd.append(r''' | awk ' !(and($2,4)) || !(and($2,8)) ' ''')

    tmp_sam_file = out_prefix + '.unsorted.sam'
    map_cmd.append(' > ' + tmp_sam_file)
    map_cmd = ' '.join(map_cmd)

    common.syscall(map_cmd, verbose=verbose, verbose_filehandle=verbose_filehandle)

    if verbose:
        print('Converting', tmp_sam_file, '->', intermediate_bam, file=verbose_filehandle)
    infile = pysam.AlignmentFile(tmp_sam_file, "r")
    outfile = pysam.AlignmentFile(intermediate_bam, "wb", template=infile)
    for x in infile:
        outfile.write(x)
    infile.close()
    outfile.close()
    os.unlink(tmp_sam_file)

    if sort:
        if verbose:
            print('Sorting', intermediate_bam, '->', final_bam, file=verbose_filehandle)
        pysam.sort('-o', final_bam, '-O', 'BAM', intermediate_bam)
        if verbose:
            print('Indexing', final_bam, file=verbose_filehandle)
        pysam.index(final_bam)
        clean_files.append(intermediate_bam)

    for fname in clean_files:
        os.unlink(fname)


def get_total_alignment_score(bam):
    '''Returns total of AS: tags in the input BAM'''
    sam_reader = pysam.Samfile(bam, "rb")
    total = 0
    for sam in sam_reader.fetch(until_eof=True):
        try:
            total += sam.opt('AS')
        except:
            pass
    return total


def sam_to_fastq(sam):
    '''Given a pysam alignment, returns the sequence a Fastq object.
       Reverse complements as required and add suffix /1 or /2 as appropriate from the flag'''
    name = sam.qname
    if sam.is_read1:
        name += '/1'
    elif sam.is_read2:
        name += '/2'
    else:
        raise Error('Read ' + name + ' must be first or second of pair according to flag. Cannot continue')

    seq = pyfastaq.sequences.Fastq(name, common.decode(sam.seq), common.decode(sam.qual))
    if sam.is_reverse:
        seq.revcomp()

    return seq


def sam_pair_to_insert(s1, s2):
    '''Returns insert size from pair of sam records, as long as their orientation is "innies".
       Otherwise returns None.'''
    if s1.is_unmapped or s2.is_unmapped or (s1.tid != s2.tid) or (s1.is_reverse == s2.is_reverse):
        return None

    # If here, reads are both mapped to the same ref, and in opposite orientations
    if s1.is_reverse:
        end = s1.reference_end - 1
        start = s2.reference_start
    else:
        end = s2.reference_end - 1
        start = s1.reference_start

    if start < end:
        return end - start + 1
    else:
        return None


