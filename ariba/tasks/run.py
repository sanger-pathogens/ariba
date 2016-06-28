import argparse
import os
import sys
import ariba


def run():
    parser = argparse.ArgumentParser(
        description = 'ARIBA: Antibiotic Resistance Identification By Assembly',
        usage = 'ariba run [options] <prepareref_dir> <reads1.fq> <reads2.fq> <outdir>')
    parser.add_argument('prepareref_dir', help='Name of output directory when "ariba prepareref" was run')
    parser.add_argument('reads_1', help='Name of fwd reads fastq file')
    parser.add_argument('reads_2', help='Name of rev reads fastq file')
    parser.add_argument('outdir', help='Output directory (must not already exist)')

    nucmer_group = parser.add_argument_group('nucmer options')
    nucmer_group.add_argument('--nucmer_min_id', type=int, help='Minimum alignment identity (delta-filter -i) [%(default)s]', default=90, metavar='INT')
    nucmer_group.add_argument('--nucmer_min_len', type=int, help='Minimum alignment length (delta-filter -i) [%(default)s]', default=20, metavar='INT')
    nucmer_group.add_argument('--nucmer_breaklen', type=int, help='Value to use for -breaklen when running nucmer [%(default)s]', default=200, metavar='INT')

    assembly_group = parser.add_argument_group('Assembly options')
    assembly_group.add_argument('--assembly_cov', type=int, help='Target read coverage when sampling reads for assembly [%(default)s]', default=50, metavar='INT')
    assembly_group.add_argument('--assembler_k', type=int, help='kmer size to use with assembler. You can use 0 to set kmer to 2/3 of the read length. Warning - lower kmers are usually better. [%(default)s]', metavar='INT', default=21)
    assembly_group.add_argument('--spades_other', help='Put options string to be used with spades in quotes. This will NOT be sanity checked. Do not use -k (see --assembler_k), or -t (use ariba option --threads) [%(default)s]', default="--only-assembler -m 4", metavar="OPTIONS")
    assembly_group.add_argument('--min_scaff_depth', type=int, help='Minimum number of read pairs needed as evidence for scaffold link between two contigs. This is also the value used for sspace -k when scaffolding [%(default)s]', default=10, metavar='INT')

    other_group = parser.add_argument_group('Other options')
    other_group.add_argument('--threads', type=int, help='Number of threads [%(default)s]', default=1, metavar='INT')
    bowtie2_presets = ['very-fast-local', 'fast-local', 'sensitive-local', 'very-sensitive-local']
    other_group.add_argument('--bowtie2_preset', choices=bowtie2_presets, help='Preset option for bowtie2 mapping [%(default)s]', default='very-sensitive-local', metavar='|'.join(bowtie2_presets))
    other_group.add_argument('--assembled_threshold', type=float, help='If proportion of gene assembled (regardless of into how many contigs) is at least this value then the flag gene_assembled is set [%(default)s]', default=0.95, metavar='FLOAT (between 0 and 1)')
    other_group.add_argument('--gene_nt_extend', type=int, help='Max number of nucleotides to extend ends of gene matches to look for start/stop codons [%(default)s]', default=30, metavar='INT')
    other_group.add_argument('--unique_threshold', type=float, help='If proportion of bases in gene assembled more than once is <= this value, then the flag unique_contig is set [%(default)s]', default=0.03, metavar='FLOAT (between 0 and 1)')
    other_group.add_argument('--noclean', action='store_true', help='Do not clean up intermediate files')
    other_group.add_argument('--tmp_dir', help='Existing directory in which to create a temporary directory used for local assemblies')
    other_group.add_argument('--verbose', action='store_true', help='Be verbose')

    options = parser.parse_args()

    reads_not_found = []

    for filename in [options.reads_1, options.reads_2]:
        if not os.path.exists(filename):
            reads_not_found.append(filename)
        elif options.verbose:
            print('Found reads file:', filename)

    if len(reads_not_found):
        print('\nThe following reads file(s) were not found:', file=sys.stderr)
        print(*reads_not_found, sep='\n', file=sys.stderr)
        print('Cannot continue', file=sys.stderr)
        sys.exit(1)

    if not os.path.exists(options.prepareref_dir):
        print('Input directory', options.prepareref_dir, 'not found. Cannot continue', file=sys.stderr)
        sys.exit(1)

    extern_progs, version_report_lines = ariba.versions.get_all_versions()
    if options.verbose:
        print(*version_report_lines, sep='\n')

    c = ariba.clusters.Clusters(
          options.prepareref_dir,
          options.reads_1,
          options.reads_2,
          options.outdir,
          extern_progs,
          version_report_lines=version_report_lines,
          assembly_kmer=options.assembler_k,
          assembly_coverage=options.assembly_cov,
          assembler='spades',
          threads=options.threads,
          verbose=options.verbose,
          min_scaff_depth=options.min_scaff_depth,
          nucmer_min_id=options.nucmer_min_id,
          nucmer_min_len=options.nucmer_min_len,
          nucmer_breaklen=options.nucmer_breaklen,
          spades_other=options.spades_other,
          assembled_threshold=options.assembled_threshold,
          unique_threshold=options.unique_threshold,
          max_gene_nt_extend=options.gene_nt_extend,
          bowtie2_preset=options.bowtie2_preset,
          clean=(not options.noclean),
          tmp_dir=options.tmp_dir,
        )
    c.run()

