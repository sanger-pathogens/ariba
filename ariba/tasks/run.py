import argparse
import os
import sys
import pyfastaq
import ariba

def get_ref_files(options):
    if options.ref_prefix is not None:
        if options.verbose:
            print('--ref_prefix used. Looking for reference input files starting with', options.ref_prefix, '...')
        d = {
            'presabs': 'presence_absence.fa',
            'varonly': 'variants_only.fa',
            'noncoding': 'noncoding.fa',
            'metadata': 'metadata.tsv',
        }

        for key in d:
            filename = options.ref_prefix + '.' + d[key]

            if os.path.exists(filename):
                if options.verbose:
                    print('Found: ', filename, '.\n    ...treating it as if this was used: --', key, ' ', filename, sep='')
                exec('options.' + key + ' = filename')
            else:
                if options.verbose:
                    print('Not found:', filename)
                exec('options.' + key + ' = None')


def run():
    parser = argparse.ArgumentParser(
        description = 'ARIBA: Antibiotic Resistance Identification By Assembly',
        usage = 'ariba run [options] <reads1.fq> <reads2.fq> <outdir>')
    parser.add_argument('reads_1', help='Name of fwd reads fastq file')
    parser.add_argument('reads_2', help='Name of rev reads fastq file')
    parser.add_argument('outdir', help='Output directory (must not already exist)')

    refdata_group = parser.add_argument_group('Reference data options')
    refdata_group.add_argument('--ref_prefix', help='Prefix of input files (same as was used with getref), to save listing --preseabs,--varonly ...etc. Will look for files called ref_prefix. followed by: metadata.tsv,presence_absence.fa,noncoding.fa,presence_absence.fa. Using this will cause these to be ignored if used: --presabs,--varonly,--noncoding,--metadata')
    refdata_group.add_argument('--presabs', help='FASTA file of presence absence genes', metavar='FILENAME')
    refdata_group.add_argument('--varonly', help='FASTA file of variants only genes', metavar='FILENAME')
    refdata_group.add_argument('--noncoding', help='FASTA file of noncoding sequences', metavar='FILENAME')
    refdata_group.add_argument('--metadata', help='tsv file of metadata about the reference sequences', metavar='FILENAME')
    refdata_group.add_argument('--min_gene_length', type=int, help='Minimum allowed length in nucleotides of reference genes [%(default)s]', metavar='INT', default=6)
    refdata_group.add_argument('--max_gene_length', type=int, help='Maximum allowed length in nucleotides of reference genes [%(default)s]', metavar='INT', default=10000)

    cdhit_group = parser.add_argument_group('cd-hit options')
    cdhit_group.add_argument('--no_cdhit', action='store_true', help='Do not run cd-hit')
    cdhit_group.add_argument('--cdhit_seq_identity_threshold', type=float, help='Sequence identity threshold (cd-hit option -c) [%(default)s]', default=0.9, metavar='FLOAT')
    cdhit_group.add_argument('--cdhit_length_diff_cutoff', type=float, help='length difference cutoff (cd-hit option -s) [%(default)s]', default=0.9, metavar='FLOAT')

    nucmer_group = parser.add_argument_group('nucmer options')
    nucmer_group.add_argument('--nucmer_min_id', type=int, help='Minimum alignment identity (delta-filter -i) [%(default)s]', default=90, metavar='INT')
    nucmer_group.add_argument('--nucmer_min_len', type=int, help='Minimum alignment length (delta-filter -i) [%(default)s]', default=50, metavar='INT')
    nucmer_group.add_argument('--nucmer_breaklen', type=int, help='Value to use for -breaklen when running nucmer [%(default)s]', default=50, metavar='INT')

    assembly_group = parser.add_argument_group('Assembly options')
    assembly_group.add_argument('--assembly_cov', type=int, help='Target read coverage when sampling reads for assembly [%(default)s]', default=100, metavar='INT')
    assembly_group.add_argument('--assembler_k', type=int, help='kmer size to use with assembler. You can use 0 to set kmer to 2/3 of the read length. Warning - lower kmers are usually better. [%(default)s]', metavar='INT', default=21)
    assembly_group.add_argument('--spades_other', help='Put options string to be used with spades in quotes. This will NOT be sanity checked. Do not use -k or -t: for these options you should use the ariba run options --assembler_k and --threads [%(default)s]', default="--only-assembler", metavar="OPTIONS")
    assembly_group.add_argument('--min_scaff_depth', type=int, help='Minimum number of read pairs needed as evidence for scaffold link between two contigs. This is also the value used for sspace -k when scaffolding [%(default)s]', default=10, metavar='INT')

    other_group = parser.add_argument_group('Other options')
    other_group.add_argument('--genetic_code', type=int, help='Number of genetic code to use. Currently supported 1,4,11 [%(default)s]', choices=[1,4,11], default=11, metavar='INT')
    other_group.add_argument('--threads', type=int, help='Number of threads [%(default)s]', default=1, metavar='INT')
    bowtie2_presets = ['very-fast-local', 'fast-local', 'sensitive-local', 'very-sensitive-local']
    other_group.add_argument('--bowtie2_preset', choices=bowtie2_presets, help='Preset option for bowtie2 mapping [%(default)s]', default='very-sensitive-local', metavar='|'.join(bowtie2_presets))
    other_group.add_argument('--assembled_threshold', type=float, help='If proportion of gene assembled (regardless of into how many contigs) is at least this value then the flag gene_assembled is set [%(default)s]', default=0.95, metavar='FLOAT (between 0 and 1)')
    other_group.add_argument('--unique_threshold', type=float, help='If proportion of bases in gene assembled more than once is <= this value, then the flag unique_contig is set [%(default)s]', default=0.03, metavar='FLOAT (between 0 and 1)')
    other_group.add_argument('--clean', type=int, choices=[0,1,2], help='Specify how much cleaning to do. 0=none, 1=some, 2=only keep the report [%(default)s]', default=1, metavar='INT')
    other_group.add_argument('--verbose', action='store_true', help='Be verbose')

    options = parser.parse_args()


    if options.verbose:
        print('{:_^79}'.format(' Reference files '), flush=True)
    get_ref_files(options)

    if {None} == {options.presabs, options.varonly, options.noncoding}:
        print('Error! Must use at least one of the options: --presabs --varonly --noncoding. Alternatively, use the option --ref_prefix. Cannot continue', file=sys.stderr)
        sys.exit(1)

    if options.verbose:
        print('\nUsing the following reference files:')
        print('Presence/absence (--presabs):', options.presabs)
        print('Variants only    (--varonly):', options.varonly)
        print('Non coding     (--noncoding):', options.noncoding)
        print('Metadata        (--metadata):', options.metadata)
        print()

    extern_progs = ariba.external_progs.ExternalProgs(verbose=options.verbose)
    pyfastaq.sequences.genetic_code = options.genetic_code

    if options.verbose:
        print('{:_^79}'.format(' Loading reference data '), flush=True)
    refdata = ariba.reference_data.ReferenceData(
        presence_absence_fa=options.presabs,
        variants_only_fa=options.varonly,
        non_coding_fa=options.noncoding,
        metadata_tsv=options.metadata,
        min_gene_length=options.min_gene_length,
        max_gene_length=options.max_gene_length,
        genetic_code=options.genetic_code,
    )

    c = ariba.clusters.Clusters(
          refdata,
          options.reads_1,
          options.reads_2,
          options.outdir,
          extern_progs,
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
          bowtie2_preset=options.bowtie2_preset,
          cdhit_seq_identity_threshold=options.cdhit_seq_identity_threshold,
          cdhit_length_diff_cutoff=options.cdhit_length_diff_cutoff,
          clean=options.clean,
          run_cd_hit=(not options.no_cdhit)
        )
    c.run()

