import sys
import argparse
from ariba import ref_preparer, external_progs, versions

def run():
    parser = argparse.ArgumentParser(
        description = 'ARIBA: Antibiotic Resistance Identification By Assembly',
        usage = 'ariba prepareref [options] <outdir>',
        epilog = 'REQUIRED: either --ref_prefix, or at least one of --presabs, --varonly, --noncoding')
    input_group = parser.add_argument_group('input files options')
    input_group.add_argument('--ref_prefix', help='Prefix of input files (same as was used with getref), to save listing --preseabs,--varonly ...etc. Will look for files called "ref_prefix." followed by: metadata.tsv,presence_absence.fa,noncoding.fa,variants_only.fa. Using this will cause these to be ignored if used: --presabs,--varonly,--noncoding,--metadata', metavar='FILENAME_PREFIX')
    input_group.add_argument('--presabs', help='FASTA file of presence absence genes', metavar='FILENAME')
    input_group.add_argument('--varonly', help='FASTA file of variants only genes', metavar='FILENAME')
    input_group.add_argument('--noncoding', help='FASTA file of noncoding sequences', metavar='FILENAME')
    input_group.add_argument('--metadata', help='tsv file of metadata about the reference sequences', metavar='FILENAME')

    cdhit_group = parser.add_argument_group('cd-hit options')
    cdhit_group.add_argument('--no_cdhit', action='store_true', help='Do not run cd-hit. Each input sequence is put into its own "cluster". Incompatible with --cdhit_clusters.')
    cdhit_group.add_argument('--cdhit_clusters', help='File specifying how the sequences should be clustered. Will be used instead of running cdhit. Format is one cluster per line. Sequence names separated by whitespace. First name in line is the cluster representative. Incompatible with --no_cdhit', metavar='FILENAME')
    cdhit_group.add_argument('--cdhit_min_id', type=float, help='Sequence identity threshold (cd-hit option -c) [%(default)s]', default=0.9, metavar='FLOAT')
    cdhit_group.add_argument('--cdhit_min_length', type=float, help='length difference cutoff (cd-hit option -s) [%(default)s]', default=0.9, metavar='FLOAT')

    other_group = parser.add_argument_group('other options')
    other_group.add_argument('--min_gene_length', type=int, help='Minimum allowed length in nucleotides of reference genes [%(default)s]', metavar='INT', default=6)
    other_group.add_argument('--max_gene_length', type=int, help='Maximum allowed length in nucleotides of reference genes [%(default)s]', metavar='INT', default=10000)
    other_group.add_argument('--genetic_code', type=int, help='Number of genetic code to use. Currently supported 1,4,11 [%(default)s]', choices=[1,4,11], default=11, metavar='INT')
    other_group.add_argument('--threads', type=int, help='Number of threads (currently only applies to cdhit) [%(default)s]', default=1, metavar='INT')
    other_group.add_argument('--verbose', action='store_true', help='Be verbose')

    parser.add_argument('outdir', help='Output directory (must not already exist)')
    options = parser.parse_args()

    if options.no_cdhit and options.cdhit_clusters is not None:
        sys.exit('Cannot use both --no_cdhit and --cdhit_clusters. Neither or exactly one of those options must be used')

    extern_progs, version_report_lines = versions.get_all_versions()
    if options.verbose:
        print(*version_report_lines, sep='\n')

    preparer = ref_preparer.RefPreparer(
        extern_progs,
        version_report_lines=version_report_lines,
        ref_prefix=options.ref_prefix,
        presabs=options.presabs,
        varonly=options.varonly,
        noncoding=options.noncoding,
        metadata=options.metadata,
        min_gene_length=options.min_gene_length,
        max_gene_length=options.max_gene_length,
        genetic_code=options.genetic_code,
        cdhit_min_id=options.cdhit_min_id,
        cdhit_min_length=options.cdhit_min_length,
        run_cdhit=not options.no_cdhit,
        clusters_file=options.cdhit_clusters,
        threads=options.threads,
        verbose=options.verbose,
    )

    preparer.run(options.outdir)
