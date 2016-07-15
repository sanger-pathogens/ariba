import sys
import argparse
from ariba import ref_preparer, external_progs, versions

def run():
    parser = argparse.ArgumentParser(
        description = 'ARIBA: Antibiotic Resistance Identification By Assembly',
        usage = 'ariba prepareref [options] <outdir>',
        epilog = 'REQUIRED: -f and -m must each be used at least once')
    input_group = parser.add_argument_group('input files options')
    input_group.add_argument('-f', '--fasta', action='append', dest='fasta_files', required=True, help='REQUIRED. Name of fasta file. Can be used more than once if your sequences are spread over more than on file', metavar='FILENAME')
    input_group.add_argument('-m', '--metadata', action='append', dest='tsv_files', required=True, help='REQUIRED. Name of tsv file of metadata about the input sequences. Can be used more than once if your metadata is spread over more than one file', metavar='FILENAME')

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
        options.fasta_files,
        options.tsv_files,
        extern_progs,
        version_report_lines=version_report_lines,
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
