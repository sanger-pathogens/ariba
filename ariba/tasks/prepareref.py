import argparse
import os
import sys
from ariba import ref_preparer


def run():
    parser = argparse.ArgumentParser(
        description = 'ARIBA: Antibiotic Resistance Identification By Assembly',
        usage = 'ariba prepareref [options] <outdir>')
    parser.add_argument('--ref_prefix', help='Prefix of input files (same as was used with getref), to save listing --preseabs,--varonly ...etc. Will look for files called ref_prefix. followed by: metadata.tsv,presence_absence.fa,noncoding.fa,presence_absence.fa. Using this will cause these to be ignored if used: --presabs,--varonly,--noncoding,--metadata')
    parser.add_argument('--presabs', help='FASTA file of presence absence genes', metavar='FILENAME')
    parser.add_argument('--varonly', help='FASTA file of variants only genes', metavar='FILENAME')
    parser.add_argument('--noncoding', help='FASTA file of noncoding sequences', metavar='FILENAME')
    parser.add_argument('--metadata', help='tsv file of metadata about the reference sequences', metavar='FILENAME')
    parser.add_argument('--min_gene_length', type=int, help='Minimum allowed length in nucleotides of reference genes [%(default)s]', metavar='INT', default=6)
    parser.add_argument('--max_gene_length', type=int, help='Maximum allowed length in nucleotides of reference genes [%(default)s]', metavar='INT', default=10000)
    parser.add_argument('--genetic_code', type=int, help='Number of genetic code to use. Currently supported 1,4,11 [%(default)s]', choices=[1,4,11], default=11, metavar='INT')
    parser.add_argument('--verbose', action='store_true', help='Be verbose')
    parser.add_argument('outdir', help='Output directory (must not already exist)')
    options = parser.parse_args()

    preparer = ref_preparer.RefPreparer(
        ref_prefix=options.ref_prefix,
        presabs=options.presabs,
        varonly=options.varonly,
        noncoding=options.noncoding,
        metadata=options.metadata,
        min_gene_length=options.min_gene_length,
        max_gene_length=options.max_gene_length,
        genetic_code=options.genetic_code,
        verbose=options.verbose,
    )

    preparer.run(options.outdir)
