import argparse
import ariba

def run():
    parser = argparse.ArgumentParser(
        description = 'Check or fix resistance genes FASTA file',
        usage = 'ariba refcheck [options] <outprefix>',
        epilog = 'Important: at least one of --presence_fa, --variants_fa, or --noncoding_fa must be specified')

    input_group = parser.add_argument_group('Input files')
    input_group.add_argument('--presence_fa', help='FASTA file of genes whose presence you want to check for', metavar='Filename')

    input_group.add_argument('--variants_fa', help='FASTA file of genes that should only be reported if they have a given variant (variants specified in the tsv file given by --metadata_tsv', metavar='Filename')
    input_group.add_argument('--noncoding_fa', help='FASTA file of generic sequences to look for', metavar='Filename')
    input_group.add_argument('--metadata_tsv', help='tsv file of metadata about the sequences/variants of interest', metavar='Filename')

    other_group = parser.add_argument_group('Other options')
    other_group.add_argument('--genetic_code', type=int, help='Number of genetic code to use. Currently supported 1,4,11 [%(default)s]', choices=[1,4,11], default=11, metavar='INT')
    other_group.add_argument('-m', '--min_gene_length', type=int, help='Minimum length in nucleotides of gene [%(default)s]', metavar='INT', default=6)
    other_group.add_argument('-n', '--max_gene_length', type=int, help='Maximum length in nucleotides of gene [%(default)s]', metavar='INT', default=10000)

    parser.add_argument('outprefix', help='Prefix of names of output files')

    options = parser.parse_args()

    ref_data = ariba.reference_data.ReferenceData(
        presence_absence_fa=options.presence_fa,
        variants_only_fa=options.variants_fa,
        non_coding_fa=options.noncoding_fa,
        metadata_tsv=options.metadata_tsv,
        min_gene_length=options.min_gene_length,
        max_gene_length=options.max_gene_length,
        genetic_code=options.genetic_code,
    )

    ref_data.sanity_check(options.outprefix)
