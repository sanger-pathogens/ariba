import argparse
from ariba import aln_to_metadata


def run():
    coding_choices = ['coding', 'noncoding']
    parser = argparse.ArgumentParser(
        description = 'Converts multi-alignment fasta and SNP info to metadata',
        usage = 'ariba aln2meta [options] <aln_fasta> <variants_tsv> <(non)coding> <cluster_rep> <outprefix>'
    )

    parser.add_argument('--genetic_code', type=int, help='Number of genetic code to use. Currently supported 1,4,11 [%(default)s]', choices=[1,4,11], default=11, metavar='INT')
    parser.add_argument('aln_fasta', help='Multi-fasta file of alignments')
    parser.add_argument('variants_tsv', help='TSV file of variants information')
    parser.add_argument('coding_or_non', help='Sequences are coding or noncoding. Must be one of: ' + ' '.join(coding_choices), choices=coding_choices, metavar='(non)coding')
    parser.add_argument('cluster_rep', help='Name of sequence to be used as cluster representative. Must exactly match a sequence in aln_fasta file')
    parser.add_argument('outprefix', help='Prefix of output filenames')
    options = parser.parse_args()

    aln_to_meta = aln_to_metadata.AlnToMetadata(
      options.aln_fasta,
      options.variants_tsv,
      options.coding_or_non == 'coding',
      options.cluster_rep,
      genetic_code=options.genetic_code
    )
    aln_to_meta.run(options.outprefix)

