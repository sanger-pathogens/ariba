import argparse
import pyfastaq
import ariba


def run():
    parser = argparse.ArgumentParser(
        description = 'ARIBA: Antibiotic Resistance Identification By Assembly',
        usage = 'ariba run [options] <db.fa> <reads1.fq> <reads2.fq> <outdir>')
    parser.add_argument('db_fasta', help='FASTA file of reference genes')
    parser.add_argument('reads_1', help='Name of fwd reads fastq file')
    parser.add_argument('reads_2', help='Name of rev reads fastq file')
    parser.add_argument('outdir', help='Output directory (must not already exist)')

    smalt_group = parser.add_argument_group('smalt options')
    smalt_group.add_argument('--smalt_k', type=int, help='kmer to use when indexing with smalt (smalt index -k) [%(default)s]', default=13, metavar='INT')
    smalt_group.add_argument('--smalt_s', type=int, help='Step length to use when indexing with smalt (see smalt index -s) [%(default)s]', default=2, metavar='INT')
    smalt_group.add_argument('--smalt_min_id', type=float, help='Minimum identity to report a match (smalt map -y) [%(default)s]', default=0.9, metavar='FLOAT')

    nucmer_group = parser.add_argument_group('nucmer options')
    nucmer_group.add_argument('--nucmer_min_id', type=int, help='Minimum alignment identity (delta-filter -i) [%(default)s]', default=90, metavar='INT')
    nucmer_group.add_argument('--nucmer_min_len', type=int, help='Minimum alignment length (delta-filter -i) [%(default)s]', default=50, metavar='INT')
    nucmer_group.add_argument('--nucmer_breaklen', type=int, help='Value to use for -breaklen when running nucmer [%(default)s]', default=50, metavar='INT')

    executables_group = parser.add_argument_group('executables locations')
    executables_group.add_argument('--bcftools', help='Path to bcftools executable [%(default)s]', default='bcftools', metavar='PATH')
    executables_group.add_argument('--gapfiller', help='Path to GapFiller.pl [%(default)s]', default='GapFiller.pl', metavar='PATH')
    executables_group.add_argument('--samtools', help='Path to samtools executable [%(default)s]', default='samtools', metavar='PATH')
    executables_group.add_argument('--sspace', help='Path to SSPACE_Basic_v2.0.pl [%(default)s]', default='SSPACE_Basic_v2.0.pl', metavar='PATH')
    executables_group.add_argument('--velvet', help='Path to prefix of velvet{g,h} [%(default)s]', default='velvet', metavar='PATH')

    other_group = parser.add_argument_group('Other options')
    other_group.add_argument('--genetic_code', type=int, help='Number of genetic code to use. Currently supported 1 (default) or 4', choices=[1,4], default=1, metavar='INT')
    other_group.add_argument('--min_scaff_depth', type=int, help='Minimum number of read pairs needed as evidence for scaffold link between two contigs. This is also the value used for sspace -k when scaffolding [%(default)s]', default=10, metavar='INT')
    allowed_assemblers = ['velvet', 'spades']
    other_group.add_argument('--assembler', help='Assembler to use. Although SPAdes can be used, it is experimental and velvet is currently recommended. Available options: ' + ','.join(allowed_assemblers) + ' [%(default)s]', choices=allowed_assemblers, default='velvet', metavar='Assembler')
    other_group.add_argument('--assembler_k', help='kmer size to use with assembler. Default is 2/3 of the read length', metavar='INT')
    other_group.add_argument('--spades_only_assembler', action='store_true', help='Use --only-assembler when running SPAdes')
    other_group.add_argument('--threads', type=int, help='Number of threads for smalt and spades [%(default)s]', default=1, metavar='INT')
    other_group.add_argument('--assembled_threshold', type=float, help='If proportion of gene assembled (regardless of into how many contigs) is at least this value then the flag gene_assembled is set [%(default)s]', default=0.95, metavar='FLOAT (between 0 and 1)')
    other_group.add_argument('--unique_threshold', type=float, help='If proportion of bases in gene assembled more than once is <= this value, then the flag unique_contig is set [%(default)s]', default=0.03, metavar='FLOAT (between 0 and 1)')
    other_group.add_argument('--verbose', action='store_true', help='Be verbose')
    options = parser.parse_args()

    pyfastaq.sequences.codon2aa = pyfastaq.genetic_codes.codes[options.genetic_code]

    c = ariba.clusters.Clusters(
          options.db_fasta,
          options.reads_1,
          options.reads_2,
          options.outdir,
          assembly_kmer=None,
          assembler=options.assembler,
          threads=options.threads,
          verbose=options.verbose,
          smalt_k=options.smalt_k,
          smalt_s=options.smalt_s,
          smalt_min_id=options.smalt_min_id,
          min_scaff_depth=options.min_scaff_depth,
          nucmer_min_id=options.nucmer_min_id,
          nucmer_min_len=options.nucmer_min_len,
          nucmer_breaklen=options.nucmer_breaklen,
          spades_only_assembler=options.spades_only_assembler,
          assembled_threshold=options.assembled_threshold,
          unique_threshold=options.unique_threshold,
          bcftools_exe=options.bcftools,
          gapfiller_exe=options.gapfiller,
          samtools_exe=options.samtools,
          sspace_exe=options.sspace,
          velvet_exe=options.velvet,
        )
    c.run()

