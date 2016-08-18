import argparse
from ariba import aln_to_metadata


def run(options):
    aln_to_meta = aln_to_metadata.AlnToMetadata(
      options.aln_fasta,
      options.variants_tsv,
      options.coding_or_non == 'coding',
      options.variant_only,
      genetic_code=options.genetic_code
    )
    aln_to_meta.run(options.outprefix)

