import argparse
from ariba import ref_genes_getter


def run(options):
    getter = ref_genes_getter.RefGenesGetter(
        options.db,
        genetic_code=options.genetic_code,
        version=options.version
    )
    getter.run(options.outprefix)

