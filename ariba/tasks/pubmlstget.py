import argparse
from ariba import pubmlst_ref_preparer


def run(options):
    preparer = pubmlst_ref_preparer.PubmlstRefPreparer(
        options.species,
        options.outdir,
        verbose=options.verbose
    )
    preparer.run()

