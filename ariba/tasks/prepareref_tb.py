import sys
import argparse
from ariba import tb

def run(options):
    tb.make_prepareref_dir(options.outdir)

