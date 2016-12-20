import argparse
from ariba import pubmlst_getter


def run(options):
    getter = pubmlst_getter.PubmlstGetter()
    getter.print_available_species()

