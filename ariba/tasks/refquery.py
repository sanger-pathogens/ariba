import sys
import argparse
from ariba import refdata_query

def run():
    parser = argparse.ArgumentParser(
        description = 'Find cluster or sequence information from prepareref directory',
        usage = 'ariba refquery <prepareref directory> <cluster|seq> <cluster name|sequence name>',
    )
    parser.add_argument('prepareref_dir', help='Name of directory output by prepareref')
    parser.add_argument('query_type', choices=['cluster', 'seq'], help='Use "cluster" to get the sequences in a cluster, or "seq" to get information about a sequence')
    parser.add_argument('search_name', help='Name of cluster or sequence to search for')
    options = parser.parse_args()

    rquery = refdata_query.RefdataQuery(options.prepareref_dir)
    rquery.query(options.query_type, options.search_name)
