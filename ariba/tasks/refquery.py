import sys
import argparse
from ariba import refdata_query

def run(options):
    rquery = refdata_query.RefdataQuery(options.prepareref_dir)
    rquery.query(options.query_type, options.search_name)
