import argparse
import sys
import ariba

def run(options):
    expander = ariba.report_flag_expander.ReportFlagExpander(options.infile, options.outfile)
    expander.run()

