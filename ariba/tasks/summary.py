import argparse
import ariba

def run():
    parser = argparse.ArgumentParser(
        description = 'Make a summry of ARIBA report files',
        usage = 'ariba summary [options] <outfile> [infiles]',
        epilog = 'Files must be listed after the output file and/or the option --fofn must be used. If both used, all files in the filename specified by --fofn AND the files listed after the output file will be used as input')
    parser.add_argument('-f', '--fofn', help='File of filenames of ariba reports to be summarised. Must b used if no input files listed after the outfile', metavar='FILENAME')
    parser.add_argument('--min_id', type=float, help='Minimum percent identity cutoff to count as assembled [%(default)s]', default=90, metavar='FLOAT')
    parser.add_argument('outfile', help='Name of output file. If file ends with ".xls", then an excel spreadsheet is written. Otherwise a tsv file is written')
    parser.add_argument('infiles', nargs='*', help='Files to be summarised')
    options = parser.parse_args()
    if len(options.infiles) == 0:
        options.infiles = None

    s = ariba.summary.Summary(
        options.outfile,
        fofn=options.fofn,
        filenames=options.infiles,
        min_id=options.min_id
    ) 
    s.run()
