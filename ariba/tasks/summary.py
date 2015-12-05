import argparse
import ariba

def run():
    parser = argparse.ArgumentParser(
        description = 'Make a summry of ARIBA report files',
        usage = 'ariba summary [options] <outfile> [report1.tsv report2.tsv ...]',
        epilog = 'Files must be listed after the output file and/or the option --fofn must be used. If both used, all files in the filename specified by --fofn AND the files listed after the output file will be used as input. The input report files must be in tsv format, not xls.')
    parser.add_argument('-f', '--fofn', help='File of filenames of ariba reports in tsv format (not xls) to be summarised. Must be used if no input files listed after the outfile.', metavar='FILENAME')
    parser.add_argument('--js_candy', help='Write files that can be used as input to JS Candy, named with the given prefix. Writes a tree and csv file, which can be drag and dropped straight into JS Candy', metavar='files_prefix')
    parser.add_argument('--min_id', type=float, help='Minimum percent identity cutoff to count as assembled [%(default)s]', default=90, metavar='FLOAT')
    parser.add_argument('--no_filter', action='store_true', help='Do not filter rows or columns of output that are all 0 (by default, they are removed from the output)')
    parser.add_argument('outfile', help='Name of output file. If file ends with ".xls", then an excel spreadsheet is written. Otherwise a tsv file is written')
    parser.add_argument('infiles', nargs='*', help='Files to be summarised')
    options = parser.parse_args()
    if len(options.infiles) == 0:
        options.infiles = None

    s = ariba.summary.Summary(
        options.outfile,
        fofn=options.fofn,
        filenames=options.infiles,
        filter_output=(not options.no_filter),
        min_id=options.min_id,
        js_candy_prefix=options.js_candy
    )
    s.run()
