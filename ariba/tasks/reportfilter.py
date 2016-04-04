import argparse
import ariba

def run():
    parser = argparse.ArgumentParser(
        description = 'Filters an ARIBA report tsv file',
        usage = 'ariba reportfilter [options] <infile> <outprefix>'
    )
    parser.add_argument('--min_pc_id', type=float, help='Minimum percent identity of nucmer match between contig and reference [%(default)s]', default=90.0, metavar='FLOAT')
    parser.add_argument('--min_ref_base_asm', type=int, help='Minimum number of reference bases matching assembly [%(default)s]', default=1, metavar='INT')
    parser.add_argument('--keep_without_known_var', action='store_true', help='Use this option to not filter out where there is a known variant, but the assembly has the wild type. By default these rows are removed.')
    parser.add_argument('infile', help='Name of input tsv file')
    parser.add_argument('outprefix', help='Prefix of output files. outprefix.tsv and outprefix.xls will be made')
    options = parser.parse_args()

    rf = ariba.report_filter.ReportFilter(
        infile=options.infile,
        min_pc_ident=options.min_pc_id,
        min_ref_base_assembled=options.min_ref_base_asm,
        ignore_not_has_known_variant=not options.keep_without_known_var,
    )
    rf.run(options.outprefix)

