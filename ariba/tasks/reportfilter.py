import argparse
import sys
import ariba

def run():
    parser = argparse.ArgumentParser(
        description = 'Filters an ARIBA report tsv file',
        usage = 'ariba reportfilter [options] <infile> <outprefix>'
    )
    parser.add_argument('--exclude_flags', help='Comma-separated list of flags to exclude. [%(default)s]', default='assembly_fail,ref_seq_choose_fail')
    parser.add_argument('--min_pc_id', type=float, help='Minimum percent identity of nucmer match between contig and reference [%(default)s]', default=90.0, metavar='FLOAT')
    parser.add_argument('--min_ref_base_asm', type=int, help='Minimum number of reference bases matching assembly [%(default)s]', default=1, metavar='INT')
    parser.add_argument('--keep_syn', action='store_true', help='Keep synonymous variants (by default they are removed')
    parser.add_argument('--discard_without_known_var', action='store_true', help='Applies to variant only genes. Filter out where there is a known variant, but the assembly has the wild type. By default these rows are kept.')
    parser.add_argument('infile', help='Name of input tsv file')
    parser.add_argument('outprefix', help='Prefix of output files. outprefix.tsv and outprefix.xls will be made')
    options = parser.parse_args()

    flags_to_exclude = options.exclude_flags.split(',')
    allowed_flags = set(ariba.flag.flags_in_order)
    bad_flags = [x for x in flags_to_exclude if x not in allowed_flags]
    if len(bad_flags):
        print('Error in option --exclude_flags. The following were not recognised:', ','.join(bad_flags), file=sys.stderr)
        print('Must choose from:', ','.join(ariba.flag.flags_in_order), file=sys.stderr)
        sys.exit(1)

    rf = ariba.report_filter.ReportFilter(
        infile=options.infile,
        min_pc_ident=options.min_pc_id,
        min_ref_base_assembled=options.min_ref_base_asm,
        ignore_not_has_known_variant=options.discard_without_known_var,
        remove_synonymous_snps=not options.keep_syn,
    )
    rf.run(options.outprefix)

