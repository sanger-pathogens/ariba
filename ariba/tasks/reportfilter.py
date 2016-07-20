import argparse
import sys
import ariba

def run(options):
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
    rf.run(options.outfile)

