import argparse
import ariba


def use_preset(options):
    if options.preset is None:
        return options

    preset_to_vals = {
        'minimal': {
            'cluster_cols': 'has_res',
            'col_filter': 'yes',
            'row_filter': 'yes',
            'known_vars': 'no',
            'novel_vars': 'no'
        },
        'cluster_small': {
            'cluster_cols': 'assembled,has_res,ref_seq,known_var',
            'col_filter': 'yes',
            'row_filter': 'yes',
            'known_vars': 'no',
            'novel_vars': 'no'
        },
        'cluster_all': {
            'cluster_cols': 'assembled,has_res,ref_seq,pct_id,known_var,novel_var',
            'col_filter': 'yes',
            'row_filter': 'yes',
            'known_vars': 'no',
            'novel_vars': 'no'
        },
        'cluster_known_vars': {
            'cluster_cols': 'assembled,has_res,ref_seq,pct_id,known_var,novel_var',
            'col_filter': 'yes',
            'row_filter': 'yes',
            'known_vars': 'yes',
            'novel_vars': 'no'
        },
        'all': {
            'cluster_cols': 'assembled,has_res,ref_seq,pct_id,known_var,novel_var',
            'col_filter': 'yes',
            'row_filter': 'yes',
            'known_vars': 'yes',
            'novel_vars': 'yes'
        },
        'all_no_filter': {
            'cluster_cols': 'assembled,has_res,ref_seq,pct_id,known_var,novel_var',
            'col_filter': 'no',
            'row_filter': 'no',
            'known_vars': 'yes',
            'novel_vars': 'yes'
        },
    }

    assert options.preset in preset_to_vals

    for key, val in preset_to_vals[options.preset].items():
        exec('options.' + key + ' = "' + val + '"')

    return options


def run():
    presets = ['minimal', 'cluster_small', 'cluster_all', 'cluster_known_vars', 'all', 'all_no_filter']

    parser = argparse.ArgumentParser(
        description = 'Make a summary of ARIBA report files, and Phandango files',
        usage = 'ariba summary [options] <outprefix> [report1.tsv report2.tsv ...]',
        epilog = 'Files must be listed after the output file and/or the option --fofn must be used. If both used, all files in the filename specified by --fofn AND the files listed after the output file will be used as input.')
    parser.add_argument('-f', '--fofn', help='File of filenames of ariba reports in tsv format (not xls) to be summarised. Must be used if no input files listed after the outfile.', metavar='FILENAME')
    parser.add_argument('--preset', choices=presets, help='Shorthand for setting --cluster_cols,--col_filter,--row_filter,--known_vars,--novel_vars. Using this overrides those options', metavar='|'.join(presets))
    parser.add_argument('--cluster_cols', help='Comma separated list of cluster columns to include. Choose from: assembled, has_res, ref_seq, pct_id, known_var, novel_var [%(default)s]', default='has_res', metavar='col1,col2,...')
    parser.add_argument('--col_filter', choices=['y', 'n'], default='y', help='Choose whether columns where all values are "no" or "NA" are removed [%(default)s]', metavar='y|n')
    parser.add_argument('--row_filter', choices=['y', 'n'], default='y', help='Choose whether rows where all values are "no" or "NA" are removed [%(default)s]', metavar='y|n')
    parser.add_argument('--known_vars', choices=['y', 'n'], default='n', help='Output a column for every known variant [%(default)s]', metavar='y|n')
    parser.add_argument('--novel_vars', choices=['y', 'n'], default='n', help='Output a column for every novel variant [%(default)s]', metavar='y|n')
    parser.add_argument('--min_id', type=float, help='Minimum percent identity cutoff to count as assembled [%(default)s]', default=90, metavar='FLOAT')
    parser.add_argument('--verbose', action='store_true', help='Be verbose')
    parser.add_argument('outprefix', help='Prefix of output files')
    parser.add_argument('infiles', nargs='*', help='Files to be summarised')
    options = parser.parse_args()
    if len(options.infiles) == 0:
        options.infiles = None

    options = use_preset(options)

    s = ariba.summary.Summary(
        options.outprefix,
        fofn=options.fofn,
        filenames=options.infiles,
        include_all_known_variant_columns=options.known_vars == 'y',
        include_all_novel_variant_columns=options.novel_vars == 'y',
        filter_rows=options.col_filter == 'y',
        filter_columns=options.row_filter == 'y',
        min_id=options.min_id,
        cluster_cols=options.cluster_cols,
        verbose=options.verbose
    )
    s.run()
