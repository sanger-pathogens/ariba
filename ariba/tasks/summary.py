import argparse
import ariba


def use_preset(options):
    if options.preset is None:
        return options

    preset_to_vals = {
        'minimal': {
            'cluster_cols': 'match',
            'variant_cols': '',
            'col_filter': 'y',
            'row_filter': 'y',
            'var_groups': 'n',
            'known_vars': 'n',
            'novel_vars': 'n'
        },
        'cluster_small': {
            'cluster_cols': 'assembled,match,ref_seq,known_var',
            'variant_cols': '',
            'col_filter': 'y',
            'row_filter': 'y',
            'var_groups': 'n',
            'known_vars': 'n',
            'novel_vars': 'n'
        },
        'cluster_all': {
            'cluster_cols': 'assembled,match,ref_seq,pct_id,known_var,novel_var',
            'variant_cols': '',
            'col_filter': 'y',
            'row_filter': 'y',
            'var_groups': 'n',
            'known_vars': 'n',
            'novel_vars': 'n'
        },
        'cluster_var_groups': {
            'cluster_cols': 'assembled,match,ref_seq,pct_id,known_var,novel_var',
            'variant_cols': 'groups',
            'col_filter': 'y',
            'row_filter': 'y',
            'var_groups': 'y',
            'known_vars': 'n',
            'novel_vars': 'n'
        },
        'cluster_known_vars': {
            'cluster_cols': 'assembled,match,ref_seq,pct_id,known_var,novel_var',
            'variant_cols': 'groups,grouped,ungrouped',
            'col_filter': 'y',
            'row_filter': 'y',
            'var_groups': 'y',
            'known_vars': 'y',
            'novel_vars': 'n'
        },
        'all': {
            'cluster_cols': 'assembled,match,ref_seq,pct_id,known_var,novel_var',
            'variant_cols': 'groups,grouped,ungrouped,novel',
            'col_filter': 'y',
            'row_filter': 'y',
            'var_groups': 'y',
            'known_vars': 'y',
            'novel_vars': 'y'
        },
        'all_no_filter': {
            'cluster_cols': 'assembled,match,ref_seq,pct_id,known_var,novel_var',
            'variant_cols': 'groups,grouped,ungrouped,novel',
            'col_filter': 'n',
            'row_filter': 'n',
            'var_groups': 'y',
            'known_vars': 'y',
            'novel_vars': 'y'
        },
    }

    assert options.preset in preset_to_vals

    for key, val in preset_to_vals[options.preset].items():
        exec('options.' + key + ' = "' + val + '"')

    return options


def run(options):
    if len(options.infiles) == 0:
        options.infiles = None

    options = use_preset(options)

    s = ariba.summary.Summary(
        options.outprefix,
        fofn=options.fofn,
        filenames=options.infiles,
        filter_rows=options.col_filter == 'y',
        filter_columns=options.row_filter == 'y',
        min_id=options.min_id,
        show_known_het=options.het,
        cluster_cols=options.cluster_cols,
        variant_cols=options.var_cols,
        verbose=options.verbose
    )
    s.run()
