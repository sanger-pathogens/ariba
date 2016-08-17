import argparse
import ariba


def use_preset(options):
    if options.preset is None:
        return options

    preset_to_vals = {
        'minimal': {
            'cluster_cols': 'match',
            'col_filter': 'y',
            'row_filter': 'y',
        },
        'cluster_small': {
            'cluster_cols': 'assembled,match,ref_seq,known_var',
            'col_filter': 'y',
            'row_filter': 'y',
        },
        'cluster_all': {
            'cluster_cols': 'assembled,match,ref_seq,pct_id,known_var,novel_var',
            'col_filter': 'y',
            'row_filter': 'y',
        },
        'cluster_var_groups': {
            'cluster_cols': 'assembled,match,ref_seq,pct_id,known_var,novel_var',
            'col_filter': 'y',
            'row_filter': 'y',
        },
        'all': {
            'cluster_cols': 'assembled,match,ref_seq,pct_id,known_var,novel_var',
            'col_filter': 'y',
            'row_filter': 'y',
        },
        'all_no_filter': {
            'cluster_cols': 'assembled,match,ref_seq,pct_id,known_var,novel_var',
            'col_filter': 'n',
            'row_filter': 'n',
        },
    }

    assert options.preset in preset_to_vals

    for key, val in preset_to_vals[options.preset].items():
        exec('options.' + key + ' = "' + val + '"')

    if options.preset in {'cluster_var_groups', 'all', 'all_no_filter'}:
        options.v_groups = True

    if options.preset in {'all', 'all_no_filter'}:
        options.known_variants = True
        options.novel_variants = True

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
        cluster_cols=options.cluster_cols,
        make_phandango_tree=(not options.no_tree),
        only_clusters=None if options.only_cluster is None else {options.only_cluster},
        show_var_groups=options.v_groups,
        show_known_vars=options.known_variants,
        show_novel_vars=options.novel_variants,
        verbose=options.verbose
    )
    s.run()
