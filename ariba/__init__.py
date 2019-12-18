from pkg_resources import get_distribution

try:
    __version__ = get_distribution('ariba').version
except:
    __version__ = 'local'


__all__ = [
    'aln_to_metadata',
    'assembly',
    'assembly_compare',
    'assembly_variants',
    'bam_parse',
    'card_record',
    'cdhit',
    'cluster',
    'clusters',
    'common',
    'external_progs',
    'faidx',
    'flag',
    'histogram',
    'link',
    'mapping',
    'megares_data_finder',
    'megares_zip_parser',
    'mic_plotter',
    'mlst_profile',
    'mlst_reporter',
    'pubmlst_getter',
    'pubmlst_ref_preparer',
    'read_filter',
    'read_store',
    'refdata_query',
    'reference_data',
    'ref_genes_getter',
    'ref_preparer',
    'ref_seq_chooser',
    'report',
    'report_filter',
    'report_flag_expander',
    'scaffold_graph',
    'samtools_variants',
    'sequence_metadata',
    'sequence_variant',
    'summary',
    'summary_cluster',
    'summary_cluster_variant',
    'summary_sample',
    'tasks',
    'tb',
    'versions',
    'vfdb_parser',
]

from ariba import *
