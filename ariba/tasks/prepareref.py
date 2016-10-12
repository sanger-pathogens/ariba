import sys
import argparse
from ariba import ref_preparer, external_progs, versions

def run(options):
    if options.no_cdhit and options.cdhit_clusters is not None:
        sys.exit('Cannot use both --no_cdhit and --cdhit_clusters. Neither or exactly one of those options must be used')

    extern_progs, version_report_lines = versions.get_all_versions()
    if options.verbose:
        print(*version_report_lines, sep='\n')

    preparer = ref_preparer.RefPreparer(
        options.fasta_files,
        extern_progs,
        metadata_tsv_files=options.tsv_files,
        all_coding=options.all_coding,
        version_report_lines=version_report_lines,
        min_gene_length=options.min_gene_length,
        max_gene_length=options.max_gene_length,
        genetic_code=options.genetic_code,
        cdhit_min_id=options.cdhit_min_id,
        cdhit_min_length=options.cdhit_min_length,
        run_cdhit=not options.no_cdhit,
        clusters_file=options.cdhit_clusters,
        threads=options.threads,
        verbose=options.verbose,
        force=options.force,
    )

    preparer.run(options.outdir)
