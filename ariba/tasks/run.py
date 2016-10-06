import argparse
import os
import shutil
import sys
import ariba


def run(options):
    reads_not_found = []

    for filename in [options.reads_1, options.reads_2]:
        if not os.path.exists(filename):
            reads_not_found.append(filename)
        elif options.verbose:
            print('Found reads file:', filename)

    if len(reads_not_found):
        print('\nThe following reads file(s) were not found:', file=sys.stderr)
        print(*reads_not_found, sep='\n', file=sys.stderr)
        print('Cannot continue', file=sys.stderr)
        sys.exit(1)

    if (options.reads_1 == options.reads_2):
        print('Same file provided for forwards and reverse reads. Cannot continue', file=sys.stderr)
        sys.exit(1)

    if not os.path.exists(options.prepareref_dir):
        print('Input directory', options.prepareref_dir, 'not found. Cannot continue', file=sys.stderr)
        sys.exit(1)

    if options.force and os.path.exists(options.outdir):
        shutil.rmtree(options.outdir)

    if os.path.exists(options.outdir):
        print('Output directory already exists. ARIBA makes the output directory. Cannot continue.', file=sys.stderr)
        sys.exit(1)

    extern_progs, version_report_lines = ariba.versions.get_all_versions()
    if options.verbose:
        print(*version_report_lines, sep='\n')

    c = ariba.clusters.Clusters(
          options.prepareref_dir,
          options.reads_1,
          options.reads_2,
          options.outdir,
          extern_progs,
          version_report_lines=version_report_lines,
          assembly_coverage=options.assembly_cov,
          assembler='fermilite',
          threads=options.threads,
          verbose=options.verbose,
          min_scaff_depth=options.min_scaff_depth,
          nucmer_min_id=options.nucmer_min_id,
          nucmer_min_len=options.nucmer_min_len,
          nucmer_breaklen=options.nucmer_breaklen,
          assembled_threshold=options.assembled_threshold,
          unique_threshold=options.unique_threshold,
          max_gene_nt_extend=options.gene_nt_extend,
          clean=(not options.noclean),
          tmp_dir=options.tmp_dir,
        )
    c.run()

