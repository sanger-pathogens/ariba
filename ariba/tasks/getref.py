import argparse
from ariba import ref_genes_getter


def run():
    allowed_dbs = ['argannot', 'card', 'resfinder','vfdb_core','vfdb_full']
    parser = argparse.ArgumentParser(
        description = 'Downloads reference data',
        usage = 'ariba getref [options] <' + '|'.join(allowed_dbs) + '> <outprefix>'
    )

    parser.add_argument('--genetic_code', type=int, help='Number of genetic code to use. Currently supported 1,4,11 [%(default)s]', choices=[1,4,11], default=11, metavar='INT')
    parser.add_argument('db', help='Database to download. Must be one of: ' + ' '.join(allowed_dbs), choices=allowed_dbs)
    parser.add_argument('outprefix', help='Prefix of output filenames')
    options = parser.parse_args()

    getter = ref_genes_getter.RefGenesGetter(options.db, genetic_code=options.genetic_code)
    getter.run(options.outprefix)

