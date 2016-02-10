import argparse
import subprocess
import shutil
import os
import sys
import ariba


def run():
    parser = argparse.ArgumentParser(
        description = 'Run ARIBA on a small test dataset',
        usage = 'ariba test [options] <outdir>')
    parser.add_argument('--threads', type=int, help='Number of threads [%(default)s]', default=1, metavar='INT')
    parser.add_argument('outdir', help='Name of output directory')
    options = parser.parse_args()

    print('Running ARIBA on test data...')

    try:
        os.mkdir(options.outdir)
        os.chdir(options.outdir)
    except:
        print('Error making output directory "', options.outdir, '". Cannot continue.', sep='', file=sys.stderr)
        sys.exit(1)

    print('Made output directory. Copying test data files into it:')

    modules_dir = os.path.dirname(os.path.abspath(ariba.__file__))
    test_data_dir = os.path.join(modules_dir, 'test_run_data')

    for filename in ['presence_absence.fa', 'non_coding.fa', 'variants_only.fa', 'metadata.tsv', 'reads_1.fq', 'reads_2.fq']:
        shutil.copy(os.path.join(test_data_dir, filename), filename)
        print('    copied', filename)

    ariba_command = ' '.join([
        sys.argv[0],
        'run',
        '--verbose',
        '--presabs presence_absence.fa',
        '--varonly variants_only.fa',
        '--noncoding non_coding.fa',
        '--metadata metadata.tsv',
        '--threads', str(options.threads),
        'reads_1.fq',
        'reads_2.fq',
        'OUT'
    ])

    print('\nRunning ARIBA with:', ariba_command, '', sep='\n')

    return_code = subprocess.call(ariba_command, shell=True)

    if return_code != 0:
        print('\nSomething went wrong. See above for error message(s). Return code was', return_code)
        sys.exit(1)

    print('-' * 79)
    print('Finished run on test data OK')
