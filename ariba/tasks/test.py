import argparse
import subprocess
import shutil
import os
import sys
import ariba


def boxymcboxface(message):
    print('-' * 79)
    print('|', '=' * 77, '|', sep='')
    print('|', '{: ^75}'.format(message), '|')
    print('|', '=' * 77, '|', sep='')
    print('-' * 79)


def run(options):
    ariba_exe = os.path.abspath(sys.argv[0])

    print('Running ARIBA on test data...')

    boxymcboxface('Preparing input data')

    try:
        os.mkdir(options.outdir)
        os.chdir(options.outdir)
    except:
        print('Error making output directory "', options.outdir, '". Cannot continue.', sep='', file=sys.stderr)
        sys.exit(1)

    print('Made output directory ', options.outdir, '. Copying test data files into it:', sep='')

    modules_dir = os.path.dirname(os.path.abspath(ariba.__file__))
    test_data_dir = os.path.join(modules_dir, 'test_run_data')

    for filename in ['ref_seqs.fa', 'metadata.tsv', 'reads_1.fq', 'reads_2.fq']:
        shutil.copy(os.path.join(test_data_dir, filename), filename)
        print('    copied', filename)


    boxymcboxface('Try running ariba prepareref')

    prepareref_command = ' '.join([
        ariba_exe,
        'prepareref',
        '--verbose',
        '-f ref_seqs.fa',
        '-m metadata.tsv',
        '--threads', str(options.threads),
        'PREPAREREF',
    ])

    print('\nRunning ariba prepareref with:', prepareref_command, '', sep='\n')
    return_code = subprocess.call(prepareref_command, shell=True)

    if return_code != 0:
        print('\nSomething went wrong. See above for error message(s). Return code was', return_code)
        sys.exit(1)

    print()
    print('ariba prepareref finished OK')


    ariba_command = ' '.join([
        ariba_exe,
        'run',
        '--verbose',
        '--threads', str(options.threads),
        'PREPAREREF',
        'reads_1.fq',
        'reads_2.fq',
        'OUT'
    ])

    boxymcboxface('Try running ariba run')
    print('\nRunning ARIBA with:', ariba_command, '', sep='\n')

    return_code = subprocess.call(ariba_command, shell=True)

    if return_code != 0:
        print('\nSomething went wrong. See above for error message(s). Return code was', return_code)
        sys.exit(1)

    print()
    print('ariba run finished OK')
    print('Finished run on test data OK')
