import shutil
import subprocess
import os
from distutils.version import LooseVersion
import re
import sys
import pyfastaq
from ariba import common

class Error (Exception): pass

def is_in_path(prog):
    return shutil.which(prog) is not None


prog_to_default = {
    'bcftools': 'bcftools',
    'bowtie2': 'bowtie2',
    'cdhit': 'cd-hit-est',
    'gapfiller': 'GapFiller.pl',
    'nucmer' : 'nucmer',
    'samtools': 'samtools',
    'smalt': 'smalt',
    'spades': 'spades.py',
    'sspace': 'SSPACE_Basic_v2.0.pl',
    'velvetg': 'velvetg',
    'velveth': 'velveth',
}


prog_to_env_var = {
    'bcftools': 'ARIBA_BCFTOOLS',
    'samtools': 'ARIBA_SAMTOOLS',
    'spades': 'ARIBA_SPADES', 
}


prog_to_version_cmd = {
    'bcftools': ('', re.compile('^Version: ([0-9\.]+)')),
    'bowtie2': ('--version', re.compile('.*bowtie2-align version (.*)$')),
    'cdhit': ('', re.compile('CD-HIT version ([0-9\.]+) \(')),
    'gapfiller': ('', re.compile('^Usage: .*pl \[GapFiller_(.*)\]')),
    'nucmer': ('--version', re.compile('^NUCmer \(NUCleotide MUMmer\) version ([0-9\.]+)')),
    'samtools': ('', re.compile('^Version: ([0-9\.]+)')),
    'smalt': ('version', re.compile('^Version: ([0-9\.]+)')),
    'spades': ('', re.compile('^SPAdes genome assembler v.([0-9\.]+)')),
    'sspace': ('', re.compile('^Usage: .*pl \[SSPACE_(.*)\]')),
    'velvetg': ('', re.compile('Version ([0-9\.]+)')),
    'velveth': ('', re.compile('Version ([0-9\.]+)')),
}


min_versions = {
    'bcftools': '1.2',
    'bowtie2': '2.1.0',
    'cd-hit': '4.6',
    'nucmer': '3.1',
    'samtools': '1.2',
    'smalt': '0.7.4',
    'spades': '3.5.0',
    'velvetg': '1.2.07',
    'velveth': '1.2.07',
}


def set_path(prog, opts):
    path_from_opts = eval('opts.' + prog)
    if path_from_opts is not None:
        return

    if prog in prog_to_env_var:
        env_var = prog_to_env_var[prog]
        if env_var in os.environ:
            exec('opts.' + prog + ' = "' + os.environ[env_var] + '"')
            return

    exec('opts.' + prog + ' = "' + prog_to_default[prog] + '"')


def get_version(prog, path=None, raise_error=True):
    assert prog in prog_to_version_cmd
    if path is None:
        path = prog

    if not is_in_path(path):
        if raise_error:
            raise Error('Error getting version of ' + path + ' - not found in path.')
        else:
            return 'Not_in_path', 'Not_in_path'

    path = shutil.which(path)

    if prog in ['sspace', 'gapfiller']:
        cmd = 'perl ' + os.path.realpath(shutil.which(path))
        regex = prog_to_version_cmd[prog][1]
    else:
        cmd, regex = prog_to_version_cmd[prog]
        cmd = path + ' ' + cmd

    cmd_output = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()
    cmd_output = common.decode(cmd_output[0]).split('\n')[:-1] + common.decode(cmd_output[1]).split('\n')[:-1]
    for line in cmd_output:
        hits = regex.search(line)
        if hits:
            return hits.group(1), path
    return 'UNKNOWN ...\n I tried running this to get the version: "' + cmd + '"\n and the output didn\'t match this regular expression: "' + regex.pattern + '"', path


def check_versions(opts, verbose=False, not_required=None):
    if not_required is None:
        not_required = set()

    if verbose:
        print('{:_^79}'.format(' Checking dependencies and their versions '))
        print('tool', 'version', 'path', sep='\t')

    to_check = [
        'bcftools',
        'bowtie2',
        'cdhit',
        'nucmer',
        'samtools',
        'sspace',
        'gapfiller',
    ]
    
    if opts.assembler == 'spades':
        to_check.append('spades')
    elif opts.assembler == 'velvet':
        to_check.append('velvetg')
        to_check.append('velveth')
    else:
        raise Error('Assembler ' + opts.assembler + ' not recognised. Cannot continue')

    errors = []
    failed_to_find = set()

    for prog in to_check:
        set_path(prog, opts)
        version, path = get_version(prog, path=eval('opts.' + prog), raise_error=prog not in not_required)
        if verbose:
            print(prog, version, path, sep='\t')
        if path == 'Not_in_path':
            print('\nWARNING:', prog, 'not found in path, so will be skipped during assembly\n', file=sys.stderr)

        if prog in min_versions and LooseVersion(version) < LooseVersion(min_versions[prog]):
            errors.append(' '.join(['Found version', version, 'of', prog, 'which is too low! Please update to at least', min_versions[prog] + '\n   Found it here:', path]))
            failed_to_find.add(prog)

    if len(errors):
        for e in errors:
            print('\n*** Error! Bad dependency! ***', file=sys.stderr)
            print(e, file=sys.stderr)
            print()
        if len(failed_to_find.difference(not_required)) > 0:
            raise Error('Cannot continue. Some dependencies need updating')
        else:
            assert failed_to_find.issubset(not_required)
            if 'sspace' in failed_to_find:
                print('WARNING: SSPACE not found. Will not run scaffolding or gap filling', file=sys.stderr)
            elif 'gapfiller' in failed_to_find:
                print('WARNING: GapFiller not found. Will not run gap filling after scaffolding', file=sys.stderr)

    if verbose:
        print('\nDependencies look OK (but check in case there are warnings about SSPACE or GapFiller)\n')
