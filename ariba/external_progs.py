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


prog_to_version_cmd = {
    'bcftools': ('', re.compile('^Version: ([0-9\.]+)')),
    'nucmer': ('--version', re.compile('^NUCmer \(NUCleotide MUMmer\) version ([0-9\.]+)')),
    'smalt': ('version', re.compile('^Version: ([0-9\.]+)')),
    'spades': ('', re.compile('^SPAdes genome assembler v.([0-9\.]+)')),
    'samtools': ('', re.compile('^Version: ([0-9\.]+)')),
    'sspace': ('', re.compile('^Usage: .*pl \[SSPACE_(.*)\]')),
    'gapfiller': ('', re.compile('^Usage: .*pl \[GapFiller_(.*)\]')),
    'velvetg': ('', re.compile('Version ([0-9\.]+)')),
    'velveth': ('', re.compile('Version ([0-9\.]+)')),
}


def get_version(prog, path=None):
    assert prog in prog_to_version_cmd
    if path is None:
        path = prog

    if not is_in_path(path):
        raise Error('Error getting version of ' + path + ' - not found in path.')

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


def check_versions(opts, verbose=False):
    if verbose:
        print('{:_^79}'.format(' Checking dependencies and their versions '))
        print('tool', 'version', 'path', sep='\t')

    l = [
        ('bcftools', opts.bcftools),
        ('nucmer', 'nucmer'),
        ('smalt', opts.smalt),
        ('samtools', opts.samtools),
        ('sspace', opts.sspace),
        ('gapfiller', opts.gapfiller),
    ]
    
    if opts.assembler == 'spades':
        l.append(('spades', opts.spades))
    elif opts.assembler == 'velvet':
        l.append(('velvetg', opts.velvet + 'g'))
        l.append(('velveth', opts.velvet + 'h'))
    else:
        raise Error('Assembler ' + opts.assembler + ' not recognised. Cannot continue')

    min_versions = {
        'bcftools': '1.2',
        'samtools': '1.2',
        'nucmer': '3.1',
        'spades': '3.5.0',
        'smalt': '0.7.4',
        'velvetg': '1.2.07',
        'velveth': '1.2.07',
    }

    errors = []

    for t in l:
        version, path = get_version(t[0], path=t[1])
        if verbose:
            print(t[0], version, path, sep='\t')

        if t[0] in min_versions and LooseVersion(version) < LooseVersion(min_versions[t[0]]):
            errors.append(' '.join(['Found version', version, 'of', t[0], 'which is too low! Please update to at least', min_versions[t[0]] + '\n   Found it here:', path]))

    if len(errors):
        for e in errors:
            print('\n*** Error! Bad dependency! ***', file=sys.stderr)
            print(e, file=sys.stderr)
            print()
        raise Error('Cannot continue. Some dependencies need updating')

    if verbose:
        print('\n... dependencies look OK.\n')
