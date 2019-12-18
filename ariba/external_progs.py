import shutil
import subprocess
import os
from distutils.version import LooseVersion
import re
import sys
from ariba import common

class Error (Exception): pass


prog_to_default = {
    'bowtie2': 'bowtie2',
    'cdhit': 'cd-hit-est',
    'nucmer' : 'nucmer',
    'spades' : 'spades.py'
}


prog_to_env_var = {x: 'ARIBA_' + x.upper() for x in prog_to_default if x not in {'nucmer'}}


# Nucmer 3.1 'nucmer --version' outputs this:
# nucmer
# NUCmer (NUCleotide MUMmer) version 3.1
#
# Numcer 4 'nucmer --version' outputs this:
# 4.0.0beta2
#
# ... make the regex permissive and hope things
# still work for later versions
prog_to_version_cmd = {
    'bowtie2': ('--version', re.compile('.*bowtie2.*version (.*)$')),
    'cdhit': ('', re.compile('CD-HIT version ([0-9\.]+) \(')),
    'nucmer': ('--version', re.compile('([0-9]+\.[0-9\.]+.*$)')),
    'spades': ('--version', re.compile('SPAdes\s+v([0-9\.]+)'))
}


min_versions = {
    'bowtie2': '2.1.0',
    'cdhit': '4.6',
    'nucmer': '3.1',
    'spades': '3.11.0'
}

prog_optional = set([
    'spades'
])

class ExternalProgs:
    def __init__(self, verbose=False, fail_on_error=True, using_spades=False):
        self.progs = {}
        self.version_report = []
        self.all_deps_ok = True
        self.versions = {}
        self.using_spades = using_spades

        if verbose:
            print('{:_^79}'.format(' Checking dependencies and their versions '))

        errors = []
        warnings = []

        for prog in sorted(prog_to_default):
            if prog == 'spades' and not self.using_spades:
                continue

            msg_sink = errors
            if prog in prog_optional:
                msg_sink = warnings

            prog_exe = self._get_exe(prog)
            self.progs[prog] = shutil.which(prog_exe)

            if self.progs[prog] is None:
                msg_sink.append(prog + ' not found in path. Looked for ' + prog_exe)

                self.version_report.append('\t'.join([prog, 'NA', 'NOT_FOUND']))
                if verbose:
                    print(self.version_report[-1])
                continue

            got_version, version = self._get_version(prog, self.progs[prog])

            if got_version:
                self.versions[prog] = version
                if prog in min_versions and LooseVersion(version) < LooseVersion(min_versions[prog]):
                    msg_sink.append(' '.join(['Found version', version, 'of', prog, 'which is too low! Please update to at least', min_versions[prog] + '. Found it here:', prog_exe]))
            else:
                self.versions[prog] = None
                msg_sink.append(version)
                version = 'ERROR'

            self.version_report.append('\t'.join([prog, version, self.progs[prog]]))
            if verbose:
                print(self.version_report[-1])


        if verbose:
            print()

        for line in warnings:
            print('WARNING:', line, file=sys.stderr)


        if len(errors):
            self.all_deps_ok = False

            for line in errors:
                print('ERROR:', line, file=sys.stderr)
            print('\nSomething wrong with at least one dependency. Please see the above error message(s)', file=sys.stderr)
            if fail_on_error:
                raise Error('Dependency error(s). Cannot continue')
        elif verbose:
            if len(warnings):
                print('\nWARNING: Required dependencies found, but at least one optional one was not. Please see previous warning(s) for more details.', file=sys.stderr)
            else:
                print('\nDependencies look OK')


    def exe(self, prog):
        return self.progs[prog]


    def version(self, prog):
        return self.versions[prog]


    @staticmethod
    def _get_exe(prog):
        '''Given a program name, return what we expect its exectuable to be called'''
        if prog in prog_to_env_var:
            env_var = prog_to_env_var[prog]
            if env_var in os.environ:
                return os.environ[env_var]

        return prog_to_default[prog]


    @staticmethod
    def _get_version(prog, path):
        '''Given a program name and expected path, tries to determine its version.
           Returns tuple (bool, version). First element True iff found version ok.
           Second element is version string (if found), otherwise an error message'''
        assert prog in prog_to_version_cmd
        args, regex = prog_to_version_cmd[prog]
        cmd = path + ' ' + args
        if prog == 'spades':
            cmd_output = subprocess.Popen(['python3', path, args], shell=False, stdout=subprocess.PIPE,
                                          stderr=subprocess.PIPE).communicate()
        else:
            cmd_output = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()

        cmd_output = common.decode(cmd_output[0]).split('\n')[:-1] + common.decode(cmd_output[1]).split('\n')[:-1]

        for line in cmd_output:
            hits = regex.search(line)
            if hits:
                return True, hits.group(1)

        return False, 'I tried to get the version of ' + prog + ' with: "' + cmd + '" and the output didn\'t match this regular expression: "' + regex.pattern + '"'

