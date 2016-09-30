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
    #'cdhit2d': 'cd-hit-est-2d',
    #'gapfiller': 'GapFiller.pl',
    'mash': 'mash',
    'nucmer' : 'nucmer',
    #'spades': 'spades.py',
    #'sspace': 'SSPACE_Basic_v2.0.pl',
}


prog_to_env_var = {x: 'ARIBA_' + x.upper() for x in prog_to_default if x not in {'nucmer'}}


prog_to_version_cmd = {
    'bowtie2': ('--version', re.compile('.*bowtie2.*version (.*)$')),
    'cdhit': ('', re.compile('CD-HIT version ([0-9\.]+) \(')),
    #'cdhit2d': ('', re.compile('CD-HIT version ([0-9\.]+) \(')),
    #'gapfiller': ('', re.compile('^Usage: .*pl \[GapFiller_(.*)\]')),
    'mash': ('', re.compile('^Mash version (.*)$')),
    'nucmer': ('--version', re.compile('^NUCmer \(NUCleotide MUMmer\) version ([0-9\.]+)')),
    #'spades': ('', re.compile('^SPAdes genome assembler v\.?([0-9\.]+)')),
    #'sspace': ('', re.compile('^Usage: .*pl \[SSPACE_(.*)\]')),
}


min_versions = {
    'bowtie2': '2.1.0',
    'cdhit': '4.6',
    #'cdhit2d': '4.6',
    'mash': '1.0.2',
    'nucmer': '3.1',
    #'spades': '3.5.0',
}


class ExternalProgs:
    def __init__(self, verbose=False, fail_on_error=True):
        optional_progs = {'sspace', 'gapfiller'}
        self.progs = {}
        self.version_report = []
        self.all_deps_ok = True

        if verbose:
            print('{:_^79}'.format(' Checking dependencies and their versions '))

        errors = []
        warnings = []

        for prog in sorted(prog_to_default):
            prog_exe = self._get_exe(prog)
            self.progs[prog] = shutil.which(prog_exe)
            # Travis is using python3.4, and actually "python" in travis means
            # python3.4, not python2. SPAdes throws an error about not being
            # compatible with python3.4.
            # This means we need to explicitly run SPAdes with python2.
            if prog == 'spades' and self.progs[prog] is not None:
                self.progs[prog] = 'python2 ' + self.progs[prog]
            if self.progs[prog] is None:
                if prog in optional_progs:
                    warnings.append(prog + ' not found in path. Looked for ' + prog_exe + '. But it is optional so will be skipped during assembly')
                else:
                    errors.append(prog + ' not found in path. Looked for ' + prog_exe)

                self.version_report.append('\t'.join([prog, 'NA', 'NOT_FOUND']))
                if verbose:
                    print(self.version_report[-1])
                continue
            elif prog in {'sspace', 'gapfiller'}:
                self.progs[prog] = os.path.realpath(self.progs[prog])

            got_version, version = self._get_version(prog, self.progs[prog])

            if got_version:
                if prog in min_versions and LooseVersion(version) < LooseVersion(min_versions[prog]):
                    errors.append(' '.join(['Found version', version, 'of', prog, 'which is too low! Please update to at least', min_versions[prog] + '. Found it here:', prog_exe]))
            else:
                errors.append(version)
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
                return True, hits.group(1)

        return False, 'I tried to get the version of ' + prog + ' with: "' + cmd + '" and the output didn\'t match this regular expression: "' + regex.pattern + '"'

