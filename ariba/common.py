import os
import sys
import subprocess
import pyfastaq

def syscall(cmd, allow_fail=False, verbose=False, verbose_filehandle=sys.stdout, print_errors=True):
    if verbose:
        print('syscall:', cmd, flush=True, file=verbose_filehandle)
    try:
        subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as error:
        errors = error.output.decode()
        if print_errors:
            print('The following command failed with exit code', error.returncode, file=sys.stderr)
            print(cmd, file=sys.stderr)
            print('\nThe output was:\n', file=sys.stderr)
            print(errors, file=sys.stderr, flush=True)

        if allow_fail:
            return False, errors
        else:
            sys.exit(1)

    return True, None


def decode(x):
    try:
        s = x.decode()
    except:
        return x
    return s


def cat_files(infiles, outfile):
    '''Cats all files in list infiles into outfile'''
    f_out = pyfastaq.utils.open_file_write(outfile)

    for filename in infiles:
        if os.path.exists(filename):
            f_in = pyfastaq.utils.open_file_read(filename)
            for line in f_in:
                print(line, end='', file=f_out)
            pyfastaq.utils.close(f_in)

    pyfastaq.utils.close(f_out)
