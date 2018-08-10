import os
import time
import sys
import subprocess
import urllib.request
import pyfastaq


class Error (Exception): pass


def syscall(cmd, allow_fail=False, verbose=False, verbose_filehandle=sys.stdout, print_errors=True, shell=True):
    if verbose:
        print('syscall:', cmd, flush=True, file=verbose_filehandle)
        if not shell:
            print('syscall string:', " ".join('"{}"'.format(_) for _ in cmd), flush=True, file=verbose_filehandle)
    try:
        subprocess.check_output(cmd, shell=shell, stderr=subprocess.STDOUT)
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
    except Exception as msg:
        print("Unexpected exception: ", msg, file=sys.stderr)
        raise
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


def download_file(url, outfile, max_attempts=3, sleep_time=2, verbose=False):
    if verbose:
        print('Downloading "', url, '" and saving as "', outfile, '" ...', end='', sep='', flush=True)

    for i in range(max_attempts):
        time.sleep(sleep_time)
        try:
            urllib.request.urlretrieve(url, filename=outfile)
        except:
            continue
        break
    else:
        raise Error('Error downloading: ' + url)

    if verbose:
        print(' done', flush=True)


def rmtree(input_dir):
    '''Does rm -r on input_dir. Meant to replace shutil.rmtree,
    which seems to be causing issues with files not getting deleted
    and the directory non-empty afterwards'''
    syscall('rm -rf ' + input_dir)
