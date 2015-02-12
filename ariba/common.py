import sys
import subprocess

version = '0.0.1'

def syscall(cmd, allow_fail=False, verbose=False):
    if verbose:
        print('syscall:', cmd, flush=True)
    try:
        subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as error:
        errors = error.output.decode()
        print('The following command failed with exit code', error.returncode, file=sys.stderr)
        print(cmd, file=sys.stderr)
        print('\nThe output was:\n', file=sys.stderr)
        print(errors, file=sys.stderr)

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
