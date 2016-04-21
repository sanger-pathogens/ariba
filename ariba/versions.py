import sys
from distutils.version import LooseVersion
from ariba import external_progs
from ariba import __version__ as ariba_version


package_min_versions = {
    'openpyxl': '1.6.2',
    'pyfastaq': '3.12.0',
    'pysam': '0.8.1',
    'pymummer' : '0.7.1',
}

package_max_versions = {
    'pysam': '0.8.3',
}


def get_all_versions(filehandle, raise_error=True):
    if filehandle is not None:
        print('ARIBA version:', ariba_version, file=filehandle)
        print('\n\nExternal dependencies:', file=filehandle)

    extern_progs = external_progs.ExternalProgs(fail_on_error=False)

    if filehandle is not None:
        print(*extern_progs.version_report, sep='\n', file=filehandle)
        print('\nExternal dependencies OK:', extern_progs.all_deps_ok, file=filehandle)
        print('\n\nPython version:', file=filehandle)
        print(sys.version, file=filehandle)
        print('\nPython packages:', file=filehandle)

    python_packages_ok = True

    for package in ['openpyxl', 'pyfastaq', 'pymummer', 'pysam']:
        try:
            exec('import ' + package)
            version = eval(package + '.__version__')
            path = eval(package + '.__file__')
        except:
            version = 'NOT_FOUND'
            path = 'NOT_FOUND'
            python_packages_ok = False

        if version != 'NOT_FOUND':
            if package in package_min_versions and LooseVersion(version) < package_min_versions[package]:
                version += '... THIS IS TOO LOW. Needs>=' + package_min_versions[package]
                python_packages_ok = False
            elif package in package_max_versions and LooseVersion(version) > package_max_versions[package]:
                version += '...THIS IS TOO HIGH. Needs <=' + package_max_versions[package]
                python_packages_ok = False

        if filehandle is not None:
            print(package + '\t' + version + '\t' + path, file=filehandle)

    all_ok = extern_progs.all_deps_ok and python_packages_ok

    if filehandle is not None:
        print('\nPython packages OK:', python_packages_ok, file=filehandle)
        print('\nEverything looks OK:', all_ok, file=filehandle)

    if raise_error and not all_ok:
        print('Some dependencies not satisfied. Cannot continue.', file=sys.stderr)
        sys.exit(1)

