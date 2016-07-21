import sys
from distutils.version import LooseVersion
from ariba import external_progs
from ariba import __version__ as ariba_version


package_min_versions = {
    'dendropy': '4.1.0',
    'pyfastaq': '3.12.0',
    'pysam': '0.8.1',
    'pymummer' : '0.7.1',
}

package_max_versions = {
    #'pysam': '0.8.3',
}


def get_all_versions(raise_error=True):
    extern_progs = external_progs.ExternalProgs(fail_on_error=False)

    report_lines = [
        'ARIBA version: ' + ariba_version,
        '\nExternal dependencies:',
        '\n'.join(extern_progs.version_report),
        '\nExternal dependencies OK: ' + str(extern_progs.all_deps_ok),
        '\nPython version:',
        str(sys.version),
        '\nPython packages:',
    ]

    python_packages_ok = True

    for package in ['ariba', 'dendropy', 'pyfastaq', 'pymummer', 'pysam']:
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

        report_lines.append(package + '\t' + version + '\t' + path)

    all_ok = extern_progs.all_deps_ok and python_packages_ok

    report_lines.extend([
        '\nPython packages OK: ' + str(python_packages_ok),
        '\nEverything looks OK: ' + str(all_ok),
    ])

    if raise_error and not all_ok:
        print(*report_lines, sep='\n', file=sys.stderr)
        print('Some dependencies not satisfied. Cannot continue.', file=sys.stderr)
        sys.exit(1)

    return extern_progs, report_lines
