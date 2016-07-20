import sys
from ariba import versions

def run(options):
    extern_progs, report_lines = versions.get_all_versions(raise_error=False)
    print(*report_lines, sep='\n')
