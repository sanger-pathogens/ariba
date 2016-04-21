import sys
from ariba import versions

def run():
    versions.get_all_versions(sys.stdout, raise_error=False)
