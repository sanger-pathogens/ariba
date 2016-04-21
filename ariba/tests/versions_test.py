import unittest
from ariba import versions

class TestVersions(unittest.TestCase):
    def test_get_all_versions(self):
        '''Test get_all_versions'''
        versions.get_all_versions(None)

