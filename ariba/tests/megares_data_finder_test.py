import unittest
import os
from ariba import megares_data_finder

modules_dir = os.path.dirname(os.path.abspath(megares_data_finder.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


class TestMegaresDataFinder(unittest.TestCase):
    def test_zips_from_index_page_string(self):
        '''test _zips_from_index_page_string'''
        html_string = r''''<!doctype html>
<html>
  <head>

  </head>

<ul>
  <li><a href="/downloads/megares_v3.00.zip">All Files</a> (.zip)</li>
  <li><a href="/downloads/megares_v2.00.zip">All Files</a> (.zip)</li>
  <li><a href="/downloads/megares_v1.01.zip">All Files</a> (.zip)</li>
  <li><a href="/downloads/megares_v1.00.zip">All Files</a> (.zip)</li>

</html>'''

        expected = {'3.00': '/downloads/megares_v3.00.zip', '2.00': '/downloads/megares_v2.00.zip', '1.01': '/downloads/megares_v1.01.zip', '1.00': '/downloads/megares_v1.00.zip'}
        got = megares_data_finder.MegaresDataFinder._zips_from_index_page_string(html_string)
        self.assertEqual(expected, got)


    def test_get_url_for_version(self):
        '''test _get_url_for_version'''
        zips = {'3.00': '/downloads/megares_v3.00.zip', '2.00': '/downloads/megares_v2.00.zip', '1.01': '/downloads/megares_v1.01.zip', '1.00': '/downloads/megares_v1.00.zip'}
        self.assertEqual('/downloads/megares_v3.00.zip', megares_data_finder.MegaresDataFinder._get_url_for_version(zips))
        self.assertEqual('/downloads/megares_v2.00.zip', megares_data_finder.MegaresDataFinder._get_url_for_version(zips, version='2.00'))
        with self.assertRaises(megares_data_finder.Error):
            self.assertEqual('megares_v1.00.zip', megares_data_finder.MegaresDataFinder._get_url_for_version(zips, version='0.42'))

