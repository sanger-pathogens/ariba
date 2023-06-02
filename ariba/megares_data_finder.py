import urllib.request
from bs4 import BeautifulSoup
from distutils.version import LooseVersion


class Error (Exception): pass


class MegaresDataFinder:
    def __init__(self, version=None):
        self.url_root = 'https://www.meglab.org/megares/download'
        self.url_root_downloads='https://www.meglab.org'
        self.index_url = self.url_root 
        self.version = version


    def _get_available_zips(self):
        try:
            response = urllib.request.urlopen(self.index_url)
            html_text = response.read()
        except:
            raise Error('Error getting megares download page ' + self.index_url)

        return html_text


    @classmethod
    def _zips_from_index_page_string(cls, html_text):
        try:
            soup = BeautifulSoup(html_text, 'html.parser')
        except:
            raise Error('Error parsing contents of megares download page. Cannot continue')

        prefix = '/downloads/megares_v'
        suffix = '.zip'
        zips = {}

        for link in soup.find_all('a'):
            href = link.get('href')
            if href.startswith(prefix) and href.endswith(suffix):
                version = href[len(prefix):-len(suffix)]
                zips[version] = href
                

        #print(zips)
        return zips


    @classmethod
    def _get_url_for_version(cls, zips, version=None):
        if version is None:
            versions = list(zips.keys())
            versions.sort(key=LooseVersion)
            #print('This is v',zips[versions[-1]])
            return zips[versions[-1]]
        else:
            try:
                return zips[version]
            except:
                versions = ', '.join(list(zips.keys()))
                raise Error('Error! version ' + version + ' of megares not found. Available versions: ' + versions)


    def run(self):
        print('Finding available megares versions from', self.index_url)
        html_text = self._get_available_zips()
        zips = MegaresDataFinder._zips_from_index_page_string(html_text)
        print('Found versions: ', ', '.join(list(zips.keys())))
        url = MegaresDataFinder._get_url_for_version(zips, version=self.version)
        return self.url_root_downloads + url


