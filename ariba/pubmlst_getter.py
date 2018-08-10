import tempfile
import re
import time
import os
import urllib.request
import xml.etree.ElementTree as ET
import pyfastaq
from ariba import common

class Error (Exception): pass


class PubmlstGetter:
    def __init__(self, debug=False, xml_file=None, verbose=False):
        self.debug = debug
        self.verbose = verbose

        if xml_file is None:
            self.xml_tree = self._get_xml_file_tree()
        else:
            self.xml_tree = ET.parse(xml_file)


    def _get_xml_file_tree(self):
        xml_url = 'http://pubmlst.org/data/dbases.xml'
        tmpdir = tempfile.mkdtemp(prefix='tmp.get_pubmlst_xml', dir=os.getcwd())
        xml_file = os.path.join(tmpdir, 'out.xml')
        self._download_file(xml_url, xml_file)
        xml_tree = ET.parse(xml_file)

        if not self.debug:
            common.rmtree(tmpdir)

        return xml_tree


    def _download_file(self, url, outfile):
        if self.verbose:
            print('Downloading "', url, '" and saving as "', outfile, '" ...', end='', sep='', flush=True)
        max_attempts = 3
        sleep_time = 3
        for i in range(max_attempts):
            time.sleep(sleep_time)
            try:
                urllib.request.urlretrieve(url, filename=outfile)
            except:
                continue
            break
        else:
            raise Error('Error downloading: ' + url)

        if self.verbose:
            print(' done', flush=True)


    def _get_species_list(self):
        return [x.text.rstrip() for x in self.xml_tree.getroot().findall('species')]


    def _get_profile_and_fasta_urls(self, species):
        species_dict = {x.text.rstrip(): x for x in self.xml_tree.getroot().findall('species')}
        if species not in species_dict:
            raise Error('Error! Species "' + species + '" not found. Cannot continue. Available species:\n' + '\n'.join(sorted(list(species_dict.keys()))))

        try:
            profile_url = species_dict[species].find('mlst').find('database').find('profiles').find('url').text
        except:
            raise Error('Error getting profile url for species ' + species + '. Cannot continue')

        locus_list = species_dict[species].find('mlst').find('database').find('loci').findall('locus')
        fasta_urls = [x.find('url').text for x in locus_list]

        if len(fasta_urls) == 0:
            raise Error('Error! No fasta files found for species ' + species + '. Cannot continue')

        return profile_url, fasta_urls


    @classmethod
    def _rename_seqs_in_fasta(cls, infile, outfile):
        f = pyfastaq.utils.open_file_write(outfile)
        file_reader = pyfastaq.sequences.file_reader(infile)
        nodot_regex = re.compile(r'^.*(?P<separator>[^.0-9])[0-9]+$')

        for seq in file_reader:
            if seq.id.startswith('Oxf.'):
                seq.id = 'Oxf_' + seq.id[4:]

            regex_match = nodot_regex.match(seq.id)
            if regex_match is not None:
                seq.id = '.'.join(seq.id.rsplit(regex_match.groupdict()['separator'], maxsplit=1))

            print(seq, file=f)

        pyfastaq.utils.close(f)


    def _download_profile_and_fastas(self, outdir, profile_url, fasta_urls):
        try:
            os.mkdir(outdir)
        except:
            raise Error('Error mkdir ' + outdir)

        profile_outfile = os.path.join(outdir, 'profile.txt')
        self._download_file(profile_url, profile_outfile)

        for fasta_url in fasta_urls:
            outfile = os.path.join(outdir, fasta_url.split('/')[-1])
            self._download_file(fasta_url, outfile + '.tmp')
            PubmlstGetter._rename_seqs_in_fasta(outfile + '.tmp', outfile)
            os.unlink(outfile + '.tmp')


    def print_available_species(self):
        species_list = self._get_species_list()
        print(*species_list, sep='\n')


    def get_species_files(self, species, outdir):
        profile_url, fasta_urls = self._get_profile_and_fasta_urls(species)
        self._download_profile_and_fastas(outdir, profile_url, fasta_urls)

