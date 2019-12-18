import unittest
import os
import filecmp
from ariba import pubmlst_getter

modules_dir = os.path.dirname(os.path.abspath(pubmlst_getter.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


class TestPubmlstGetter(unittest.TestCase):
    def setUp(self):
        xml_file = os.path.join(data_dir, 'pubmlst_getter.dbases.xml')
        self.pmlst_getter = pubmlst_getter.PubmlstGetter(xml_file=xml_file)


    def test_get_species_list(self):
        '''test _get_species_list'''
        expected = ['Achromobacter spp.', 'Acinetobacter baumannii#1', 'Acinetobacter baumannii#2', 'Aeromonas spp.']
        got = self.pmlst_getter._get_species_list()
        self.assertEqual(expected, got)


    def test_get_profile_and_fasta_urls(self):
        '''test _get_profile_and_fasta_urls'''
        expected_profile = 'http://pubmlst.org/data/profiles/abaumannii_2.txt'
        expected_fasta = [
            'http://pubmlst.org/data/alleles/abaumannii_2/Pas_cpn60.tfa',
            'http://pubmlst.org/data/alleles/abaumannii_2/Pas_fusA.tfa',
            'http://pubmlst.org/data/alleles/abaumannii_2/Pas_gltA.tfa',
            'http://pubmlst.org/data/alleles/abaumannii_2/Pas_pyrG.tfa',
            'http://pubmlst.org/data/alleles/abaumannii_2/Pas_recA.tfa',
            'http://pubmlst.org/data/alleles/abaumannii_2/Pas_rplB.tfa',
            'http://pubmlst.org/data/alleles/abaumannii_2/Pas_rpoB.tfa',
        ]

        got_profile, got_fasta = self.pmlst_getter._get_profile_and_fasta_urls('Acinetobacter baumannii#2')
        self.assertEqual(expected_profile, got_profile)
        self.assertEqual(expected_fasta, got_fasta)

    def test_rename_seqs_in_fasta(self):
        '''test _rename_seqs_in_fasta'''
        infile = os.path.join(data_dir, 'pubmlst_rename_seqs.in.fa')
        expected_file = os.path.join(data_dir, 'pubmlst_rename_seqs.expected.fa')
        tmp_out = 'tmp.test.pubmlst_rename_seqs.out.fa'
        pubmlst_getter.PubmlstGetter._rename_seqs_in_fasta(infile, tmp_out)
        self.assertTrue(filecmp.cmp(expected_file, tmp_out, shallow=False))
        os.unlink(tmp_out)

