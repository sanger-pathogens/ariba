import unittest
import os
from ariba import mlst_profile

modules_dir = os.path.dirname(os.path.abspath(mlst_profile.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


class TestMlstProfile(unittest.TestCase):
    def test_init(self):
        '''test init'''
        infile = os.path.join(data_dir, 'mlst_profile_test.init.profile.tsv')
        profile = mlst_profile.MlstProfile(infile)
        expected_genes = ['nusA', 'rpoB', 'eno', 'gltB', 'lepA', 'nuoL', 'nrdA']
        self.assertEqual(expected_genes, profile.genes_list)
        self.assertEqual(set(expected_genes), profile.genes_set)
        expected_dict = {
            (1, 26, 2, 2, 59, 8, 1): 1,
            (1, 26, 2, 4, 59, 2, 5): 2
        }
        self.assertEqual(expected_dict, profile.profile_to_type)
		
    def test_init_multiple_extra_columns(self):
        '''test init'''
        infile = os.path.join(data_dir, 'mlst_profile_test.init_multiple_extra_columns.profile.tsv')
        profile = mlst_profile.MlstProfile(infile)
        expected_genes = ['nusA', 'rpoB', 'eno', 'gltB', 'lepA', 'nuoL', 'nrdA']
        self.assertEqual(expected_genes, profile.genes_list)
        self.assertEqual(set(expected_genes), profile.genes_set)
        expected_dict = {
            (1, 26, 2, 2, 59, 8, 1): 1,
            (1, 26, 2, 4, 59, 2, 5): 2
        }
        self.assertEqual(expected_dict, profile.profile_to_type)


    def test_has_gene(self):
        '''Test has_gene'''
        infile = os.path.join(data_dir, 'mlst_profile_test.profile.tsv')
        profile = mlst_profile.MlstProfile(infile)
        for gene in ['nusA', 'rpoB', 'eno', 'gltB', 'lepA', 'nuoL', 'nrdA']:
            self.assertTrue(profile.has_gene(gene))
        self.assertFalse(profile.has_gene('not there'))


    def test_get_sequence_type(self):
        '''Test get_sequence_type'''
        infile = os.path.join(data_dir, 'mlst_profile_test.profile.tsv')
        profile = mlst_profile.MlstProfile(infile)
        tests = [
            ('ND', {}),
            ('ND', {'foo': 1}),
            (1, {'nusA': 1, 'rpoB': 26, 'eno': 2, 'gltB': 2, 'lepA': 59, 'nuoL': 8, 'nrdA': 1}),
            (2, {'nusA': 1, 'rpoB': 26, 'eno': 2, 'gltB': 4, 'lepA': 59, 'nuoL': 2, 'nrdA': 5}),
            ('ND', {'nusA': 'ND', 'rpoB': 26, 'eno': 2, 'gltB': 4, 'lepA': 59, 'nuoL': 2, 'nrdA': 5}),
            ('Novel', {'nusA': 1000, 'rpoB': 26, 'eno': 2, 'gltB': 4, 'lepA': 59, 'nuoL': 2, 'nrdA': 5}),
        ]
        for expected, profile_dict in tests:
            self.assertEqual(expected, profile.get_sequence_type(profile_dict))

