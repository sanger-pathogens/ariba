import unittest
import sys
import os
import shutil
import filecmp
import pyfastaq
from ariba import ref_preparer

modules_dir = os.path.dirname(os.path.abspath(ref_preparer.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')


class TestRefPreparer(unittest.TestCase):
    def test_get_ref_files_ref_prefix(self):
        '''test _get_ref_files using ref_prefix'''
        ref_prefix = os.path.abspath('tmp.test.ref_preparer')
        presabs = ref_prefix + '.presence_absence.fa'
        varonly = ref_prefix + '.variants_only.fa'
        noncoding = ref_prefix + '.noncoding.fa'
        metadata = ref_prefix + '.metadata.tsv'

        for filename in (presabs, varonly, noncoding, metadata):
            with open(filename, 'w') as f:
                pass

        got = ref_preparer.RefPreparer._get_ref_files(ref_prefix, None, None, None, None)
        expected = {
            'presabs': presabs,
            'varonly': varonly,
            'noncoding': noncoding,
            'metadata': metadata,
        }

        self.assertEqual(expected, got)

        os.unlink(metadata)
        got = ref_preparer.RefPreparer._get_ref_files(ref_prefix, None, None, None, None)
        expected['metadata'] = None
        self.assertEqual(expected, got)

        os.unlink(presabs)
        got = ref_preparer.RefPreparer._get_ref_files(ref_prefix, None, None, None, None)
        expected['presabs'] = None
        self.assertEqual(expected, got)

        os.unlink(varonly)
        got = ref_preparer.RefPreparer._get_ref_files(ref_prefix, None, None, None, None)
        expected['varonly'] = None
        self.assertEqual(expected, got)

        os.unlink(noncoding)
        with self.assertRaises(ref_preparer.Error):
            ref_preparer.RefPreparer._get_ref_files(ref_prefix, None, None, None, None)


    def test_get_ref_files_naming_each_file(self):
        ref_prefix = os.path.abspath('tmp.test.ref_preparer')
        presabs = ref_prefix + '.presence_absence.fa'
        varonly = ref_prefix + '.variants_only.fa'
        noncoding = ref_prefix + '.noncoding.fa'
        metadata = ref_prefix + '.metadata.tsv'
        not_a_file = 'notafile'
        self.assertFalse(os.path.exists(not_a_file))

        for filename in (presabs, varonly, noncoding, metadata):
            with open(filename, 'w') as f:
                pass

        got = ref_preparer.RefPreparer._get_ref_files(None, presabs, varonly, noncoding, metadata)
        expected = {
            'presabs': presabs,
            'varonly': varonly,
            'noncoding': noncoding,
            'metadata': metadata,
        }

        self.assertEqual(expected, got)

        got = ref_preparer.RefPreparer._get_ref_files(None, presabs, varonly, noncoding, None)
        expected['metadata'] = None
        self.assertEqual(expected, got)

        got = ref_preparer.RefPreparer._get_ref_files(None, None, varonly, noncoding, None)
        expected['presabs'] = None
        self.assertEqual(expected, got)

        got = ref_preparer.RefPreparer._get_ref_files(None, None, None, noncoding, None)
        expected['varonly'] = None
        self.assertEqual(expected, got)

        with self.assertRaises(ref_preparer.Error):
            got = ref_preparer.RefPreparer._get_ref_files(None, None, None, None, None)

        with self.assertRaises(ref_preparer.Error):
            got = ref_preparer.RefPreparer._get_ref_files(None, not_a_file, None, noncoding, None)

        with self.assertRaises(ref_preparer.Error):
            got = ref_preparer.RefPreparer._get_ref_files(None, presabs, not_a_file, None, None)

        with self.assertRaises(ref_preparer.Error):
            got = ref_preparer.RefPreparer._get_ref_files(None, presabs, None, not_a_file, None)

        os.unlink(presabs)
        os.unlink(varonly)
        os.unlink(noncoding)
        os.unlink(metadata)
