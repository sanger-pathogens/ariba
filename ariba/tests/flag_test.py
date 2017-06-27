import unittest
import os
from ariba import flag

modules_dir = os.path.dirname(os.path.abspath(flag.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')

class TestFlag(unittest.TestCase):
    def test_init_and_to_number(self):
        '''Test __init__ and to_number'''
        for i in range(512):
            f = flag.Flag(i)
            self.assertEqual(f.to_number(), i)


    def test_set_flag(self):
        '''Test set_flag'''
        for i in range(512):
            f = flag.Flag()
            f.set_flag(i)
            self.assertEqual(f.to_number(), i)


    def test_add(self):
        '''Test add'''
        f = flag.Flag()
        expected = [1, 3, 7, 15, 31, 63, 127, 255, 511, 1023, 2047]
        for i in range(len(flag.flags_in_order)):
            f.add(flag.flags_in_order[i])
            self.assertEqual(f.to_number(), expected[i])


    def test_str(self):
        '''Test __str__'''
        for i in range(512):
            f = flag.Flag(i)
            self.assertEqual(str(f), str(i))


    def test_to_long_str(self):
        '''Test to_long_str'''
        f = flag.Flag(13)
        expected = '\n'.join([
            '[X] assembled',
            '[ ] assembled_into_one_contig',
            '[X] region_assembled_twice',
            '[X] complete_gene',
            '[ ] unique_contig',
            '[ ] scaffold_graph_bad',
            '[ ] assembly_fail',
            '[ ] variants_suggest_collapsed_repeat',
            '[ ] hit_both_strands',
            '[ ] has_variant',
            '[ ] ref_seq_choose_fail',
        ])

        self.assertEqual(expected, f.to_long_string())


    def test_has(self):
        '''Test has'''
        for x in flag.flags_in_order:
            f = flag.Flag(0)
            self.assertFalse(f.has(x))
            f.add(x)
            self.assertTrue(f.has(x))


    def test_to_comma_separated_string(self):
        '''Test to_comma_separated_string'''
        f = flag.Flag(27)
        expected = 'assembled,assembled_into_one_contig,complete_gene,unique_contig'
        self.assertEqual(expected, f.to_comma_separated_string())

