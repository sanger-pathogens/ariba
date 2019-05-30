import unittest
import json
import shutil
import os
import pickle
import pyfastaq
import filecmp
from ariba import clusters, common, external_progs, histogram, sequence_metadata

modules_dir = os.path.dirname(os.path.abspath(clusters.__file__))
data_dir = os.path.join(modules_dir, 'tests', 'data')
extern_progs = external_progs.ExternalProgs()


def file_to_list(infile):
    f = pyfastaq.utils.open_file_read(infile)
    lines = [x for x in f.readlines()]
    pyfastaq.utils.close(f)
    return lines


class TestClusters(unittest.TestCase):
    def setUp(self):
        self.cluster_dir = 'tmp.Cluster'
        self.refdata_dir = 'tmp.RefData'
        os.mkdir(self.refdata_dir)
        shutil.copyfile(os.path.join(data_dir, 'clusters_test_dummy_db.fa'), os.path.join(self.refdata_dir, '02.cdhit.all.fa'))
        shutil.copyfile(os.path.join(data_dir, 'clusters_test_dummy_db.tsv'), os.path.join(self.refdata_dir, '01.filter.check_metadata.tsv'))
        with open(os.path.join(self.refdata_dir, '00.info.txt'), 'w') as f:
            print('genetic_code\t11', file=f)

        with open(os.path.join(self.refdata_dir, '02.cdhit.clusters.pickle'), 'wb') as f:
            pickle.dump({'0': {'x'}}, f)

        reads1 = os.path.join(data_dir, 'clusters_test_dummy_reads_1.fq')
        reads2 = os.path.join(data_dir, 'clusters_test_dummy_reads_2.fq')
        self.clusters = clusters.Clusters(self.refdata_dir, reads1, reads2, self.cluster_dir, extern_progs, clean=False)


    def tearDown(self):
        common.rmtree(self.cluster_dir)
        common.rmtree(self.refdata_dir)


    def test_load_reference_data_info_file(self):
        '''test _load_reference_data_info_file'''
        infile = os.path.join(data_dir, 'clusters_test_load_data_info_file')
        expected = {'genetic_code': 11}
        got = clusters.Clusters._load_reference_data_info_file(infile)
        self.assertEqual(expected, got)


    def test_load_ref_data_from_dir(self):
        '''test _load_reference_data_from_dir'''
        indir = os.path.join(data_dir, 'clusters_load_ref_data_from_dir')
        got_refdata, got_clusters = clusters.Clusters._load_reference_data_from_dir(indir)
        expected_seq_dict = {
            'variants_only1': pyfastaq.sequences.Fasta('variants_only1', 'atggcgtgcgatgaataa'),
            'presabs1': pyfastaq.sequences.Fasta('presabs1', 'atgatgatgagcccggcgatggaaggcggctag'),
            'noncoding1': pyfastaq.sequences.Fasta('noncoding1', 'ACGTA'),
        }
        self.assertEqual(expected_seq_dict, got_refdata.sequences)
        self.assertEqual(11, got_refdata.genetic_code)

        expected_metadata = {
            'presabs1': {
                'seq_type': 'p',
                'variant_only': False,
                '.': {sequence_metadata.SequenceMetadata('presabs1\t1\t0\t.\t.\tpresabs1 description')},
                'n': {},
                'p': {}
            },
            'variants_only1': {
                'seq_type': 'p',
                'variant_only': True,
                '.': set(),
                'n': {},
                'p': {1: {sequence_metadata.SequenceMetadata('variants_only1\t1\t1\tC2I\t.\tdescription of variants_only1 C2I')}}
            },
            'noncoding1': {
                'seq_type': 'n',
                'variant_only': False,
                '.': {sequence_metadata.SequenceMetadata('noncoding1\t0\t0\t.\t.\t.')},
                'n': {},
                'p': {},
            }
        }
        self.assertEqual(expected_metadata, got_refdata.metadata)

        expected_clusters = {'0': {'presabs1'}, '1': {'variants_only1'}, '2': {'noncoding1'}}
        self.assertEqual(expected_clusters, got_clusters)

        self.assertEqual({'a': 'b'}, got_refdata.extra_parameters)

    def test_minimap_reads_to_all_ref_seqs(self):
        '''test test_minimap_reads_to_all_ref_seqs'''
        clusters_tsv = os.path.join(data_dir, 'clusters_minimap_reads_to_all_refs.clstrs.tsv')
        ref_fasta = os.path.join(data_dir, 'clusters_minimap_reads_to_all_refs.ref.fa')
        reads_1 = os.path.join(data_dir, 'clusters_minimap_reads_to_all_refs.reads_1.fq')
        reads_2 = os.path.join(data_dir, 'clusters_minimap_reads_to_all_refs.reads_2.fq')
        tmp_outprefix = 'tmp.clusters_test_minimap_reads_to_all_ref_seqs'
        clusters.Clusters._minimap_reads_to_all_ref_seqs(clusters_tsv, ref_fasta, reads_1, reads_2, tmp_outprefix)
        expected_cluster2rep = os.path.join(data_dir, 'clusters_minimap_reads_to_all_refs.out.clstr2rep')
        expected_cluster_counts = os.path.join(data_dir, 'clusters_minimap_reads_to_all_refs.out.clstr_count')
        expected_proper_pairs = os.path.join(data_dir, 'clusters_minimap_reads_to_all_refs.out.pairs')
        expected_insert_hist = os.path.join(data_dir, 'clusters_minimap_reads_to_all_refs.out.hist')

        # not sure that the reads order is preserved, so just check read store file exists
        self.assertTrue(os.path.exists(os.path.join(tmp_outprefix + '.reads')))

        self.assertTrue(filecmp.cmp(expected_cluster2rep, tmp_outprefix + '.cluster2representative', shallow=False))
        self.assertTrue(filecmp.cmp(expected_cluster_counts, tmp_outprefix + '.clusterCounts', shallow=False))
        self.assertTrue(filecmp.cmp(expected_proper_pairs, tmp_outprefix + '.properPairs', shallow=False))
        self.assertTrue(filecmp.cmp(expected_insert_hist, tmp_outprefix + '.insertHistogram', shallow=False))
        os.unlink(tmp_outprefix + '.cluster2representative')
        os.unlink(tmp_outprefix + '.clusterCounts')
        os.unlink(tmp_outprefix + '.properPairs')
        os.unlink(tmp_outprefix + '.insertHistogram')
        os.unlink(tmp_outprefix + '.reads')


    def test_load_minimap_out_cluster2representative(self):
        '''test _load_minimap_out_cluster2representative'''
        infile = os.path.join(data_dir, 'clusters_test_load_minimap_out_cluster2representative.in')
        got = clusters.Clusters._load_minimap_out_cluster2representative(infile)
        expected = {'1': 'ref2', '2': 'ref42'}
        self.assertEqual(expected, got)


    def test_load_minimap_out_cluster_counts(self):
        '''test _load_minimap_out_cluster_counts'''
        infile = os.path.join(data_dir, 'clusters_test_load_minimap_out_cluster_counts.in')
        got_read_count, got_base_count = clusters.Clusters._load_minimap_out_cluster_counts(infile)
        expected_read_count = {'1': 42, '2': 43}
        expected_base_count = {'1': 4242, '2': 4343}
        self.assertEqual(got_read_count, expected_read_count)
        self.assertEqual(got_base_count, expected_base_count)


    def test_load_minimap_insert_histogram(self):
        '''test _load_minimap_insert_histogram'''
        infile = os.path.join(data_dir, 'clusters_test_load_minimap_insert_histogram.in')
        bin_size = 10
        got = clusters.Clusters._load_minimap_insert_histogram(infile, bin_size)
        expected = histogram.Histogram(bin_size)
        expected.add(85, count=1)
        expected.add(86, count=2)
        expected.add(90, count=4)
        expected.add(91, count=6)
        expected.add(97, count=10)
        expected.add(100, count=7)
        expected.add(111, count=3)
        self.assertEqual(expected, got)


    def test_load_minimap_proper_pairs(self):
        '''test _load_minimap_proper_pairs'''
        infile = os.path.join(data_dir, 'clusters_test_load_minimap_proper_pairs.in')
        got = clusters.Clusters._load_minimap_proper_pairs(infile)
        self.assertEqual(42424242, got)


    def test_load_minimap_files(self):
        '''test _load_minimap_files'''
        bin_size = 10
        inprefix = os.path.join(data_dir, 'clusters_test_load_minimap_files')
        got_clster2rep, got_cluster_read_count, got_cluster_base_count, got_insert_hist, got_proper_pairs = clusters.Clusters._load_minimap_files(inprefix, bin_size)
        expected_clster2rep = {'1': 'ref2', '2': 'ref42'}
        expected_cluster_read_count = {'1': 42, '2': 43}
        expected_cluster_base_count = {'1': 4242, '2': 4343}
        expected_insert_hist_bins = {80: 3, 90: 20, 100: 7, 110: 3}
        expected_proper_pairs = 42424242
        self.assertEqual(expected_clster2rep, got_clster2rep)
        self.assertEqual(expected_cluster_read_count, got_cluster_read_count)
        self.assertEqual(expected_cluster_base_count, got_cluster_base_count)
        self.assertEqual(expected_insert_hist_bins, got_insert_hist.bins)
        self.assertEqual(expected_proper_pairs, got_proper_pairs)


    def test_set_insert_size_data(self):
        '''test _set_insert_size_data'''
        self.clusters.insert_hist.bins = {
            1: 1,
            2: 1,
            3: 3,
            4: 3,
            5: 5,
            6: 3,
            7: 2,
            8: 2,
            9: 1,
            10: 1,
        }
        self.clusters.insert_hist.bin_width=1

        self.clusters._set_insert_size_data()
        self.assertEqual(self.clusters.insert_size, 5.5)
        self.assertEqual(self.clusters.insert_sspace_sd, 0.91)


    def test_write_report(self):
        class FakeCluster:
            def __init__(self, lines):
                self.report_lines = lines

        clusters_dict = {
            'gene1': FakeCluster(['gene1\tline1']),
            'gene2': FakeCluster(['gene2\tline2'])
        }

        tmp_tsv = 'tmp.test_write_report.tsv'
        clusters.Clusters._write_report(clusters_dict, tmp_tsv)

        expected = os.path.join(data_dir, 'clusters_test_write_report.tsv')
        self.assertTrue(filecmp.cmp(expected, tmp_tsv, shallow=False))
        os.unlink(tmp_tsv)


    def test_write_catted_assemblies_fasta(self):
        '''test _write_catted_assemblies_fasta'''
        seq1 = pyfastaq.sequences.Fasta('seq1', 'ACGT')
        seq2 = pyfastaq.sequences.Fasta('seq2', 'TTTT')
        seq3 = pyfastaq.sequences.Fasta('seq3', 'AAAA')
        class FakeAssembly:
            def __init__(self, seqs):
                if seqs is not None:
                    self.sequences = {x.id: x for x in seqs}

        class FakeCluster:
            def __init__(self, seqs):
                self.assembly = FakeAssembly(seqs)

        self.clusters.clusters = {
            'cluster1': FakeCluster([seq1, seq2]),
            'cluster2': FakeCluster([seq3]),
            'cluster3': FakeCluster(None),
        }

        tmp_file = 'tmp.test_write_catted_assemblies_fasta.fa'
        self.clusters._write_catted_assemblies_fasta(tmp_file)
        expected = os.path.join(data_dir, 'clusters_test_write_catted_assemblies_fasta.expected.out.fa')
        self.assertTrue(filecmp.cmp(expected, tmp_file, shallow=False))
        os.unlink(tmp_file)


    def test_write_catted_assembled_seqs_fasta(self):
        '''test _write_catted_assembled_seqs_fasta'''
        seq1 = pyfastaq.sequences.Fasta('seq1', 'ACGT')
        seq2 = pyfastaq.sequences.Fasta('seq2', 'TTTT')
        seq3 = pyfastaq.sequences.Fasta('seq3', 'AAAA')
        class FakeAssemblyCompare:
            def __init__(self, assembled_seqs):
                if assembled_seqs is not None:
                    self.assembled_reference_sequences = {x.id: x for x in assembled_seqs}

        class FakeCluster:
            def __init__(self, assembled_seqs):
                self.assembly_compare = FakeAssemblyCompare(assembled_seqs)

        self.clusters.clusters = {
            'gene1': FakeCluster([seq1, seq2]),
            'gene2': FakeCluster([seq3]),
            'gene3': FakeCluster(None),
        }

        tmp_file = 'tmp.test_write_catted_assembled_seqs_fasta.fa'
        self.clusters._write_catted_assembled_seqs_fasta(tmp_file)
        expected = os.path.join(data_dir, 'clusters_test_write_catted_assembled_genes_fasta.expected.out.fa')
        self.assertTrue(filecmp.cmp(expected, tmp_file, shallow=False))
        os.unlink(tmp_file)


    def test_cat_genes_match_ref(self):
        '''test _write_catted_genes_matching_refs_fasta'''
        seq1 = pyfastaq.sequences.Fasta('seq1', 'ACGT')
        seq3 = pyfastaq.sequences.Fasta('seq3', 'AAAA')
        class FakeAssemblyCompare:
            def __init__(self, seq, seq_type, start, end):
                self.gene_matching_ref = seq
                self.gene_matching_ref_type = seq_type
                self.gene_start_bases_added = start
                self.gene_end_bases_added = end

        class FakeCluster:
            def __init__(self, seq, seq_type, start, end):
                self.assembly_compare = FakeAssemblyCompare(seq, seq_type, start, end)

        self.clusters.clusters = {
            'gene1': FakeCluster(seq1, 'TYPE1', 1, 3),
            'gene2': FakeCluster(None, None, None, None),
            'gene3': FakeCluster(seq3, 'TYPE3', 4, 5),
        }

        tmp_file = 'tmp.test_write_catted_genes_matching_refs_fasta.fa'
        self.clusters._write_catted_genes_matching_refs_fasta(tmp_file)
        expected = os.path.join(data_dir, 'clusters_cat_genes_match_ref.fa')
        self.assertTrue(filecmp.cmp(expected, tmp_file, shallow=False))
        os.unlink(tmp_file)


    def test_write_mlst_reports(self):
        '''test _write_mlst_reports'''
        ariba_report = os.path.join(data_dir, 'clusters_test_write_mlst_reports.ariba.report.tsv')
        mlst_file = os.path.join(data_dir, 'clusters_test_write_mlst_reports.mlst_profile.tsv')
        expected_short = os.path.join(data_dir, 'clusters_test_write_mlst_reports.out.tsv')
        expected_long = os.path.join(data_dir, 'clusters_test_write_mlst_reports.out.details.tsv')
        outprefix = 'tmp.test_clusters__write_mlst_reports'
        got_short = outprefix + '.tsv'
        got_long = outprefix + '.details.tsv'

        clusters.Clusters._write_mlst_reports('filenotthere', ariba_report, outprefix)
        self.assertFalse(os.path.exists(got_short))
        self.assertFalse(os.path.exists(got_long))

        clusters.Clusters._write_mlst_reports(mlst_file, ariba_report, outprefix)
        self.assertTrue(filecmp.cmp(expected_short, got_short, shallow=False))
        self.assertTrue(filecmp.cmp(expected_long, got_long, shallow=False))
        os.unlink(got_short)
        os.unlink(got_long)


    def test_write_tb_resistance_calls_json(self):
        '''test _write_tb_resistance_calls_json'''
        ariba_report = os.path.join(data_dir, 'clusters_write_tb_resistance_calls_json.in.tsv')
        tmp_out = 'tmp.write_tb_resistance_calls_json.out.json'
        clusters.Clusters._write_tb_resistance_calls_json(ariba_report, tmp_out)
        expected = os.path.join(data_dir, 'clusters_write_tb_resistance_calls_json.out.json')
        self.assertTrue(filecmp.cmp(expected, tmp_out, shallow=False))
        os.unlink(tmp_out)


    def test_run_with_tb(self):
        '''test complete run with TB amr calling'''
        tmp_out = 'tmp.clusters_run_with_tb'
        c = clusters.Clusters(
            os.path.join(data_dir, 'clusters_run_with_tb.ref'),
            os.path.join(data_dir, 'clusters_run_with_tb.reads_1.fq.gz'),
            os.path.join(data_dir, 'clusters_run_with_tb.reads_2.fq.gz'),
            tmp_out,
            extern_progs,
        )
        c.run()

        expect_json = {"Pyrazinamide": [[ "pncA", "A3E"]]}
        json_file = os.path.join(tmp_out, 'tb.resistance.json')
        self.assertTrue(os.path.exists(json_file))
        with open(json_file) as f:
            got_json = json.load(f)
        self.assertEqual(expect_json, got_json)
        common.rmtree(tmp_out)
