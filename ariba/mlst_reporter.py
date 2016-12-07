import os
import pyfastaq
from ariba import mlst_profile, summary_sample


class MlstReporter:
    def __init__(self, report_tsv, mlst_file, outprefix):
        self.summary_sample = summary_sample.SummarySample(report_tsv)
        self.mlst_profile = mlst_profile.MlstProfile(mlst_file)
        self.outprefix = outprefix
        self.allele_calls = {}
        self.any_allele_unsure = False


    def _call_gene(self, gene):
        results = {
            'allele': 'ND',
            'sure': None,
            'cov': '.',
            'ctgs': '.',
            'gene': None,
            'pc': '.',
            'depth': '.',
            'hetmax': '.',
        }

        if gene not in self.summary_sample.clusters:
            return results

        cluster = self.summary_sample.clusters[gene]

        best_hit = None
        for d in cluster.data:
            if best_hit is None or best_hit['ref_base_assembled'] < d['best_hit']:
                best_hit = d

        assert best_hit is not None
        results['cov'] = round(100.0 * best_hit['ref_base_assembled'] / best_hit['ref_len'], 2)
        results['ctgs'] = len({x['ctg'] for x in cluster.data})
        results['depth'] = best_hit['ctg_cov']
        results['pc'] = best_hit['pc_ident']
        results['allele'] = int(best_hit['ref_name'].split('.')[-1])

        return results


    def _call_genes(self):
        self.gene_results = {gene: self._call_gene(gene) for gene in self.mlst_profile.genes_list}


    def _call_sequence_type(self):
        type_dict = {gene: self.gene_results[gene]['allele'] for gene in self.mlst_profile.genes_list if self.gene_results[gene]['allele'] is not None}
        self.sequence_type = self.mlst_profile.get_sequence_type(type_dict)


    @classmethod
    def _report_strings(cls, results):
        allele_str = str(results['allele'])
        if results['sure'] is not None and not results['sure']:
            self.any_allele_unsure = True
            allele_str += '*'

        details_list = [x + ':' + str(results[x]) for x in ['cov', 'pc', 'ctgs', 'depth', 'hetmax']]

        if 'hets' in details_list:
            details_list.append('hets:' + ','.join([str(x) for x in details_list['hets']]))

        return allele_str, ';'.join(details_list)


    def _write_reports(self):
        f_out_simple = pyfastaq.utils.open_file_write(self.outprefix + '.tsv')
        f_out_all = pyfastaq.utils.open_file_write(self.outprefix + '.all.tsv')
        print('ST', *self.mlst_profile.genes_list, sep='\t', file=f_out_simple)
        print('ST', *[x + '\t' + x + '_details' for x in self.mlst_profile.genes_list], sep='\t', file=f_out_all)

        if self.sequence_type != 'ND' and self.any_allele_unsure:
            st_string = self.sequence_type + '*'
        else:
            st_string = self.sequence_type

        print(st_string, end='', file=f_out_simple)
        print(st_string, end='', file=f_out_all)

        for gene in self.mlst_profile.genes_list:
            allele_str, detail_str = MlstReporter._report_strings(self.gene_results[gene])
            print('\t', allele_str, sep='', end='', file=f_out_simple)
            print('\t', allele_str, '\t', detail_str, sep='', end='', file=f_out_all)

        print('', file=f_out_simple)
        print('', file=f_out_all)
        pyfastaq.utils.close(f_out_simple)
        pyfastaq.utils.close(f_out_all)


    def run(self):
        self.summary_sample.run()
        self._call_genes()
        self._call_sequence_type()
        self._write_reports()

