import pyfastaq
from ariba import mlst_profile, summary_sample


class MlstReporter:
    def __init__(self, report_tsv, mlst_file, outprefix):
        self.summary_sample = summary_sample.SummarySample(report_tsv)
        self.mlst_profile = mlst_profile.MlstProfile(mlst_file, duplicate_warnings=False)
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
            'hetmin': '.',
            'hets': '.',
        }

        if gene not in self.summary_sample.clusters:
            return results

        cluster = self.summary_sample.clusters[gene]

        best_hit = None

        for d in cluster.data:
            if best_hit is None or best_hit['ref_base_assembled'] < d['ref_base_assembled']:
                best_hit = d

        assert best_hit is not None

        het_data = []

        for d in cluster.data:
            if d['ctg'] == best_hit['ctg']:
                if d['var_type'] == 'HET':
                    het_data.append(d['smtls_nts_depth'])
                    depths = [int(x) for x in d['smtls_nts_depth'].split(',')]
                    depths.sort()
                    het_pc = round(100.0 * depths[-1] / sum(depths), 2)
                    if results['hetmin'] == '.' or results['hetmin'] > het_pc:
                        results['hetmin'] = het_pc
        if len(het_data):
            results['hets'] = '.'.join(het_data)

        results['cov'] = round(100.0 * best_hit['ref_base_assembled'] / best_hit['ref_len'], 2)
        results['ctgs'] = len({x['ctg'] for x in cluster.data})
        results['depth'] = best_hit['ctg_cov']
        results['pc'] = best_hit['pc_ident']
        results['allele'] = int(best_hit['ref_name'].split('.')[-1])
        results['sure'] = True

        if results['hetmin'] != '.' or results['cov'] < 100 or results['ctgs'] != 1 or results['pc'] < 100:
            results['sure'] = False

        return results


    def _call_genes(self):
        self.gene_results = {gene: self._call_gene(gene) for gene in self.mlst_profile.genes_list}


    def _call_sequence_type(self):
        type_dict = {gene: self.gene_results[gene]['allele'] for gene in self.mlst_profile.genes_list if self.gene_results[gene]['allele'] is not None}
        self.sequence_type = self.mlst_profile.get_sequence_type(type_dict)


    def _report_strings(self, results):
        allele_str = str(results['allele'])
        if results['sure'] is not None and not results['sure']:
            self.any_allele_unsure = True
            allele_str += '*'

        details_list = [str(results[x]) for x in ['cov', 'pc', 'ctgs', 'depth', 'hetmin', 'hets']]
        return allele_str, '\t'.join(details_list)


    def _write_reports(self):
        f_out_all = pyfastaq.utils.open_file_write(self.outprefix + '.details.tsv')
        print('gene\tallele\tcov\tpc\tctgs\tdepth\thetmin\thets', file=f_out_all)
        out_simple = [str(self.sequence_type)]

        for gene in self.mlst_profile.genes_list:
            allele_str, detail_str = self._report_strings(self.gene_results[gene])
            out_simple.append(allele_str)
            print(gene, allele_str, detail_str, sep='\t', file=f_out_all)

        pyfastaq.utils.close(f_out_all)

        if self.sequence_type != 'ND' and self.any_allele_unsure:
            out_simple[0] += '*'

        f_out_simple = pyfastaq.utils.open_file_write(self.outprefix + '.tsv')
        print('ST', *self.mlst_profile.genes_list, sep='\t', file=f_out_simple)
        print(*out_simple, sep='\t', file=f_out_simple)
        pyfastaq.utils.close(f_out_simple)


    def run(self):
        self.summary_sample.run()
        self._call_genes()
        self._call_sequence_type()
        self._write_reports()

