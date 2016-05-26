import pyfastaq
from ariba import report, summary_cluster

class Error (Exception): pass

class SummarySample:
    def __init__(self, report_tsv, min_pc_id=90):
        self.report_tsv = report_tsv
        self.min_pc_id = min_pc_id
        self.clusters = {}


    def __eq__(self, other):
       return type(other) is type(self) and self.__dict__ == other.__dict__


    @staticmethod
    def _load_file(filename, min_pc_id):
        f = pyfastaq.utils.open_file_read(filename)
        clusters = {}

        for line in f:
            if line.startswith('#'):
                if line.rstrip()[1:].split('\t') != report.columns:
                    pyfastaq.utils.close(f)
                    raise Error('Error parsing the following line.\n' + line)
                continue

            data_dict = summary_cluster.SummaryCluster.line2dict(line)
            cluster = data_dict['cluster']
            if cluster not in clusters:
                clusters[cluster] = summary_cluster.SummaryCluster(min_pc_id=min_pc_id)
            clusters[cluster].add_data_dict(data_dict)

        pyfastaq.utils.close(f)
        return clusters


    def _column_summary_data(self):
        return {c: self.clusters[c].column_summary_data() for c in self.clusters}


    def _var_groups(self):
        return {c: self.clusters[c].has_var_groups() for c in self.clusters}


    def _variant_column_names_tuples(self):
        variants = {}
        for cluster_name, cluster in self.clusters.items():
            cluster_vars = cluster.non_synon_variants()
            if len(cluster_vars):
                variants[cluster_name] = cluster_vars
        return variants


    def run(self):
        self.clusters = self._load_file(self.report_tsv, self.min_pc_id)
        self.column_summary_data = self._column_summary_data()
        self.variant_column_names_tuples = self._variant_column_names_tuples()
        self.var_groups = self._var_groups()

