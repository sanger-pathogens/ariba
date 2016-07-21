import os
import pickle
from ariba import reference_data


class Error (Exception): pass

class RefdataQuery:
    def __init__(self,
        prepareref_dir
    ):
        if not os.path.exists(prepareref_dir):
            raise Error('Error querying refdata. Input directory "' + prepareref_dir + '" not found.')

        self.prepareref_dir = prepareref_dir
        self.refdata_fa = os.path.join(self.prepareref_dir, '02.cdhit.all.fa')
        self.refdata_tsv = os.path.join(self.prepareref_dir, '01.filter.check_metadata.tsv')
        self.clusters_pickle = os.path.join(self.prepareref_dir, '02.cdhit.clusters.pickle')


    @staticmethod
    def _load_clusters(pickle_file):
        with open(pickle_file, 'rb') as f:
            clusters = pickle.load(f)
        return clusters


    @staticmethod
    def _seq2cluster(clusters, seqname):
        for cluster in clusters:
            if seqname in clusters[cluster]:
                return cluster
        return None


    def _cluster2seqs(self, cluster_name):
        clusters = self._load_clusters(self.clusters_pickle)
        if cluster_name not in clusters:
            return ['Cluster name "' + cluster_name + '" not found']
        else:
            return ['Sequences belonging to cluster ' + cluster_name + ':'] + sorted(list(clusters[cluster_name]))


    def _seqinfo(self, seqname):
        refdata = reference_data.ReferenceData([self.refdata_fa], [self.refdata_tsv])
        if seqname not in refdata.sequences:
            return ['Sequence "' + seqname + '" not found']

        assert seqname in refdata.metadata

        clusters = RefdataQuery._load_clusters(self.clusters_pickle)
        cluster = RefdataQuery._seq2cluster(clusters, seqname)
        assert cluster is not None

        if refdata.metadata[seqname]['seq_type'] == 'p':
            is_gene = '1'
            var_dict = refdata.metadata[seqname]['p']
        else:
            is_gene = '0'
            var_dict = refdata.metadata[seqname]['n']

        description_lines = ['Description\t' + x.free_text for x in sorted(list(refdata.metadata[seqname]['.'])) if x.free_text != '.']

        var_lines = []
        for position in sorted(var_dict):
            for seq_meta in sorted(var_dict[position]):
                ident = '.' if seq_meta.variant.identifier is None else seq_meta.variant.identifier
                var_lines.append('\t'.join(['Variant', str(seq_meta.variant), ident, seq_meta.free_text]))

        return [
            'Name\t' + seqname,
            'Is gene\t' + is_gene,
            'Variant only\t' + ('1' if refdata.metadata[seqname]['variant_only'] else '0'),
            'Cluster\t' + cluster,
        ] + description_lines + var_lines + ['Sequence\t' + refdata.sequences[seqname].seq]


    def query(self, query_type, query_string):
        queries = {
            'cluster': self._cluster2seqs,
            'seq': self._seqinfo,
        }

        if query_type not in queries:
            raise Error('Unknown query type "' + query_type + '". Choices are:\n' + ','.join(sorted(queries.keys())) + '\n')

        lines = queries[query_type](query_string)
        print(*lines, sep='\n')

