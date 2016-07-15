import os
import sys
import pyfastaq
from ariba import reference_data


class Error (Exception): pass

class RefdataQuery:
    def __init__(self,
        prepareref_dir
    ):
        if not os.path.exists(prepareref_dir):
            raise Error('Error querying refdata. Input directory "' + prepareref_dir + '" not found.')

        self.prepareref_dir = prepareref_dir


    def _cluster2seqs(self, cluster_name):
        lines = ['Sequences belonging to cluster ' + cluster_name + ':', '1', '2']
        return lines


    def query(self, query_type, query_string):
        queries = {
            'cluster2seqs': self._cluster2seqs,
        }

        if query_type not in queries:
            raise Error('Unknown query type "' + query_type + '". Choices are:\n' + ','.join(sorted(queries.keys())) + '\n')

        lines = queries[query_type](query_string)
        print(*lines, sep='\n')

