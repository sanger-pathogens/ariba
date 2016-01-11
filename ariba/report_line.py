from ariba import flag

class Error (Exception): pass


columns = [
    'ref_name',
    'ref_type',
    'flag',
    'reads',
    'cluster',
    'gene_len',
    'assembled',
    'pc_ident',
    'known_snp_type',
    'known_snp',
    'var_type',
    'var_effect',
    'new_aa',
    'gene_start',
    'gene_end',
    'gene_nt',
    'scaffold',
    'scaff_len',
    'scaff_start',
    'scaff_end',
    'scaff_nt',
    'read_depth',
    'alt_bases',
    'ref_alt_depth',
    'sequence_description',
    'free_text',
]


class ReportLine:
    def __init__(self,
        ref_name,
        ref_data,
        flag,
        total_reads,
        cluster,
    ):
        self.data = {column: '.' for column in columns}
        self.ref_data = ref_data
        self.data['ref_name'] = ref_name
        self.data['ref_type'] = self.ref_data.sequence_type(ref_name)
        self.data['flag'] = flag
        self.data['total_reads'] = total_reads
        self.data['cluster'] = cluster
        self.data['gene_length'] = self.ref_data.sequence_length(ref_name)


    def __str__(self):
        return '\t'.join([str(self.data[column]) for column in columns])


