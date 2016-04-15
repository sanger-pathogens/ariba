import pyfastaq
import pprint

class Error (Exception): pass

class CardRecord:
    def __init__(self, data_dict):
        self.data_dict = data_dict


    @staticmethod
    def _ARO_id(data_dict):
        return data_dict.get('ARO_id', None)


    @staticmethod
    def _ARO_accession(data_dict):
        return data_dict.get('ARO_accession', None)


    @staticmethod
    def _ARO_name(data_dict):
        return data_dict.get('ARO_name', None)


    @staticmethod
    def _ARO_description(data_dict):
        return data_dict.get('ARO_description', None)


    @staticmethod
    def _dna_seqs_and_genbank_ids(gene_dict):
        try:
            seq_dict = gene_dict['model_sequences']['sequence']
        except:
            return []

        if len(seq_dict) == 0:
            return []

        seqs_and_ids = []

        for key, seq_dict in sorted(gene_dict['model_sequences']['sequence'].items()):
            try:
                dna_seq = seq_dict['dna_sequence']['sequence']
                genbank_id = seq_dict['dna_sequence']['accession']
                strand = seq_dict['dna_sequence']['strand']
                gi = seq_dict['protein_sequence']['GI']
            except:
                continue

            seqs_and_ids.append((key, gi, genbank_id, strand, dna_seq))

        return seqs_and_ids


    @staticmethod
    def _snps(data_dict):
        try:
            snps_dict = data_dict['model_param']['snp']
        except:
            return set()

        try:
            snps_set = set(snps_dict['param_value'].values())
        except:
            return set()

        return snps_set


    def get_data(self):
        data = {
            'ARO_id': self._ARO_id(self.data_dict),
            'ARO_accession': self._ARO_accession(self.data_dict),
            'ARO_name': self._ARO_name(self.data_dict),
            'ARO_description': self._ARO_description(self.data_dict),
            'dna_seqs_and_ids': self._dna_seqs_and_genbank_ids(self.data_dict),
        }

        if None in (data.values()) or data['dna_seqs_and_ids'] == 0:
            pprint(self.data_dict)
            raise Error('Error getting ARO_acccesion, ARO_id, ARO_name, dna_sequence(s) or genbank accession(s) from the above dictionary. Cannot continue.')

        data['snps'] = self._snps(self.data_dict)
        return data

