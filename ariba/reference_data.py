import os
import sys
import copy
import pyfastaq
from ariba import sequence_metadata


class Error (Exception): pass


class ReferenceData:
    def __init__(self,
        presence_absence_fa=None,
        variants_only_fa=None,
        non_coding_fa=None,
        metadata_tsv=None,
        min_gene_length=6,
        max_gene_length=10000,
        genetic_code=11,
    ):
        self.seq_filenames = {}
        self.seq_dicts = {}
        self.min_gene_length = min_gene_length
        self.max_gene_length = max_gene_length


        for x in ['presence_absence', 'variants_only', 'non_coding']:
            exec('self.seq_filenames[x] = self._get_filename(' + x + '_fa)')
            self.seq_dicts[x] = self._load_fasta_file(self.seq_filenames[x])

        if {None} == set(self.seq_filenames.values()):
            raise Error('Error! Must supply at least one of presence_absence_fa, variants_only_fa, non_coding_fa. Cannot continue')

        self.metadata = self._load_metadata_tsv(metadata_tsv)
        self.genetic_code = genetic_code
        pyfastaq.sequences.genetic_code = self.genetic_code
        common_names = self._dict_keys_intersection(list(self.seq_dicts.values()))
        if len(common_names):
            raise Error('Error! Non-unique names found in input fasta files:\n' + '\n'.join(common_names))


    @staticmethod
    def _dict_keys_intersection(dicts):
        dicts = [x for x in dicts if x is not None]
        if len(dicts) == 0:
            return set()

        inter = set(dicts[0].keys())

        for d in dicts[1:]:
            inter = inter.intersection(set(d.keys()))
        return inter


    @staticmethod
    def _get_filename(filename):
        if filename is None:
            return None
        else:
            if os.path.exists(filename):
                return os.path.abspath(filename)
            else:
                raise Error('Error! File not found: ' + filename)


    @staticmethod
    def _load_metadata_tsv(filename):
        if filename is None:
            return {}

        f = pyfastaq.utils.open_file_read(filename)
        metadata_dict = {}

        for line in f:
            try:
                metadata = sequence_metadata.SequenceMetadata(line)
            except:
                print('Problem with this line of metadata, which will be ignored:', line.rstrip(), file=sys.stderr)

            if metadata.name not in metadata_dict:
                metadata_dict[metadata.name] = []

            metadata_dict[metadata.name].append(metadata)

        pyfastaq.utils.close(f)
        return metadata_dict


    @staticmethod
    def _load_fasta_file(filename):
        d = {}

        if filename is not None:
            seq_reader = pyfastaq.sequences.file_reader(filename)
            for seq in seq_reader:
                seq.id = seq.id.split()[0]
                if seq.id in d:
                    raise Error('Duplicate name "' + seq.id + '" found in file ' + filename + '. Cannot continue)')
                d[seq.id] = copy.copy(seq)

        return d


    @staticmethod
    def _find_gene_in_seqs(name, dicts):
        for dict_name, this_dict in dicts.items():
            if this_dict is None:
                continue
            elif name in this_dict:
                return dict_name

        return None


    @staticmethod
    def _write_metadata_tsv(metadata, filename):
        f = pyfastaq.utils.open_file_write(filename)

        for gene_name, variants in sorted(metadata.items()):
            for variant in variants:
                print(variant, sep='\t', file=f)

        pyfastaq.utils.close(f)


    @staticmethod
    def _write_dict_of_sequences(seq_dict, filename):
        f = pyfastaq.utils.open_file_write(filename)
        for seq in sorted(seq_dict):
            print(seq_dict[seq], file=f)
        pyfastaq.utils.close(f)


    def _write_sequences(self, filename, sequences_to_write):
        assert sequences_to_write in self.seq_dicts and sequences_to_write in self.seq_filenames
        if self.seq_filenames[sequences_to_write] is not None:
            self._write_dict_of_sequences(self.seq_dicts[sequences_to_write], filename)


    def _filter_bad_variant_data(self, out_prefix, presence_absence_removed, variants_only_removed):
        genes_to_remove = set()
        variants_only_genes_not_found = set(self.seq_dicts['variants_only'].keys())
        log_file = out_prefix + '.log'
        tsv_file = out_prefix + '.tsv'
        new_variants_fa_file = out_prefix + '.variants_only.fa'
        log_fh = pyfastaq.utils.open_file_write(log_file)

        for gene_name, metadata_list in sorted(self.metadata.items()):
            if gene_name in presence_absence_removed:
                for metadata in metadata_list:
                    print(gene_name, 'was removed from presence/absence fasta, so removing its metadata. Line of tsv was:', metadata, file=log_fh)
                genes_to_remove.add(gene_name)
                continue
            elif gene_name in variants_only_removed:
                for metadata in metadata_list:
                    print(gene_name, 'was removed from variants only fasta, so removing its metadata. Line of tsv was:', metadata, file=log_fh)
                genes_to_remove.add(gene_name)
                continue

            gene_in_seq_dict = self._find_gene_in_seqs(gene_name, self.seq_dicts)
            if gene_in_seq_dict is None:
                for metadata in metadata_list:
                    print(gene_name, 'is in input tsv file, but not found in any input sequence files. Removing. Line of tsv file was:', metadata, file=log_fh)
                genes_to_remove.add(gene_name)
                continue

            to_remove = []

            for i in range(len(metadata_list)):
                metadata = metadata_list[i]
                if metadata.variant_type == '.':
                    if metadata.free_text is None or gene_in_seq_dict == 'variants_only':
                        print(gene_name, 'metadata has no info. Just gene name given. Removing. Line of file was:', metadata, file=log_fh)
                        to_remove.append(i)
                    continue

                if gene_in_seq_dict == 'non_coding' and metadata.variant_type == 'p':
                    print(gene_name, 'variant of type "p" for protein, but sequence is non-coding. Removing. Line of file was:', metadata, file=log_fh)
                    to_remove.append(i)
                    continue

                to_translate = metadata.variant_type == 'p'
                if not metadata.variant.sanity_check_against_seq(self.seq_dicts[gene_in_seq_dict][gene_name], translate_seq=to_translate):
                    print(gene_name, 'variant does not match reference. Removing. Line of file was:', metadata, file=log_fh)
                    to_remove.append(i)
                    continue

                if gene_in_seq_dict == 'variants_only':
                    variants_only_genes_not_found.discard(gene_name)

            for i in sorted(to_remove, reverse=True):
                metadata_list.pop(i)

            if len(metadata_list) == 0 and gene_in_seq_dict == 'variants_only':
                print(gene_name, 'No remaining data after checks. Removing this sequence because it is in the variants only file', file=log_fh)
                genes_to_remove.add(gene_name)

        for gene_name in genes_to_remove:
            self.metadata.pop(gene_name)

        for gene_name in variants_only_genes_not_found:
            print(gene_name, 'is in variants only gene file, but no variants found. Removing.', file=log_fh)
            self.seq_dicts['variants_only'].pop(gene_name)

        pyfastaq.utils.close(log_fh)
        self._write_metadata_tsv(self.metadata, tsv_file)
        self._write_sequences(out_prefix + '.presence_absence.fa', 'presence_absence')
        self._write_sequences(out_prefix + '.non_coding.fa', 'non_coding')
        self._write_sequences(out_prefix + '.variants_only.fa', 'variants_only')


    @staticmethod
    def _gene_seq_is_ok(seq, min_length, max_length, genetic_code):
        seq.seq = seq.seq.upper()
        if len(seq) < min_length:
            return False, 'Remove: too short. Length: ' + str(len(seq))
        elif len(seq) > max_length:
            return False, 'Remove: too long. Length: ' + str(len(seq))
        elif not seq.looks_like_gene(translation_table=genetic_code):
            length_over_three = round(len(seq) / 3, 2)
            return False, 'Does not look like a gene (does not start with start codon, length (' + str(len(seq)) + ') is not a multiple of 3 (length/3=' + str(length_over_three) + '), or contains internal stop codons). Translation: ' + seq.translate().seq

        return True, None


    def _remove_bad_genes(self, seqs_dict, log_file):
        to_remove = set()

        if len(seqs_dict) == 0:
            return to_remove

        log_fh = pyfastaq.utils.open_file_write(log_file)

        for name, sequence in sorted(seqs_dict.items()):
            ok, message = self._gene_seq_is_ok(sequence, self.min_gene_length, self.max_gene_length, self.genetic_code)
            if message is not None:
                print(name, message, file=log_fh)
            if not ok:
                to_remove.add(name)

        pyfastaq.utils.close(log_fh)

        for name in to_remove:
            seqs_dict.pop(name)

        return to_remove


    def sanity_check(self, outprefix):
        variants_only_removed = self._remove_bad_genes(self.seq_dicts['variants_only'], outprefix + '.00.check_fasta_variants_only.log')
        presence_absence_removed = self._remove_bad_genes(self.seq_dicts['presence_absence'], outprefix + '.00.check_fasta_presence_absence.log')
        self._filter_bad_variant_data(outprefix + '.01.check_variants', variants_only_removed, presence_absence_removed)


    def make_catted_fasta(self, outfile):
        f = pyfastaq.utils.open_file_write(outfile)

        for key in ['presence_absence', 'variants_only', 'non_coding']:
            filename = self.seq_filenames[key]
            if filename is not None:
                file_reader = pyfastaq.sequences.file_reader(filename)
                for seq in file_reader:
                    print(seq, file=f)

        pyfastaq.utils.close(f)


    def sequence_type(self, sequence_name):
        return self._find_gene_in_seqs(sequence_name, self.seq_dicts)


    def sequence(self, sequence_name):
        d = self._find_gene_in_seqs(sequence_name, self.seq_dicts)
        if d is None:
            return None
        else:
            return self.seq_dicts[d][sequence_name]


    def sequence_length(self, sequence_name):
        seq = self.sequence(sequence_name)
        assert seq is not None
        return len(seq)
