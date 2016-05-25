import os
import sys
import re
import copy
import pyfastaq
from ariba import sequence_metadata, cdhit


class Error (Exception): pass

rename_sub_regex = re.compile(r'[^\w.-]')


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

        total_ref_seqs_loaded = 0

        for x in ['presence_absence', 'variants_only', 'non_coding']:
            exec('self.seq_filenames[x] = self._get_filename(' + x + '_fa)')
            self.seq_dicts[x] = self._load_fasta_file(self.seq_filenames[x])
            total_ref_seqs_loaded += len(self.seq_dicts[x])

        if {None} == set(self.seq_filenames.values()):
            raise Error('Error! Must supply at least one of presence_absence_fa, variants_only_fa, non_coding_fa. Cannot continue')

        if total_ref_seqs_loaded == 0:
            raise Error('Error! No sequences found in input file(s). Maybe they were empty? Cannot continue.')

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
                continue

            if metadata.name not in metadata_dict:
                metadata_dict[metadata.name] = {'n': {}, 'p': {}, '.': set()}

            if metadata.variant_type == '.':
                metadata_dict[metadata.name]['.'].add(metadata)
            else:
                if metadata.variant.position not in metadata_dict[metadata.name][metadata.variant_type]:
                    metadata_dict[metadata.name][metadata.variant_type][metadata.variant.position] = set()

                metadata_dict[metadata.name][metadata.variant_type][metadata.variant.position].add(metadata)

        pyfastaq.utils.close(f)
        return metadata_dict


    @staticmethod
    def _load_fasta_file(filename):
        d = {}

        if filename is not None:
            seq_reader = pyfastaq.sequences.file_reader(filename)
            for seq in seq_reader:
                #seq.id = seq.id.split()[0]
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

        for gene_name, data_dict in sorted(metadata.items()):
            for meta in sorted([str(x) for x in data_dict['.']]):
                print(meta, file=f)

            variants = []

            for variant_type in ['n', 'p']:
                for position in data_dict[variant_type]:
                    for meta in data_dict[variant_type][position]:
                        variants.append(meta)

            variants.sort()
            for v in variants:
                print(v, file=f)

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
        log_fh = pyfastaq.utils.open_file_write(log_file)

        for gene_name, metadata_dict in sorted(self.metadata.items()):
            if gene_name in presence_absence_removed:
                print(gene_name, 'was removed from presence/absence fasta, so removing its metadata', file=log_fh)
                genes_to_remove.add(gene_name)
                continue
            elif gene_name in variants_only_removed:
                print(gene_name, 'was removed from variants only fasta, so removing its metadata', file=log_fh)
                genes_to_remove.add(gene_name)
                continue

            gene_in_seq_dict = self._find_gene_in_seqs(gene_name, self.seq_dicts)
            if gene_in_seq_dict is None:
                print(gene_name, 'is in input tsv file, but not found in any input sequence files. Removing', file=log_fh)
                genes_to_remove.add(gene_name)
                continue

            # take out any metadata that is not a variant and has no extra info.
            to_remove = []

            for metadata in metadata_dict['.']:
                if metadata.free_text == '.':
                    print(gene_name, 'metadata has no info. Just gene name given. Removing. Line of file was:', metadata, file=log_fh)
                    to_remove.append(metadata)

            for metadata in to_remove:
                metadata_dict['.'].remove(metadata)


            # if this is non_coding, we shouldn't have any amino acid variants
            if gene_in_seq_dict == 'non_coding':
                for position in metadata_dict['p']:
                    for metadata in metadata_dict['p'][position]:
                        print(gene_name, 'variant of type "p" for protein, but sequence is non-coding. Removing. Line of file was:', metadata, file=log_fh)

                metadata_dict['p'] = {}


            # take out variant metadata that doesn't make sense (eg bases not matching ref sequence)
            for variant_type in ['n', 'p']:
                positions_to_remove = []
                for position in metadata_dict[variant_type]:
                    meta_to_remove = []
                    for metadata in metadata_dict[variant_type][position]:
                        to_translate = variant_type == 'p'

                        if not metadata.variant.sanity_check_against_seq(self.seq_dicts[gene_in_seq_dict][gene_name], translate_seq=to_translate):
                            print(gene_name, 'variant does not match reference. Removing. Line of file was:', metadata, file=log_fh)
                            meta_to_remove.append(metadata)
                            continue

                        if gene_in_seq_dict == 'variants_only':
                            variants_only_genes_not_found.discard(gene_name)

                    for metadata in meta_to_remove:
                        metadata_dict[variant_type][position].remove(metadata)
                    if len(metadata_dict[variant_type][position]) == 0:
                        positions_to_remove.append(position)

                for position in positions_to_remove:
                    del metadata_dict[variant_type][position]


            if gene_in_seq_dict == 'variants_only' and len(metadata_dict['n']) == len(metadata_dict['p']) == len(metadata_dict['.']) == 0:
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
    def _try_to_get_gene_seq(seq, min_length, max_length):
        seq.seq = seq.seq.upper()
        if len(seq) < min_length:
            return None, 'Remove: too short. Length: ' + str(len(seq))
        elif len(seq) > max_length:
            return None, 'Remove: too long. Length: ' + str(len(seq))
        else:
            got = seq.make_into_gene()
            if got is None:
                return None, 'Does not look like a gene (tried both strands and all reading frames) ' + seq.seq
            else:
                return got[0], 'Made ' + seq.id + ' into gene. strand=' + got[1] + ', frame=' + str(got[2])


    def _remove_bad_genes(self, seqs_dict, log_file):
        to_remove = set()

        if len(seqs_dict) == 0:
            return to_remove

        log_fh = pyfastaq.utils.open_file_write(log_file)

        for name in sorted(seqs_dict):
            new_seq, message = self._try_to_get_gene_seq(seqs_dict[name], self.min_gene_length, self.max_gene_length)
            if new_seq is None:
                to_remove.add(name)
            else:
                seqs_dict[name] = new_seq

            if message is not None:
                print(name, message, file=log_fh)

        pyfastaq.utils.close(log_fh)

        for name in to_remove:
            seqs_dict.pop(name)

        return to_remove


    def sanity_check(self, outprefix):
        variants_only_removed = self._remove_bad_genes(self.seq_dicts['variants_only'], outprefix + '.00.check_fasta_variants_only.log')
        presence_absence_removed = self._remove_bad_genes(self.seq_dicts['presence_absence'], outprefix + '.00.check_fasta_presence_absence.log')
        self._filter_bad_variant_data(outprefix + '.01.check_variants', presence_absence_removed, variants_only_removed)


    @classmethod
    def _new_seq_name(cls, name):
        name = name.split()[0]
        return re.sub(rename_sub_regex, '_', name)


    @classmethod
    def _seq_names_to_rename_dict(cls, names):
        used_names = set()
        old_name_to_new = {}

        for old_name in sorted(names):
            new_name = ReferenceData._new_seq_name(old_name)
            if new_name in used_names:
                i = 1
                new_name_prefix = new_name
                while new_name in used_names:
                    new_name = new_name_prefix + '_' + str(i)
                    i += 1

            assert new_name not in used_names
            if new_name != old_name:
                old_name_to_new[old_name] = new_name

            used_names.add(new_name)

        return old_name_to_new


    @classmethod
    def _rename_names_in_seq_dicts(cls, seq_dicts, rename_dict):
        '''Changes seq_dicts in place'''
        for seq_type in ['presence_absence', 'variants_only', 'non_coding']:
            new_dict = {}
            while len(seq_dicts[seq_type]):
                old_name, seq = seq_dicts[seq_type].popitem()
                if old_name in rename_dict:
                    seq.id = rename_dict[old_name]

                new_dict[seq.id] = seq
            seq_dicts[seq_type] = new_dict


    @classmethod
    def _rename_metadata_set(cls, metadata_set, new_name):
        new_set = set()
        for meta in metadata_set:
            new_meta = copy.copy(meta)
            new_meta.name = new_name
            new_set.add(new_meta)
        return new_set


    @classmethod
    def _rename_names_in_metadata(cls, meta_dict, rename_dict):
        new_dict = {}

        while len(meta_dict):
            old_name, gene_dict = meta_dict.popitem()
            if old_name in rename_dict:
                new_name = rename_dict[old_name]
                for seq_type in ['n', 'p']:
                    for position, metaset in gene_dict[seq_type].items():
                        gene_dict[seq_type][position] = ReferenceData._rename_metadata_set(metaset, new_name)

                gene_dict['.'] = ReferenceData._rename_metadata_set(gene_dict['.'], new_name)
            else:
                new_name = old_name

            new_dict[new_name] = gene_dict

        return new_dict


    def rename_sequences(self, outfile):
        presabs_names = set(self.seq_dicts['presence_absence'].keys())
        noncoding_names = set(self.seq_dicts['non_coding'].keys())
        varonly_names = set(self.seq_dicts['variants_only'].keys())
        # we should have already checked that all the names are unique, but let's do it again!
        all_names = presabs_names.union(noncoding_names).union(varonly_names)
        if len(all_names) != len(presabs_names) + len(noncoding_names) + len(varonly_names):
            raise Error('Got a non-unique name in input data. Cannot continue')

        rename_dict = ReferenceData._seq_names_to_rename_dict(all_names)
        if len(rename_dict):
            print('Had to rename some sequences. See', outfile, 'for old -> new names', file=sys.stderr)
            with open(outfile, 'w') as f:
                for old_name, new_name in sorted(rename_dict.items()):
                    print(old_name, new_name, sep='\t', file=f)

            ReferenceData._rename_names_in_seq_dicts(self.seq_dicts, rename_dict)
            self.metadata = ReferenceData._rename_names_in_metadata(self.metadata, rename_dict)


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


    def all_non_wild_type_variants(self, ref_name):
        ref_seq = self.sequence(ref_name)
        variants = {'n': {}, 'p': {}}

        if ref_seq is None or ref_name not in self.metadata:
            return variants

        for variant_type in ['n', 'p']:
            for position, metadata_set in self.metadata[ref_name][variant_type].items():
                for metadata in metadata_set:
                    if position not in variants[variant_type]:
                        variants[variant_type][position] = set()

                    variants[variant_type][position].add(metadata)

        return variants


    @staticmethod
    def write_cluster_allocation_file(clusters, outfile):
        f_out = pyfastaq.utils.open_file_write(outfile)

        for seq_type in ['presence_absence', 'variants_only', 'non_coding']:
            if clusters[seq_type] is None:
                continue

            for seq_name in sorted(clusters[seq_type]):
                other_seqs = clusters[seq_type][seq_name].difference({seq_name})
                if len(other_seqs) > 0:
                    other_seq_string = '\t'.join(sorted(list(other_seqs)))
                    print(seq_name, other_seq_string, sep='\t', file=f_out)
                else:
                    print(seq_name, file=f_out)

        pyfastaq.utils.close(f_out)


    def cluster_with_cdhit(self, inprefix, outprefix, seq_identity_threshold=0.9, threads=1, length_diff_cutoff=0.9, nocluster=False, verbose=False, cd_hit_est='cd-hit-est', clusters_file=None):
        files_to_cat = []
        clusters = {}

        for seqs_type in ['presence_absence', 'variants_only', 'non_coding']:
            if len(self.seq_dicts[seqs_type]) > 0:
                outfile = outprefix + '.' + seqs_type + '.cdhit'
                files_to_cat.append(outfile)
                cdhit_runner = cdhit.Runner(
                  inprefix + '.' + seqs_type + '.fa',
                  outfile,
                  seq_identity_threshold=seq_identity_threshold,
                  threads=threads,
                  length_diff_cutoff=length_diff_cutoff,
                  verbose=verbose,
                  cd_hit_est=cd_hit_est,
                  rename_suffix = seqs_type[0]
                )

                if clusters_file is not None:
                    new_clusters = cdhit_runner.run_get_clusters_from_file(clusters_file)
                elif nocluster:
                    new_clusters = cdhit_runner.fake_run()
                else:
                    new_clusters = cdhit_runner.run()

                clusters[seqs_type] = new_clusters
            else:
                clusters[seqs_type] = None

        assert len(files_to_cat) > 0
        f_out = pyfastaq.utils.open_file_write(outprefix + '.cluster_representatives.fa')

        for filename in files_to_cat:
            for seq in pyfastaq.sequences.file_reader(filename):
                print(seq, file=f_out)

        pyfastaq.utils.close(f_out)
        self.write_cluster_allocation_file(clusters, outprefix + '.clusters.tsv')
        return clusters


    def write_seqs_to_fasta(self, outfile, names):
        f_out = pyfastaq.utils.open_file_write(outfile)

        for name in sorted(names):
            print(self.sequence(name), file=f_out)

        pyfastaq.utils.close(f_out)
