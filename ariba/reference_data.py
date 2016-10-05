import os
import sys
import re
import copy
import pyfastaq
from ariba import sequence_metadata, cdhit


class Error (Exception): pass

rename_sub_regex = re.compile(r'''[^a-zA-Z0-9_.]''')


class ReferenceData:
    def __init__(self,
        fasta_files,
        metadata_tsv_files,
        rename_file=None,
        min_gene_length=6,
        max_gene_length=10000,
        genetic_code=11,
    ):
        self.seq_filenames = {}
        self.seq_dicts = {}
        self.min_gene_length = min_gene_length
        self.max_gene_length = max_gene_length

        self.sequences, self.metadata = ReferenceData._load_input_files_and_check_seq_names(fasta_files, metadata_tsv_files)
        if len(self.sequences) == 0:
            raise Error('Error. No sequences found in input file(s):' + '\n'.join(fasta_files) + '\nCannot continue')

        self.genetic_code = genetic_code
        pyfastaq.sequences.genetic_code = self.genetic_code
        self.rename_dict = None

        if rename_file is None or not os.path.exists(rename_file):
            self.ariba_to_original_name = {}
        else:
            self.ariba_to_original_name = ReferenceData._load_rename_file(rename_file)


    @classmethod
    def _load_rename_file(cls, filename):
        ariba_name_to_original_name = {}
        f = pyfastaq.utils.open_file_read(filename)
        for line in f:
            original_name, ariba_name = line.rstrip().split('\t')
            ariba_name_to_original_name[ariba_name] = original_name
        pyfastaq.utils.close(f)
        return ariba_name_to_original_name


    @classmethod
    def _load_metadata_tsv(cls, filename, metadata_dict):
        if filename is None:
            return {}

        f = pyfastaq.utils.open_file_read(filename)

        for line in f:
            try:
                metadata = sequence_metadata.SequenceMetadata(line)
            except:
                print('Problem with this line of metadata, which will be ignored:', line.rstrip(), file=sys.stderr)
                continue

            if metadata.name not in metadata_dict:
                metadata_dict[metadata.name] = {
                    'seq_type': metadata.seq_type,
                    'variant_only': metadata.variant_only,
                    'n': {},
                    'p': {},
                    '.': set()
                }
            elif metadata.seq_type != metadata_dict[metadata.name]['seq_type'] or metadata.variant_only != metadata_dict[metadata.name]['variant_only']:
                raise Error('Inconsistent metadata for sequence ' + metadata.name + '. Cannot continue')

            if metadata.variant is None:
                metadata_dict[metadata.name]['.'].add(metadata)
            else:
                if metadata.variant.position not in metadata_dict[metadata.name][metadata.seq_type]:
                    metadata_dict[metadata.name][metadata.seq_type][metadata.variant.position] = set()

                metadata_dict[metadata.name][metadata.seq_type][metadata.variant.position].add(metadata)

        pyfastaq.utils.close(f)
        return metadata_dict


    @classmethod
    def _load_all_metadata_tsvs(cls, filenames):
        metadata_dict = {}
        assert len(filenames) > 0
        for filename in filenames:
            ReferenceData._load_metadata_tsv(filename, metadata_dict)
        return metadata_dict


    @classmethod
    def _load_fasta_file(cls, filename, seq_dict):
        if filename is not None:
            seq_reader = pyfastaq.sequences.file_reader(filename)
            for seq in seq_reader:
                seq.id = seq.id.split()[0]

                if seq.id in seq_dict:
                    raise Error('Duplicate name "' + seq.id + '" found in file ' + filename + '. Cannot continue)')
                seq_dict[seq.id] = copy.copy(seq)


    @classmethod
    def _load_all_fasta_files(cls, filenames):
        seq_dict = {}
        assert len(filenames) > 0
        for filename in filenames:
            ReferenceData._load_fasta_file(filename, seq_dict)
        return seq_dict



    @staticmethod
    def _load_input_files_and_check_seq_names(fasta_files, metadata_files):
        metadata = ReferenceData._load_all_metadata_tsvs(metadata_files)
        all_seqs = ReferenceData._load_all_fasta_files(fasta_files)

        for seq_name in all_seqs:
            if seq_name not in metadata:
                raise Error('Sequence "' + seq_name + '" found in input fasta file but not in metadata file. Cannot continue')

        to_remove = set()

        for seq_name in metadata:
            if seq_name not in all_seqs:
                to_remove.add(seq_name)
                print('WARNING: sequence "', seq_name, '" found in metadata, but not in fasta file. Ignoring it.', sep='', file=sys.stderr)

        for key in to_remove:
            del metadata[key]

        return all_seqs, metadata


    @classmethod
    def _write_metadata_tsv(cls, metadata, filename):
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


    @classmethod
    def _write_sequences_to_files(cls, sequences, metadata, outprefix):
        filenames = {
            ('n', False): outprefix + '.noncoding.fa',
            ('n', True): outprefix + '.noncoding.varonly.fa',
            ('p', False): outprefix + '.gene.fa',
            ('p', True): outprefix + '.gene.varonly.fa',
        }

        filename2filehandle = {}
        for key in filenames:
            filename2filehandle[filenames[key]] = pyfastaq.utils.open_file_write(filenames[key])

        all_fh = pyfastaq.utils.open_file_write(outprefix + '.all.fa')

        for sequence in sorted(sequences):
            key = metadata[sequence]['seq_type'], metadata[sequence]['variant_only']
            filehandle = filename2filehandle[filenames[key]]
            print(sequences[sequence], file=filehandle)
            print(sequences[sequence], file=all_fh)

        for filehandle in filename2filehandle.values():
            pyfastaq.utils.close(filehandle)

        pyfastaq.utils.close(all_fh)


    @classmethod
    def _filter_bad_variant_data(cls, sequences, all_metadata, out_prefix, removed_sequences):
        genes_to_remove = set()
        variants_only_genes_found_variant = set()
        log_file = out_prefix + '.check_metadata.log'
        tsv_file = out_prefix + '.check_metadata.tsv'
        log_fh = pyfastaq.utils.open_file_write(log_file)

        for sequence_name, metadata_dict in sorted(all_metadata.items()):
            if sequence_name in removed_sequences:
                print(sequence_name, 'was removed because does not look like a gene, so removing its metadata', file=log_fh)
                del all_metadata[sequence_name]
                continue

            assert sequence_name in sequences

            # if this is non_coding, we shouldn't have any amino acid variants
            if metadata_dict['seq_type'] != 'p':
                for position in metadata_dict['p']:
                    for metadata in metadata_dict['p'][position]:
                        print(sequence_name, 'variant is an amino acid change, but sequence is non-coding. Removing. Line of file was:', metadata, file=log_fh)

                metadata_dict['p'] = {}

            # take out variant metadata that doesn't make sense (eg bases not matching ref sequence)
            for variant_type in ['n', 'p']:
                positions_to_remove = []
                for position in metadata_dict[variant_type]:
                    meta_to_remove = []
                    for metadata in metadata_dict[variant_type][position]:
                        to_translate = variant_type == 'p'

                        if not metadata.variant.sanity_check_against_seq(sequences[sequence_name], translate_seq=to_translate):
                            print(sequence_name, 'variant does not match reference. Removing. Line of file was:', metadata, file=log_fh)
                            meta_to_remove.append(metadata)
                            continue

                        if metadata_dict['variant_only']:
                            variants_only_genes_found_variant.add(sequence_name)

                    for metadata in meta_to_remove:
                        metadata_dict[variant_type][position].remove(metadata)
                    if len(metadata_dict[variant_type][position]) == 0:
                        positions_to_remove.append(position)

                for position in positions_to_remove:
                    del metadata_dict[variant_type][position]

            if metadata_dict['variant_only'] and len(metadata_dict['n']) == len(metadata_dict['p']) == len(metadata_dict['.']) == 0:
                print(sequence_name, 'No remaining data after checks. Removing this sequence because it is variants only', file=log_fh)
                genes_to_remove.add(sequence_name)

        for sequence_name in genes_to_remove:
            del all_metadata[sequence_name]
            del sequences[sequence_name]

        pyfastaq.utils.close(log_fh)
        ReferenceData._write_metadata_tsv(all_metadata, tsv_file)


    @classmethod
    def _try_to_get_gene_seq(cls, seq, min_length, max_length):
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


    @classmethod
    def _remove_bad_genes(cls, sequences, metadata, log_file, min_gene_length, max_gene_length):
        to_remove = set()

        if len(sequences) == 0:
            return to_remove

        log_fh = pyfastaq.utils.open_file_write(log_file)

        for name in sorted(sequences):
            if metadata[name]['seq_type'] != 'p':
                continue

            new_seq, message = ReferenceData._try_to_get_gene_seq(sequences[name], min_gene_length, max_gene_length)
            if new_seq is None:
                to_remove.add(name)
            else:
                sequences[name] = new_seq

            if message is not None:
                print(name, message, file=log_fh)

        pyfastaq.utils.close(log_fh)

        for name in to_remove:
            sequences.pop(name)

        return to_remove


    def sanity_check(self, outprefix):
        removed_seqs = self._remove_bad_genes(self.sequences, self.metadata, outprefix + '.check_genes.log', self.min_gene_length, self.max_gene_length)
        ReferenceData._filter_bad_variant_data(self.sequences, self.metadata, outprefix, removed_seqs)


    @classmethod
    def _new_seq_name(cls, name):
        assert len(name.split()) == 1 and name.strip() == name
        return re.sub(rename_sub_regex, '_', name)


    @classmethod
    def _seq_names_to_rename_dict(cls, names):
        used_names = set()
        old_name_to_new = {}

        for old_name in sorted(names):
            assert len(old_name.split()) == 1 and old_name.strip() == old_name
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
    def _rename_names_in_seq_dict(cls, seq_dict, rename_dict):
        new_dict = {}
        for name, seq in seq_dict.items():
            if name in rename_dict:
                seq.id = rename_dict[name]

            new_dict[seq.id] = seq
        return new_dict

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
        self.rename_dict = ReferenceData._seq_names_to_rename_dict(self.sequences.keys())
        if len(self.rename_dict):
            with open(outfile, 'w') as f:
                for old_name, new_name in sorted(self.rename_dict.items()):
                    print(old_name, new_name, sep='\t', file=f)

            self.sequences = ReferenceData._rename_names_in_seq_dict(self.sequences, self.rename_dict)
            self.metadata = ReferenceData._rename_names_in_metadata(self.metadata, self.rename_dict)


    def sequence_type(self, sequence_name):
        assert sequence_name in self.metadata
        return self.metadata[sequence_name]['seq_type'], self.metadata[sequence_name]['variant_only']


    def sequence(self, sequence_name):
        return self.sequences.get(sequence_name, None)


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

        for cluster, name_set in sorted(clusters.items()):
            seq_names = sorted(list(name_set))
            print(cluster, *seq_names, sep='\t', file=f_out)

        pyfastaq.utils.close(f_out)


    def cluster_with_cdhit(self, outprefix, seq_identity_threshold=0.9, threads=1, length_diff_cutoff=0.9, nocluster=False, verbose=False, clusters_file=None):
        clusters = {}
        ReferenceData._write_sequences_to_files(self.sequences, self.metadata, outprefix)
        ref_types = ('noncoding', 'noncoding.varonly', 'gene', 'gene.varonly')

        for ref_type in ref_types:
            ref_file = outprefix + '.' + ref_type + '.fa'
            if os.path.getsize(ref_file) == 0:
                continue

            if len(clusters) == 0:
                min_cluster_number = 0
            else:
                min_cluster_number = 1 + max([int(x) for x in clusters.keys()])

            cdhit_runner = cdhit.Runner(
              ref_file,
              seq_identity_threshold=seq_identity_threshold,
              threads=threads,
              length_diff_cutoff=length_diff_cutoff,
              verbose=verbose,
              min_cluster_number = min_cluster_number,
            )

            if clusters_file is not None:
                new_clusters = cdhit_runner.run_get_clusters_from_file(clusters_file, self.sequences, rename_dict=self.rename_dict)
            elif nocluster:
                new_clusters = cdhit_runner.fake_run()
            else:
                new_clusters = cdhit_runner.run()

            clusters.update(new_clusters)

        self.write_cluster_allocation_file(clusters, outprefix + '.clusters.tsv')
        return clusters


    def write_seqs_to_fasta(self, outfile, names):
        f_out = pyfastaq.utils.open_file_write(outfile)

        for name in sorted(names):
            print(self.sequence(name), file=f_out)

        pyfastaq.utils.close(f_out)
