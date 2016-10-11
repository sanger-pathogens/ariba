import copy
import re
import sys
import pymummer
from ariba import sequence_variant

class Error (Exception): pass

columns = [
    'ariba_ref_name',        # 0  ariba (renamed) name of reference sequence
    'ref_name',              # 1  original name of ref sequence
    'gene',                  # 2  is a gene 0|1
    'var_only',              # 3  is variant only 0|1
    'flag',                  # 4  cluster flag
    'reads',                 # 5  number of reads in this cluster
    'cluster',               # 6  name of cluster
    'ref_len',               # 7  length of reference sequence
    'ref_base_assembled',    # 8  number of reference nucleotides assembled by this contig
    'pc_ident',              # 9  %identity between ref sequence and contig
    'ctg',                   # 10  name of contig matching reference
    'ctg_len',               # 11  length of contig matching reference
    'ctg_cov',               # 12 mean mapped read depth of this contig
    'known_var',             # 13 is this a known SNP from reference metadata? 1|0
    'var_type',              # 14 The type of variant. Currently only SNP supported
    'var_seq_type',          # 15 if known_var=1, n|p for nucleotide or protein
    'known_var_change',      # 16 if known_var=1, the wild/variant change, eg I42L
    'has_known_var',         # 17 if known_var=1, 1|0 for whether or not the assembly has the variant
    'ref_ctg_change',        # 18 amino acid or nucleotide change between reference and contig, eg I42L
    'ref_ctg_effect',        # 19 effect of change between reference and contig, eg SYS, NONSYN (amino acid changes only)
    'ref_start',             # 20 start position of variant in contig
    'ref_end',               # 21 end position of variant in contig
    'ref_nt',                # 22 nucleotide(s) in contig at variant position
    'ctg_start',             # 23 start position of variant in contig
    'ctg_end',               # 24 end position of variant in contig
    'ctg_nt',                # 25 nucleotide(s) in contig at variant position
    'smtls_total_depth',     # 26 total read depth at variant start position in contig, reported by mpileup
    'smtls_nts',             # 27 alt nucleotides on contig, reported by mpileup
    'smtls_nts_depth',       # 28 alt depth on contig, reported by mpileup
    'var_description',       # 29 description of variant from reference metdata
    'free_text',             # 30 other free text about reference sequence, from reference metadata
]


var_columns = [
    'known_var',
    'var_type',
    'var_seq_type',
    'known_var_change',
    'has_known_var',
    'ref_ctg_change',
    'ref_ctg_effect',
    'ref_start',
    'ref_end',
    'ref_nt',
    'ctg_start',
    'ctg_end',
    'ctg_nt',
    'smtls_total_depth',
    'smtls_nts',
    'smtls_nts_depth',
    'var_description',
]

int_columns = [
    'reads',
    'ref_len',
    'ref_base_assembled',
    'ctg_len',
    'ref_start',
    'ref_end',
    'ctg_start',
    'ctg_end',
]


float_columns = [
    'ctg_cov',
    'pc_ident',
]


def header_line():
    return '\t'.join(columns)


def _samtools_depths_at_known_snps_all_wild(sequence_meta, contig_name, cluster, variant_list):
    '''Input is a known variants, as sequence_metadata object. The
       assumption is that both the reference and the assembly have the
       variant type, not wild type. The list variant_list should be a list
       of pymummer.variant.Variant objects, only contaning variants to the
       relevant query contig'''
    ref_nuc_range = sequence_meta.variant.nucleotide_range()

    if ref_nuc_range is None:
        return None

    bases = []
    ctg_nts = []
    ref_nts = []
    smtls_total_depths = []
    smtls_nts = []
    smtls_depths = []
    contig_positions = []

    for ref_position in range(ref_nuc_range[0], ref_nuc_range[1]+1, 1):
        nucmer_match = cluster.assembly_compare.nucmer_hit_containing_reference_position(cluster.assembly_compare.nucmer_hits, cluster.ref_sequence.id, ref_position)

        if nucmer_match is not None:
            # work out contig position. Needs indels variants to correct the position
            ref_nts.append(cluster.ref_sequence[ref_position])
            contig_position, in_indel = nucmer_match.qry_coords_from_ref_coord(ref_position, variant_list)
            contig_positions.append(contig_position)
            bases, total_depth, base_depths = cluster.samtools_vars.get_depths_at_position(contig_name, contig_position)
            #ctg_nts.append(ref)
            #samtools_nts.append(bases)
            ctg_nts.append(cluster.assembly.sequences[contig_name][contig_position])
            smtls_nts.append(bases)
            smtls_total_depths.append(total_depth)
            smtls_depths.append(base_depths)

    ctg_nts = ';'.join(ctg_nts) if len(ctg_nts) else '.'
    ref_nts = ';'.join(ref_nts) if len(ref_nts) else '.'
    smtls_nts = ';'.join(smtls_nts) if len(smtls_nts) else '.'
    smtls_total_depths = ';'.join([str(x)for x in smtls_total_depths]) if len(smtls_total_depths) else '.'
    smtls_depths = ';'.join([str(x)for x in smtls_depths]) if len(smtls_depths) else '.'
    ctg_start = str(min(contig_positions) + 1) if contig_positions is not None else '.'
    ctg_end = str(max(contig_positions) + 1) if contig_positions is not None else '.'

    return [str(x) for x in [
        ref_nuc_range[0] + 1,
        ref_nuc_range[1] + 1,
        ref_nts,
        ctg_start,
        ctg_end,
        ctg_nts,
        smtls_total_depths,
        smtls_nts,
        smtls_depths
    ]]


def _report_lines_for_one_contig(cluster, contig_name, ref_cov_per_contig, pymummer_variants):
    lines = []
    reported_known_vars = set()
    contig_length = len(cluster.assembly.sequences[contig_name])
    assert contig_length != 0

    if contig_name in ref_cov_per_contig:
        if contig_name == cluster.assembly_compare.scaff_name_matching_ref and cluster.assembly_compare.gene_matching_ref_type == 'GENE_FOUND':
            ref_cov = len(cluster.ref_sequence)
        else:
            ref_cov = ref_cov_per_contig[contig_name]
    else:
        ref_cov = 0

    common_first_columns = [
        cluster.ref_sequence.id,
        cluster.refdata.ariba_to_original_name.get(cluster.ref_sequence.id, cluster.ref_sequence.id),
        cluster.is_gene,
        cluster.is_variant_only,
        str(cluster.status_flag),
        str(cluster.total_reads),
        cluster.name,
        str(len(cluster.ref_sequence)),
        str(ref_cov),
        str(cluster.assembly_compare.percent_identities[contig_name]) if contig_name in cluster.assembly_compare.percent_identities else '0',
        contig_name,
        str(contig_length),  # 9 length of scaffold matching reference
    ]

    # it's possible that there is no read depth on an assembled contig
    if contig_name in cluster.total_contig_depths:
        common_first_columns.append(str(round(cluster.total_contig_depths[contig_name] / contig_length, 1)))
    else:
        common_first_columns.append('0')

    if cluster.ref_sequence.id in cluster.refdata.metadata and  len(cluster.refdata.metadata[cluster.ref_sequence.id]['.']) > 0:
        free_text_column = ';'.join([x.free_text for x in cluster.refdata.metadata[cluster.ref_sequence.id]['.']])
    else:
        free_text_column = ';'.join(['.'])

    remaining_samtools_variants = copy.copy(cluster.variants_from_samtools)
    if cluster.assembled_ok and contig_name in cluster.assembly_variants and len(cluster.assembly_variants[contig_name]) > 0:
        for (position, var_seq_type, ref_ctg_change, var_effect, contributing_vars, matching_vars_set, metainfo_set) in cluster.assembly_variants[contig_name]:
            if len(matching_vars_set) > 0:
                is_known_var = '1'
                known_var_change = 'unknown'
                var_type = 'SNP'
                has_known_var = '1'
                matching_vars_column = ';;;'.join([x.to_string(separator=':') for x in matching_vars_set])
            else:
                is_known_var = '0'
                known_var_change = '.'
                has_known_var = '0'
                var_type = '.'
                matching_vars_column = '.'

            variant_columns = ['.' if x is None else str(x) for x in [is_known_var, var_type, var_seq_type, known_var_change, has_known_var, ref_ctg_change, var_effect]]

            if contributing_vars is None:
                samtools_columns = [['.'] * 9]
            else:
                if var_effect in ['INDELS', 'MULTIPLE']:
                    ref_start_pos = min([x.ref_start for x in contributing_vars])
                    ref_end_pos = max([x.ref_start for x in contributing_vars])
                    ctg_start_pos = min([x.qry_start for x in contributing_vars])
                    ctg_end_pos = max([x.qry_start for x in contributing_vars])
                else:
                    ref_start_pos = 3 * position if cluster.is_gene == '1' else position
                    assert contig_name in cluster.assembly_compare.nucmer_hits
                    ref_start_hit = None
                    for hit in cluster.assembly_compare.nucmer_hits[contig_name]:
                        if hit.ref_name == cluster.ref_sequence.id and hit.ref_coords().distance_to_point(ref_start_pos) == 0:
                            ref_start_hit = copy.copy(hit)
                            break

                    assert ref_start_hit is not None
                    ctg_start_pos, ctg_start_in_indel = ref_start_hit.qry_coords_from_ref_coord(ref_start_pos, pymummer_variants)

                    if known_var_change not in  ['.', 'unknown']:
                        regex = re.match('^([^0-9]+)([0-9]+)([^0-9]+)$', known_var_change)
                        try:
                            ref_var_string, ref_var_position, ctg_var_string = regex.group(1, 2, 3)
                        except:
                            raise Error('Error parsing variant ' + known_var_change)
                    elif ref_ctg_change != '.':
                        if '_' in ref_ctg_change:
                            regex = re.match('^([^0-9]+)([0-9]+)_[^0-9]+[0-9]+([^0-9]+)$', ref_ctg_change)
                            try:
                                ref_var_string, ref_var_position, ctg_var_string = regex.group(1, 2, 3)
                            except:
                                raise Error('Error parsing variant ' + ref_ctg_change)
                        else:
                            regex = re.match('^([^0-9]+)([0-9]+)([^0-9]+)$', ref_ctg_change)
                            try:
                                ref_var_string, ref_var_position, ctg_var_string = regex.group(1, 2, 3)
                            except:
                                raise Error('Error parsing variant ' + ref_ctg_change)
                    else:
                        assert var_effect == 'SYN'

                    if var_effect == 'SYN':
                        ref_end_pos = ref_start_pos + 2
                        ctg_end_pos = ctg_start_pos + 2
                    elif ref_var_string == '.' or var_effect in {'INS', 'DEL', 'FSHIFT', 'TRUNC', 'INDELS', 'UNKNOWN'}:
                        ref_end_pos = ref_start_pos
                        ctg_end_pos = ctg_start_pos
                    elif cluster.is_gene == '1':
                        ref_end_pos = ref_start_pos + 3 * len(ref_var_string) - 1
                        ctg_end_pos = ctg_start_pos + 3 * len(ctg_var_string) - 1
                    else:
                        ref_end_pos = ref_start_pos + len(ref_var_string) - 1
                        ctg_end_pos = ctg_start_pos + len(ctg_var_string) - 1

                smtls_total_depth = []
                smtls_alt_nt = []
                smtls_alt_depth = []

                for qry_pos in range(ctg_start_pos, ctg_end_pos + 1, 1):
                    if contig_name in remaining_samtools_variants:
                        try:
                            remaining_samtools_variants[contig_name].discard(qry_pos)
                        except:
                            pass

                    depths_tuple = cluster.samtools_vars.get_depths_at_position(contig_name, qry_pos)

                    if depths_tuple is not None:
                        smtls_alt_nt.append(depths_tuple[0])
                        smtls_total_depth.append(str(depths_tuple[1]))
                        smtls_alt_depth.append(str(depths_tuple[2]))

                smtls_total_depth = ';'.join(smtls_total_depth) if len(smtls_total_depth) else '.'
                smtls_alt_nt = ';'.join(smtls_alt_nt) if len(smtls_alt_nt) else '.'
                smtls_alt_depth = ';'.join(smtls_alt_depth) if len(smtls_alt_depth) else '.'
                samtools_columns = [
                        str(ref_start_pos + 1), #ref_start
                        str(ref_end_pos + 1), # ref_end
                        cluster.ref_sequence[ref_start_pos:ref_end_pos+1],
                        str(ctg_start_pos + 1),  # ctg_start
                        str(ctg_end_pos + 1),  #ctg_end
                        cluster.assembly.sequences[contig_name][ctg_start_pos:ctg_end_pos + 1], # ctg_nt
                        smtls_total_depth,
                        smtls_alt_nt,
                        smtls_alt_depth,
                ]


            if len(matching_vars_set) > 0:
                for matching_var in matching_vars_set:
                    if contributing_vars is None:
                        samtools_columns = _samtools_depths_at_known_snps_all_wild(matching_var, contig_name, cluster, pymummer_variants)
                        samtools_columns[2] = samtools_columns[2].replace(';', '')
                        samtools_columns[5] = samtools_columns[5].replace(';', '')
                    reported_known_vars.add(str(matching_var.variant))
                    variant_columns[3] = str(matching_var.variant)

                    if matching_var.has_variant(cluster.ref_sequence) == (ref_ctg_change is not None):
                        variant_columns[4] = '0'
                    else:
                        variant_columns[4] = '1'

                    if samtools_columns is None:
                        samtools_columns = [['.'] * 9]
                    else:
                        for ctg_pos in range(int(samtools_columns[3]) - 1, int(samtools_columns[4]), 1):
                            if contig_name in remaining_samtools_variants:
                                remaining_samtools_variants[contig_name].discard(ctg_pos)

                    lines.append('\t'.join(common_first_columns + variant_columns + samtools_columns + [matching_vars_column] + [free_text_column]))
            else:
                lines.append('\t'.join(
                    common_first_columns + variant_columns + \
                    samtools_columns + \
                    [matching_vars_column] + [free_text_column]
                ))

    if contig_name in remaining_samtools_variants:
        for var_position in remaining_samtools_variants[contig_name]:
            depths_tuple = cluster.samtools_vars.get_depths_at_position(contig_name, var_position)

            if depths_tuple is not None:
                ref_coord, in_indel = None, None
                if contig_name in cluster.assembly_compare.nucmer_hits:
                    for hit in cluster.assembly_compare.nucmer_hits[contig_name]:
                        if hit.qry_coords().distance_to_point(var_position) == 0:
                            ref_coord, in_indel = hit.ref_coords_from_qry_coord(var_position, pymummer_variants)
                            break

                if ref_coord is None:
                    ref_coord = '.'
                    ref_nt = '.'
                    var_string = None
                else:
                    ref_nt = cluster.ref_sequence[ref_coord]
                    ctg_nt = cluster.assembly.sequences[contig_name][var_position]
                    alt_strings = [x for x in depths_tuple[0].split(',') if x != ctg_nt]
                    var_string = ctg_nt + str(ref_coord + 1) + ','.join(alt_strings)
                    ref_coord = str(ref_coord + 1)

                if var_string not in reported_known_vars:
                    new_cols = [
                        '0',  # known_var column
                        'HET', # var_type
                        '.', '.', '.', var_string, '.', ref_coord, ref_coord, ref_nt, # var_seq_type ... ref_nt
                        str(var_position + 1), str(var_position + 1), # ctg_start, ctg_end
                        ctg_nt, # ctg_nt
                        str(depths_tuple[1]), # smtls_total_depth
                        depths_tuple[0], # smtls_alt_nt
                        str(depths_tuple[2]), # smtls_alt_depth
                        '.',
                        free_text_column,
                    ]
                    lines.append('\t'.join(common_first_columns + new_cols))

    if len(lines) == 0:
        lines.append('\t'.join(common_first_columns + ['.'] * (len(columns) - len(common_first_columns) - 1) + [free_text_column]))

    return lines


def report_lines(cluster):
    if cluster.status_flag.has('ref_seq_choose_fail'):
        fields = ['.', '.', '.', '.', str(cluster.status_flag), str(cluster.total_reads), cluster.name] + ['.'] * (len(columns) - 7)
        assert len(fields) == len(columns)
        return ['\t'.join(fields)]
    elif cluster.status_flag.has('assembly_fail'):
        fields = ['.', '.', '.', '.', str(cluster.status_flag), str(cluster.total_reads), cluster.name] + ['.'] * (len(columns) - 7)
        assert len(fields) == len(columns)
        return ['\t'.join(fields)]

    ref_cov_per_contig = cluster.assembly_compare.ref_cov_per_contig(cluster.assembly_compare.nucmer_hits)
    lines = []
    pymummer_variants = pymummer.snp_file.get_all_variants(cluster.assembly_compare.nucmer_snps_file)

    for contig_name in sorted(cluster.assembly.sequences):
        contig_pymummer_variants = [x for x in pymummer_variants if x.qry_name == contig_name]
        lines.extend(_report_lines_for_one_contig(cluster, contig_name, ref_cov_per_contig, contig_pymummer_variants))

    lines_ok = True

    for line in lines:
        if len(line.split('\t')) != len(columns):
            cols = line.split('\t')
            print('Error making report - wrong number of columns. Expected', len(columns), 'but got', len(cols), file=sys.stderr)
            for i in range(len(cols)):
                print(i, cols[i], sep='\t', file=sys.stderr)
            lines_ok = False

    if not lines_ok:
        raise Error('Error making report. Cannot continue')

    return lines if len(lines) > 0 else None

