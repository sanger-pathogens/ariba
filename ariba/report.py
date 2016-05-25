import sys
import pymummer

class Error (Exception): pass

columns = [
    'ref_name',              # 0  name of reference sequence
    'ref_type',              # 1  type of reference sequence (presence/absence, variants only, noncoding)
    'flag',                  # 2  cluster flag
    'reads',                 # 3  number of reads in this cluster
    'cluster',               # 4  name of cluster
    'ref_len',               # 5  length of reference sequence
    'ref_base_assembled',    # 6  number of reference nucleotides assembled by this contig
    'pc_ident',              # 7  %identity between ref sequence and contig
    'ctg',                   # 8  name of contig matching reference
    'ctg_len',               # 9  length of contig matching reference
    'ctg_cov',               # 10 mean mapped read depth of this contig
    'known_var',             # 11 is this a known SNP from reference metadata? 1|0
    'var_type',              # 12 The type of variant. Currently only SNP supported
    'var_seq_type',          # 13 if known_var=1, n|p for nucleotide or protein
    'known_var_change',      # 14 if known_var=1, the wild/variant change, eg I42L
    'has_known_var',         # 15 if known_var=1, 1|0 for whether or not the assembly has the variant
    'ref_ctg_change',        # 16 amino acid or nucleotide change between reference and contig, eg I42L
    'ref_ctg_effect',        # 17 effect of change between reference and contig, eg SYS, NONSYN (amino acid changes only)
    'ref_start',             # 18 start position of variant in contig
    'ref_end',               # 19 end position of variant in contig
    'ref_nt',                # 20 nucleotide(s) in contig at variant position
    'ctg_start',             # 21 start position of variant in contig
    'ctg_end',               # 22 end position of variant in contig
    'ctg_nt',                # 23 nucleotide(s) in contig at variant position
    'smtls_total_depth',     # 24 total read depth at variant start position in contig, reported by mpileup
    'smtls_alt_nt',          # 25 alt nucleotides on contig, reported by mpileup
    'smtls_alt_depth',       # 26 alt depth on contig, reported by mpileup
    'var_description',       # 27 description of variant from reference metdata
    'free_text',             # 28 other free text about reference sequence, from reference metadata
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
    'smtls_alt_nt',
    'smtls_alt_depth',
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

    ctg_nts = []
    ref_nts = []
    smtls_total_depths = []
    smtls_alt_nts = []
    smtls_alt_depths = []
    contig_positions = []

    for ref_position in range(ref_nuc_range[0], ref_nuc_range[1]+1, 1):
        nucmer_match = cluster.assembly_compare.nucmer_hit_containing_reference_position(cluster.assembly_compare.nucmer_hits, cluster.ref_sequence.id, ref_position)

        if nucmer_match is not None:
            # work out contig position. Needs indels variants to correct the position
            ref_nts.append(cluster.ref_sequence[ref_position])
            contig_position, in_indel = nucmer_match.qry_coords_from_ref_coord(ref_position, variant_list)
            contig_positions.append(contig_position)
            ref, alt, total_depth, alt_depths = cluster.samtools_vars.get_depths_at_position(contig_name, contig_position)
            ctg_nts.append(ref)
            smtls_alt_nts.append(alt)
            smtls_total_depths.append(total_depth)
            smtls_alt_depths.append(alt_depths)

    ctg_nts = ';'.join(ctg_nts) if len(ctg_nts) else '.'
    ref_nts = ';'.join(ref_nts) if len(ref_nts) else '.'
    smtls_alt_nts = ';'.join(smtls_alt_nts) if len(smtls_alt_nts) else '.'
    smtls_total_depths = ';'.join([str(x)for x in smtls_total_depths]) if len(smtls_total_depths) else '.'
    smtls_alt_depths = ';'.join([str(x)for x in smtls_alt_depths]) if len(smtls_alt_depths) else '.'
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
        smtls_alt_nts,
        smtls_alt_depths
    ]]


def _report_lines_for_one_contig(cluster, contig_name, ref_cov_per_contig, pymummer_variants):
    lines = []
    contig_length = len(cluster.assembly.sequences[contig_name])
    assert contig_length != 0

    common_first_columns = [
        cluster.ref_sequence.id,
        cluster.ref_sequence_type,
        str(cluster.status_flag),
        str(cluster.total_reads),
        cluster.name,
        str(len(cluster.ref_sequence)),
        str(ref_cov_per_contig[contig_name]) if contig_name in ref_cov_per_contig else '0', # 6 ref bases assembled
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
                contributing_vars.sort(key = lambda x: x.qry_start)

                smtls_total_depth = []
                smtls_alt_nt = []
                smtls_alt_depth = []

                for var in contributing_vars:
                    depths_tuple = cluster.samtools_vars.get_depths_at_position(contig_name, var.qry_start)
                    if depths_tuple is not None:
                        smtls_alt_nt.append(depths_tuple[1])
                        smtls_total_depth.append(str(depths_tuple[2]))
                        smtls_alt_depth.append(str(depths_tuple[3]))

                smtls_total_depth = ';'.join(smtls_total_depth) if len(smtls_total_depth) else '.'
                smtls_alt_nt = ';'.join(smtls_alt_nt) if len(smtls_alt_nt) else '.'
                smtls_alt_depth = ';'.join(smtls_alt_depth) if len(smtls_alt_depth) else '.'
                samtools_columns = [
                        str(contributing_vars[0].ref_start), #ref_start
                        str(contributing_vars[0].ref_end), # ref_end
                        ';'.join([x.ref_base for x in contributing_vars]), # ref_nt
                        str(contributing_vars[0].qry_start),  # ctg_start
                        str(contributing_vars[0].qry_end),  #ctg_end
                        ';'.join([x.qry_base for x in contributing_vars]), #ctg_nt
                        smtls_total_depth,
                        smtls_alt_nt,
                        smtls_alt_depth,
                ]


            if len(matching_vars_set) > 0:
                for matching_var in matching_vars_set:
                    if contributing_vars is None:
                        samtools_columns = _samtools_depths_at_known_snps_all_wild(matching_var, contig_name, cluster, pymummer_variants)
                    variant_columns[3] = str(matching_var.variant)

                    if matching_var.has_variant(cluster.ref_sequence) == (ref_ctg_change is not None):
                        variant_columns[4] = '0'
                    else:
                        variant_columns[4] = '1'

                    if samtools_columns is None:
                        samtools_columns = [['.'] * 9]

                    lines.append('\t'.join(common_first_columns + variant_columns + samtools_columns + [matching_vars_column] + [free_text_column]))
            else:
                lines.append('\t'.join(
                    common_first_columns + variant_columns + \
                    samtools_columns + \
                    [matching_vars_column] + [free_text_column]
                ))
    else:
        lines.append('\t'.join(common_first_columns + ['.'] * (len(columns) - len(common_first_columns) - 1) + [free_text_column]))

    return lines


def report_lines(cluster):
    if cluster.status_flag.has('ref_seq_choose_fail'):
        return ['\t'.join(['.', '.', str(cluster.status_flag), str(cluster.total_reads), cluster.name] + ['.'] * (len(columns) - 5))]
    elif cluster.status_flag.has('assembly_fail'):
        return ['\t'.join([cluster.ref_sequence.id, cluster.ref_sequence_type, str(cluster.status_flag), str(cluster.total_reads), cluster.name] + ['.'] * (len(columns) - 5))]


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

