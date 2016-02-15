import pymummer

columns = [
    'ref_name',              # 0  name of reference sequence
    'ref_type',              # 1  type of reference sequence (presence/absence, variants only, noncoding)
    'flag',                  # 2  cluster flag
    'reads',                 # 3  number of reads in this cluster
    'cluster',               # 4  cluster number
    'ref_len',               # 5  length of reference sequence
    'ref_base_assembled',    # 6  number of reference nucleotides assembled by this scaffold
    'pc_ident',              # 7  %identity between ref sequence and scaffold
    'scaffold',              # 8  name of scaffold matching reference
    'scaff_len',             # 9  length of scaffold matching reference
    'known_var',             # 10 is this a known SNP from reference metadata? T|F
    'known_var_type',        # 11 if known_var=T, the type of variant. Currently only SNP supported
    'var_type',              # 12 n|p for nucleotide or protein
    'var_effect',            # 13 Effect of variant (SYN, NONSYN, FSHIFT... etc)
    'var_change',            # 14 amino acid or nucleotide change, eg I42L
    'var_start_scaff',       # 15 start position of variant in scaffold
    'var_end_scaff',         # 16 end position of variant in scaffold
    'var_scaff',             # 17 nucleotide in scaffold at variant position
    'var_scaff_read_depth',  # 18 total read depth at variant start position in scaffold, reported by mpileup
    'alt_scaff_bases',       # 19 alt bases on scaffold, reported by mpileup
    'alt_scaff_depth',       # 20 alt depth on scaffold, reported by mpileup
    'var_description',       # 21 description of variant from reference metdata
    'free_text',             # 22 other free text about reference sequence, from reference metadata
]


def header_line():
    return '\t'.join(columns)


def _samtools_depths_at_known_snps_all_wild(sequence_meta, contig_name, cluster, variant_list):
    '''Input is a known variants, as sequence_metadata object. The
       assumption is that both the reference and the assembly have the
       variant type, not wild type. The list variant_list should be a list
       of pymummer.variant.Variant objects, only contaning variants to the
       relevant query contig'''
    nuc_range = sequence_meta.variant.nucleotide_range()

    if nuc_range is None:
        return None

    depths = []

    for ref_position in range(nuc_range[0], nuc_range[1]+1, 1):
        nucmer_match = cluster.assembly_compare.nucmer_hit_containing_reference_position(cluster.assembly_compare.nucmer_hits, cluster.ref_sequence.id, ref_position)

        if nucmer_match is not None:
            # work out contig position. Needs indels variants to correct the position
            contig_position, in_indel = nucmer_match.qry_coords_from_ref_coord(ref_position, variant_list)
            ref, alt, total_depth, alt_depths = cluster.samtools_vars.get_depths_at_position(contig_name, contig_position)
            depths.append([str(x) for x in (contig_position + 1, contig_position + 1, ref, alt, total_depth, alt_depths)])

    return depths


def _report_lines_for_one_contig(cluster, contig_name, ref_cov_per_contig, pymummer_variants):
    lines = []

    common_first_columns = [
        cluster.ref_sequence.id,
        cluster.ref_sequence_type,
        str(cluster.status_flag),
        str(cluster.total_reads),
        cluster.name,
        str(len(cluster.ref_sequence)),
    ]

    if cluster.ref_sequence.id in cluster.refdata.metadata and  len(cluster.refdata.metadata[cluster.ref_sequence.id]['.']) > 0:
        free_text_columns = [x.free_text for x in cluster.refdata.metadata[cluster.ref_sequence.id]['.']]
    else:
        free_text_columns = ['.']

    if cluster.assembled_ok and contig_name in cluster.assembly_variants and len(cluster.assembly_variants[contig_name]) > 0:
        more_common_columns = [
            str(ref_cov_per_contig[contig_name]) if contig_name in ref_cov_per_contig else '0', # 6 ref bases assembled
            str(cluster.assembly_compare.percent_identities[contig_name]) if contig_name in cluster.assembly_compare.percent_identities else '0',
            contig_name,
            str(len(cluster.assembly.sequences[contig_name])),  # 9 length of scaffold matching reference
        ]

        for (position, var_type, var_change, var_effect, contributing_vars, matching_vars_set, metainfo_set) in cluster.assembly_variants[contig_name]:
            if len(matching_vars_set) > 0:
                known_var = 'T'
                known_var_type = 'SNP'
            else:
                known_var = 'F'
                known_var_type = '.'

            var_columns = ['.' if x is None else str(x) for x in [known_var, known_var_type, var_type, var_effect, var_change]]



            if len(matching_vars_set) > 0:
                matching_vars_column = ';;;'.join([x.to_string(separator='_') for x in matching_vars_set])
            else:
                matching_vars_column = '.'

            if contributing_vars is None:
                if len(matching_vars_set) > 0:
                    for matching_var in matching_vars_set:
                        depths_list = _samtools_depths_at_known_snps_all_wild(matching_var, contig_name, cluster, pymummer_variants)
                        if depths_list is None:
                            depths_list = [['.'] * 6]

                        for depths in depths_list:
                            lines.append('\t'.join(common_first_columns + more_common_columns + var_columns + depths + [matching_vars_column] + free_text_columns))
                else:
                    lines.append('\t'.join(
                        common_first_columns + more_common_columns + var_columns + \
                        ['.'] * (len(columns) - len(common_first_columns) - len(more_common_columns) - len(var_columns) - 2) + \
                        [matching_vars_column] + free_text_columns
                    ))
            else:
                for var in contributing_vars:
                    ref, alt, total_depth, alt_depths = cluster.samtools_vars.get_depths_at_position(contig_name, var.qry_start)
                    lines.append('\t'.join(
                        common_first_columns + more_common_columns + var_columns + \
                        [str(x) for x in (var.qry_start, var.qry_end, ref, alt, total_depth, alt_depths)] + \
                        [matching_vars_column] + free_text_columns
                    ))
    else:
        lines.append('\t'.join(common_first_columns + ['.'] * (len(columns) - len(common_first_columns) - 1) + free_text_columns))

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

    return lines if len(lines) > 0 else None

