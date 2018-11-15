import csv
import json
import os
import re
import sys

from Bio import SeqIO

from ariba import flag

data_dir = os.path.join(os.path.dirname(__file__), 'tb_data')
assert os.path.exists(data_dir)


def report_to_resistance_dict(infile):
    '''Takes final ariba report.tsv file, and extracts
    resistance calls, returning a dict of
    drug name -> list of mutations.
    each "mutation" in the list is a tuple of (gene name, mutation).
    Mutation is of the form X42Y, or "incomplete_gene" for katG and
    pncA when they are not complete.
    This all assumes that the reference data are in the "correct"
    form, where the variant descriptions in the var_description column of the
    TSV file ends with a comma-separated list of the drug names'''
    complete_genes = {'katG': 'Isoniazid', 'pncA': 'Pyrazinamide'}
    res_calls = {}
    incomplete_genes = set()
    with open(infile) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for d in reader:
            if d['ref_name'] in complete_genes and d['gene'] == '1':
                f = flag.Flag(int(d['flag']))
                if not f.has('complete_gene'):
                    incomplete_genes.add(d['ref_name'])

            if d['has_known_var'] == '1':
                drugs = d['var_description'].split()[-1].split(',')
                for drug in drugs:
                    if drug not in res_calls:
                        res_calls[drug] = []
                    res_calls[drug].append((d['ref_name'], d['known_var_change']))

    for gene in incomplete_genes:
        drug = complete_genes[gene]
        if drug not in res_calls:
            res_calls[drug] = []
        res_calls[drug].append((gene, 'Incomplete_gene'))

    return res_calls


def genbank_to_gene_coords(infile, genes):
    '''Input file in genbank format. genes = list of gene names to find.
    Returns dict of gene name -> {start: x, end: y}, where x and y are
    zero-based. x<y iff gene is on forwards strand'''
    coords = {}

    for seq_record in SeqIO.parse(infile, "genbank"):
        for feature in seq_record.features:
            if feature.type == 'gene':
                gene_name = feature.qualifiers.get('gene', [None])[0]
                if gene_name not in genes:
                    continue

                if feature.location.strand == 1:
                    coords[gene_name] = {'start': int(feature.location.start), 'end': int(feature.location.end) - 1}
                else:
                    coords[gene_name] = {'end': int(feature.location.start), 'start': int(feature.location.end) - 1}

    return coords


def load_mutations(gene_coords, mutation_to_drug_json, variants_txt, upstream_before=100):
    '''Load mutations from "mykrobe-style" files. mutation_to_drug_json is json file
    of mutation -> list of drugs. variants_txt is text file of variants used my mykrobe's
    make probes. gene_coords should be dict of gene coords made by the function
    genbank_to_gene_coords'''
    with open(mutation_to_drug_json) as f:
        drug_data = json.load(f)

    mutations = []
    genes_with_indels = set()
    genes_need_upstream = set()
    genes_non_upstream = set()

    with open(variants_txt) as f:
        for line in f:
            gene, variant, d_or_p = line.rstrip().split('\t')
            coding = 0 if gene == 'rrs' else 1
            d = {'gene': gene, 'var': variant, 'coding': coding, 'upstream': False}
            drug_data_key = d['gene'] + '_' + d['var']
            if drug_data_key not in drug_data:
                print('KEY', drug_data_key, 'NOT FOUND', file=sys.stderr)
            else:
                d['drugs'] = ','.join(sorted(drug_data[drug_data_key]))

            if d_or_p == 'DNA' and gene != 'rrs':
                assert gene != 'rrs'
                re_match = re.match('([ACGT]+)(-?[0-9]+)([ACGTX]+)', d['var'])
                try:
                    ref, pos, alt = re_match.groups()
                except:
                    print('regex error:', d['var'], file=sys.stderr)
                    continue

                pos = int(pos)
                if len(ref) != len(alt):
                    genes_with_indels.add(d['gene'])
                    continue
                elif pos > 0:
                    print('ignoring:', d, file=sys.stderr)
                    continue
                elif pos < 0:
                    this_gene_coords = gene_coords[d['gene']]
                    d['upstream'] = True
                    if this_gene_coords['start'] < this_gene_coords['end']:
                        variant_pos_in_output_seq = upstream_before + pos + 1
                    else:
                        variant_pos_in_output_seq = upstream_before + pos + 1
                    assert variant_pos_in_output_seq > 0
                    d['var'] = ref + str(variant_pos_in_output_seq) + alt
                    d['original_mutation'] = variant
                    genes_need_upstream.add(d['gene'])
                elif pos == 0:
                    print('Zero coord!', d, file=sys.stderr)
                    continue
                else:
                    print('deal with?', d, file=sys.stderr)
                    continue

            mutations.append(d)
            if not d['upstream']:
                genes_non_upstream.add(d['gene'])

    return mutations, genes_with_indels, genes_need_upstream, genes_non_upstream

