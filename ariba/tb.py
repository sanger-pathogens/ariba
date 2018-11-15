import csv
import json

from Bio import SeqIO

from ariba import flag


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
