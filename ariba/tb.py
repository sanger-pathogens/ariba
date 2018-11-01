import csv
import json

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

