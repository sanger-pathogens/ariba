import csv
import json
import os
import re
import sys
import tempfile

from Bio import SeqIO
import pyfastaq

from ariba import common, flag, ref_preparer

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
                if 'Original mutation' in d['var_description']:
                    drugs = d['var_description'].split(':')[-1].split('.')[0].split()[-1].split(',')
                    change = d['var_description'].split()[-1]
                else:
                    drugs = d['var_description'].split()[-1].split(',')
                    change = d['known_var_change']
                for drug in drugs:
                    if drug not in res_calls:
                        res_calls[drug] = []
                    res_calls[drug].append((d['ref_name'], change))

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
                    #print('ignoring synonymous change (not implemented):', d['gene'], d['var'], d['drugs'], file=sys.stderr)
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


def write_prepareref_fasta_file(outfile, gene_coords, genes_need_upstream, genes_non_upstream, upstream_before=100, upstream_after=100):
    '''Writes fasta file to be used with -f option of prepareref'''
    tmp_dict = {}
    fasta_in = os.path.join(data_dir, 'NC_000962.3.fa.gz')
    pyfastaq.tasks.file_to_dict(fasta_in, tmp_dict)
    ref_seq = tmp_dict['NC_000962.3']

    with open(outfile, 'w') as f:
        for gene in genes_non_upstream:
            start = gene_coords[gene]['start']
            end = gene_coords[gene]['end']
            if start < end:
                gene_fa = pyfastaq.sequences.Fasta(gene, ref_seq[start:end+1])
            else:
                gene_fa = pyfastaq.sequences.Fasta(gene, ref_seq[end:start+1])
                gene_fa.revcomp()

            print(gene_fa, file=f)

        for gene in genes_need_upstream:
            start = gene_coords[gene]['start']
            end = gene_coords[gene]['end']
            if start < end:
                gene_fa = pyfastaq.sequences.Fasta(gene, ref_seq[start - upstream_before:start + upstream_after])
            else:
                gene_fa = pyfastaq.sequences.Fasta(gene, ref_seq[start - upstream_after + 1:start + upstream_before + 1])
                gene_fa.revcomp()

            gene_fa.id += '_upstream'
            print(gene_fa, file=f)


def write_prepareref_metadata_file(mutations, outfile):
    aa_letters = {'G', 'P', 'A', 'V', 'L', 'I', 'M', 'C', 'F', 'Y',
        'W', 'H', 'K', 'R', 'Q', 'N', 'E', 'D', 'S', 'T'}

    with open(outfile, 'w') as f:
        for d in mutations:
            if d['upstream']:
                gene = d['gene'] + '_upstream'
                coding = '0'
            else:
                gene = d['gene']
                coding = d['coding']

            if 'original_mutation' in d:
                original_mutation_string = '. Original mutation ' + d['original_mutation']
            else:
                original_mutation_string = ''

            if d['var'].endswith('X'):
                ref = d['var'][0]
                if d['coding'] == 1 and not d['upstream']:
                    letters = aa_letters
                else:
                    letters = {'A', 'C', 'G', 'T'}

                assert ref in letters

                for x in letters:
                    if x == ref:
                        continue
                    variant = d['var'][:-1] + x
                    print(gene, coding, 1, variant, '.', 'Resistance to ' + d['drugs'] + original_mutation_string, sep='\t', file=f)
            else:
                print(gene, coding, 1, d['var'], '.', 'Resistance to ' + d['drugs'] + original_mutation_string, sep='\t', file=f)


def make_prepareref_files(outprefix):
    genbank_file = os.path.join(data_dir, 'NC_000962.3.gb')
    mut_to_drug_json = os.path.join(data_dir, 'panel.20181115.json')
    panel_txt_file = os.path.join(data_dir, 'panel.20181115.txt')
    fasta_out = outprefix + '.fa'
    tsv_out = outprefix + '.tsv'

    with open(panel_txt_file) as f:
        genes = set([x.split()[0] for x in f])

    ref_gene_coords = genbank_to_gene_coords(genbank_file, genes)
    mutations, genes_with_indels, genes_need_upstream, genes_non_upstream = load_mutations(ref_gene_coords, mut_to_drug_json, panel_txt_file)
    write_prepareref_fasta_file(fasta_out, ref_gene_coords, genes_need_upstream, genes_non_upstream)
    write_prepareref_metadata_file(mutations, tsv_out)


def make_prepareref_dir(outdir):
    if os.path.exists(outdir):
        raise Exception('Output directory ' + outdir + ' already exists. Cannot continue')

    tmpdir = tempfile.mkdtemp(prefix=outdir + '.tmp', dir=os.getcwd())
    tmp_prefix = os.path.join(tmpdir, 'out')
    make_prepareref_files(tmp_prefix)
    ref_prep = ref_preparer.RefPreparer(
        [tmp_prefix + '.fa'],
        None,
        metadata_tsv_files=[tmp_prefix + '.tsv'],
        run_cdhit=False,
        threads=1,
    )
    ref_prep.run(outdir)
    common.rmtree(tmpdir)

    json_data = {'tb': True}
    json_file = os.path.join(outdir, '00.params.json')
    with open(json_file, 'w') as f:
        json.dump(json_data, f)
