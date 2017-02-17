import csv
import re
import os

class Error (Exception): pass

regex_string_to_float = re.compile(r'\s*(?P<lt_or_gt>[<>]?)\s*(?P<equals>=?)\s*(?P<number>[0-9.]+)\s*$')

class MicPlotter:
    def __init__(self, mic_file, summary_file):
        self.mic_file = mic_file
        self.summary_file = summary_file


    @classmethod
    def _mic_string_to_float(cls, s):
        regex_match = regex_string_to_float.match(s)

        if regex_match is None or regex_match.group('number') == '.':
            if s.strip() in {'NA', 'na', '', '.'}:
                return 'NA'
            else:
                return None

        try:
            flt = float(regex_match.group('number'))
        except:
            return None

        if regex_match.group('equals') == '':
            if regex_match.group('lt_or_gt') == '<':
                return 0.5 * flt
            elif regex_match.group('lt_or_gt') == '>':
                return 2 * flt

        return flt


    @classmethod
    def _load_mic_file(cls, infile):
        mic_data = {}

        with open(infile) as f:
            reader = csv.DictReader(f, delimiter='\t')
            if reader.fieldnames[0] != 'Sample':
                raise Error('Error. Expected first column of MIC file "' + infile + '" to be "Sample"')
                
            for row in reader:
                mic_data[row['Sample']] = {x: MicPlotter._mic_string_to_float(row[x]) for x in reader.fieldnames[1:]}

        return mic_data


    @classmethod
    def _load_summary_file(cls, infile):
        data = {}

        with open(infile) as f:
            reader = csv.DictReader(f, delimiter='\t')
            if reader.fieldnames[0] != 'name':
                raise Error('Error. Expected first column of summary file "' + infile + '" to be "name"')

            clusters = [x.split('.', maxsplit=1)[0] for x in reader.fieldnames[1:]]

            for row in reader:
                data[row['name']] = {}

                for field in row:
                    if field == 'name':
                        continue 

                    cluster, col = field.split('.', maxsplit=1)
                    if cluster not in clusters:
                        raise Error('Cluster "' + cluster + '" not recognised. Cannot continue')
                    if cluster not in data[row['name']]:
                        data[row['name']][cluster] = {}

                    try:
                        value = float(row[field])
                    except:
                        value = row[field]
                    data[row['name']][cluster][col] = value

        return data


    @classmethod
    def _to_boxplot_tsv(cls, summary_data, mic_data, antibiotic, outfile):
        ignore_columns = {'assembled', 'match', 'ref_seq', 'pct_id', 'known_var', 'novel_var'}
        with open(outfile, 'w') as f:
            print('Sample\tMIC\tMutations', file=f)

            for sample in sorted(summary_data):
                if sample not in mic_data:
                    raise Error('No MIC data found for sample "' + sample + '". Cannot continue')

                if antibiotic not in mic_data[sample]:
                    raise Error('Antibiotic "' + antibiotic + '" not found. Cannot continue')

                if mic_data[sample][antibiotic] == 'NA':
                    continue

                mutations = set()

                for cluster in summary_data[sample]:
                    if summary_data[sample][cluster]['assembled'] == 'interrupted':
                        mutations.add(cluster + '.interrupted')

                    for column, value in summary_data[sample][cluster].items():
                        if column in ignore_columns or column.endswith('.%'):
                            continue

                        if value == 'yes':
                            mutations.add(cluster + '.' + column)

                if len(mutations):
                    mutations = list(mutations)
                    mutations.sort()
                    mutations = '+'.join(mutations)
                else:
                    mutations = 'without_mutation'

                print(sample, mic_data[sample][antibiotic], mutations, sep='\t', file=f)

