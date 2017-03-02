import csv
import re
import os
import itertools
import collections
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as cmx
import math
import pyfastaq
from ariba import common

class Error (Exception): pass

regex_string_to_float = re.compile(r'\s*(?P<lt_or_gt>[<>]?)\s*(?P<equals>=?)\s*(?P<number>[0-9.]+)\s*$')

regex_position_from_var = re.compile(r'^[^0-9]*(?P<coord>[0-9]+)[^0-9]*$')

class MicPlotter:
    def __init__(self,
      antibiotic,
      mic_file,
      summary_file,
      outprefix,
      use_hets='yes',
      main_title=None,
      plot_height=7,
      plot_width=7,
      log_y=2,
      plot_types="points,violin",
      jitter_width=0.1,
      no_combinations=False,
      hlines='0.25,2',
      point_size=4,
      point_range='2,15',
      point_break='10,50,100,200,300',
      dot_size=100,
      dot_outline=False,
      dot_y_text_size=18,
      panel_heights='5,1',
      colourmap='Accent',
      number_of_colours=0,
      colour_skip=None,
      interrupted=False,
      violin_width=0.75,
      xkcd=False,
      min_samples=1
    ):
        self.antibiotic = antibiotic
        self.mic_file = mic_file
        self.summary_file = summary_file
        self.outprefix = outprefix

        allowed_use_hets = {'yes', 'no', 'exclude'}
        if not use_hets in allowed_use_hets:
            raise Error('Error in use_hets option. Allowed options are: ' + str(allowed_use_hets) + '. Got: ' + use_hets)
        self.use_hets = use_hets

        self.main_title = self.antibiotic if main_title is None else main_title
        self.plot_height = plot_height
        self.plot_width = plot_width
        self.log_y = log_y
        self.plot_types = set(plot_types.split(','))

        allowed_plot_types = {'point', 'violin', 'boxplot'}
        if not self.plot_types.issubset(allowed_plot_types):
            raise Error('Error in plot_types option. Allowed types are: ' + str(allowed_plot_types) + '. Got: ' +  str(self.plot_types))

        self.jitter_width = jitter_width
        self.no_combinations = no_combinations

        try:
            if len(hlines) == 0:
                self.hlines = []
            else:
                self.hlines = [float(x) for x in hlines.split(',')]
        except:
            raise Error('Error in hlines option. Needs to be a list of numbers separated by commas, or empty. Got this:\n' + hlines)

        self.point_size = point_size

        try:
            self.point_range = [int(x) for x in point_range.split(',')]
        except:
            raise Error('Error in point_range option. Needs to be of the form integer1,integer2. Got this:\n' + point_range)

        try:
            self.point_break = [int(x) for x in point_break.split(',')]
        except:
            raise Error('Error in point_break option. Needs to be comma-sparated list of integers. Got this:\n' + point_break)

        self.dot_size = dot_size
        self.dot_outline = dot_outline
        self.dot_y_text_size = dot_y_text_size

        try:
            self.panel_heights = [int(x) for x in panel_heights.split(',')]
        except:
            raise Error('Error in panel_heights option. Needs to be of the form integer1,integer2. Got this:\n' + panel_heights)

        self.colourmap = colourmap
        self.number_of_colours = number_of_colours

        if colour_skip is None:
            self.colour_skip = None
        else:
            try:
                self.colour_skip = [float(x) for x in colour_skip.split(',')]
            except:
                raise Error('Error in colour_skip option. Needs to be of the form a,b where 0 <= a < b <= 1. Got this:\n' + colour_skip)

        self.interrupted = interrupted
        self.violin_width = violin_width
        if xkcd:
            plt.xkcd()
        self.min_samples = min_samples


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
            reader = csv.DictReader(f, delimiter=',')

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
    def _get_colours(cls, total_length, number_of_colours, colormap, skip=None):
        if number_of_colours == 1:
            return ["black"] * total_length
        elif number_of_colours == 0:
            cmap = cmx.get_cmap(colormap)
            if skip is None:
                vals = [1.0 * x / (total_length - 1) for x in range(total_length)]
            else:
                assert len(skip) == 2 and 0 <= skip[0] <= 1 and 0 <= skip[1] <= 1
                if skip[-1] == 1:
                    vals = [skip[0] * x / (total_length - 1) for x in range(total_length)]
                elif skip[0] == 0:
                    vals = [skip[1] + (1 - skip[1]) * x / (total_length - 1) for x in range(total_length)]
                else:
                    length = 1 - (skip[1] - skip[0])
                    vals = [(length) * x / (total_length - 1) for x in range(total_length)]
                    vals = [x if x < skip[0] else x + (1-length) for x in vals]
                    
            return [cmap(x) for x in vals]
        else:
            cmap = cmx.get_cmap(colormap)
            colours = []
            for i in itertools.cycle(range(number_of_colours)):
                colours.append(cmap(i))
                if len(colours) >= total_length:
                    break
            return colours


    @classmethod
    def _get_top_plot_data(cls, summary_data, mic_data, antibiotic, use_hets, no_combinations=False, interrupted=False, outfile=None):
        assert use_hets in {'yes', 'no', 'exclude'}
        if outfile is not None:
            f = pyfastaq.utils.open_file_write(outfile)
            print('Sample\tMIC\tMutations', file=f)

        ignore_columns = {'assembled', 'match', 'ref_seq', 'pct_id', 'known_var', 'novel_var'}
        all_mutations = set()
        all_mutations_seen_combinations = set()
        top_plot_data = {} # cluster combination -> list of y coords (MIC values)

        for sample in sorted(summary_data):
            if sample not in mic_data:
                raise Error('No MIC data found for sample "' + sample + '". Cannot continue')
            if antibiotic not in mic_data[sample]:
                raise Error('Antibiotic "' + antibiotic + '" not found. Cannot continue')

            if mic_data[sample][antibiotic] == 'NA':
                continue

            mutations = set()
            found_het_and_exclude = False

            for cluster in summary_data[sample]:
                if summary_data[sample][cluster]['assembled'] == 'interrupted' and interrupted:
                    mutations.add(cluster + '.interrupted')

                for column, value in summary_data[sample][cluster].items():
                    if column in ignore_columns or column.endswith('.%'):
                        continue

                    if value == 'yes' or (use_hets == 'yes' and value == 'het'):
                        mutations.add(cluster + '.' + column.strip())
                    elif use_hets == 'exclude' and value == 'het':
                        found_het_and_exclude = True
                        break

                if found_het_and_exclude:
                    break

            if found_het_and_exclude:
                continue


            if len(mutations) == 0:
                mutations.add('without_mutation')

            all_mutations.update(mutations)
            mutations = list(mutations)
            mutations.sort()
            if no_combinations:
                for mutation in mutations:
                    all_mutations_seen_combinations.add((mutation,))
                    if mutation not in top_plot_data:
                        top_plot_data[mutation] = []
                    top_plot_data[mutation].append(mic_data[sample][antibiotic])
                    if outfile is not None:
                        print(sample, mic_data[sample][antibiotic], mutation, sep='\t', file=f)
            else:
                all_mutations_seen_combinations.add(tuple(mutations))
                mutations = '.'.join(mutations)
                if mutations not in top_plot_data:
                    top_plot_data[mutations] = []
                top_plot_data[mutations].append(mic_data[sample][antibiotic])
                if outfile is not None:
                    print(sample, mic_data[sample][antibiotic], mutations, sep='\t', file=f)


        if outfile is not None:
            pyfastaq.utils.close(f)

        return top_plot_data, all_mutations, all_mutations_seen_combinations


    @classmethod
    def _filter_top_plot_data(cls, top_plot_data, all_mutations, seen_combinations, min_samples):
        if min_samples == 1:
            return top_plot_data, all_mutations, seen_combinations

        new_top_plot_data = {}
        new_all_mutations = set()
        new_seen_combinations = set()

        for mutation_tuple in seen_combinations:
            mutation_string = '.'.join(mutation_tuple)
            mics = top_plot_data[mutation_string]

            if len(mics) >= min_samples:
                new_top_plot_data[mutation_string] = mics
                new_seen_combinations.add(mutation_tuple)
                new_all_mutations.update(mutation_tuple)

        return new_top_plot_data, new_all_mutations, new_seen_combinations


    @classmethod
    def _top_plot_y_ticks(cls, mic_data, antibiotic, log_y):
        mic_values = set()
        for sample in mic_data:
            mic = mic_data[sample][antibiotic]
            if mic not in [None, 'NA']:
                mic_values.add(mic)

        max_mic = max(mic_values)
        min_mic = min(mic_values)
        new_mic_values = []
        i = 1
        while i < max_mic * 2:
            new_mic_values.append(i)
            i *= 2

        i = 0.5
        while i > min_mic / 2:
            new_mic_values.append(i)
            i *= 0.5

        new_mic_values.sort()
        new_mic_values = [round(x, 4) for x in new_mic_values]

        if log_y > 0:
            tick_positions = [math.log(x, log_y) for x in new_mic_values]
        else:
            tick_positions = new_mic_values

        return tick_positions, new_mic_values


    @classmethod
    def _top_plot_scatter_counts(cls, mutations, top_plot_data, colours, log_y):
        x_coords = []
        y_coords = []
        sizes = []
        colour_list = []

        for i, mutation in enumerate(mutations):
            counts = collections.Counter(top_plot_data[mutation])
            for mic in sorted(counts):
                x_coords.append(i + 1)
                if log_y > 0:
                    y_coords.append(math.log(mic, log_y))
                else:
                    y_coords.append(mic)
                sizes.append(counts[mic])
                colour_list.append(colours[i])
                
        return x_coords, y_coords, sizes, colour_list


    @classmethod
    def _top_plot_scatter_data(cls, mutations, top_plot_data, colours, log_y):
        x_coords = []
        y_coords = []
        colour_list = []

        for i, mutation in enumerate(mutations):
            for mic in top_plot_data[mutation]:
                x_coords.append(i + 1)
                if log_y > 0:
                    y_coords.append(math.log(mic, log_y))
                else:
                    y_coords.append(mic)
                colour_list.append(colours[i])
                
        return x_coords, y_coords, colour_list


    @classmethod
    def _top_plot_violin_data(cls, mutations, top_plot_data, log_y):
        violin_data = []
        violin_pos = []

        for i, mutation in enumerate(mutations):
            if log_y > 0:
                violin_data.append([math.log(x, log_y) for x in top_plot_data[mutation]])
            else:
                violin_data.append(top_plot_data[mutation])
            violin_pos.append(i + 1)

        return violin_data, violin_pos

            
    @classmethod
    def _ordered_bottom_plot_rows(cls, mutations):
        l = []
        infinity = float('inf')

        for x in mutations:
            try:
                cluster, variant = x.split('.', maxsplit=1)
            except:
                l.append((x, infinity, x))
                continue

            if '.' in variant:
                try:
                    var_group, var = variant.split('.', maxsplit=1)
                except:
                    var_group = None
                    var = variant

                variant = var

            regex_match = regex_position_from_var.match(variant)
            if regex_match is not None and regex_match.group('coord') != '':
                coord = int(regex_match.group('coord'))
            else:
                coord = infinity

            l.append((cluster, coord, x))

        l.sort()
        return [x[-1] for x in l]


    @classmethod
    def _ordered_columns(cls, mutations, top_plot_data):
        # FIXME
        return sorted(list(mutations))


    @classmethod
    def _bottom_scatter_data(cls, bottom_plot_rows, columns, colours):
        x_coords = []
        y_coords = []
        colour_list = []

        for i, row in enumerate(bottom_plot_rows):
            for j, col in enumerate(columns):
                if row in col:
                    x_coords.append(j + 1)
                    y_coords.append(len(bottom_plot_rows) - i)
                    colour_list.append(colours[j])

        return x_coords, y_coords, colour_list


    @classmethod
    def _right_plot_data(cls, scatter_count_sizes, number_of_circles, x_pos):
        y_max = max(scatter_count_sizes)
        if y_max > 100:
            y_max = int(math.ceil(y_max / 100.0)) * 100
            sizes = [5, 50] + [x for x in range(100, y_max, 100)]
        else:
            y_max = int(math.ceil(y_max / 10.0)) * 10
            sizes = [5] + [x for x in range(10, y_max, 10)]
        x_coords = [x_pos] * len(sizes)
        y_coords = [x + 1 for x in range(len(sizes))]
        y_coords.reverse()
        return x_coords, y_coords, sizes


    def _make_plot(self, mic_data, top_plot_data, all_mutations, mut_combinations):
        bottom_plot_rows = MicPlotter._ordered_bottom_plot_rows(all_mutations)
        columns = MicPlotter._ordered_columns(mut_combinations, top_plot_data)
        colours = MicPlotter._get_colours(len(columns), self.number_of_colours, self.colourmap, self.colour_skip)
        bottom_scatter_x, bottom_scatter_y, bottom_colours = MicPlotter._bottom_scatter_data(bottom_plot_rows, columns, colours)
        columns = ['.'.join(x) for x in columns]
        assert len(colours) == len(columns)
        max_x = len(colours) + 1
    
        scatter_count_x, scatter_count_y, scatter_count_sizes, scatter_count_colours = MicPlotter._top_plot_scatter_counts(columns, top_plot_data, colours, self.log_y)
        scatter_data_x, scatter_data_y, scatter_data_colours = MicPlotter._top_plot_scatter_data(columns, top_plot_data, colours, self.log_y)
        violin_data, violin_positions = MicPlotter._top_plot_violin_data(columns, top_plot_data, self.log_y)

        # -------------------- SET UP GRID & PLOTS -----------------
        fig=plt.figure(figsize=(self.plot_width, self.plot_height))
        if self.point_size == 0:
            gs = gridspec.GridSpec(2, 2, height_ratios=self.panel_heights, width_ratios=[5,1])
        else:
            gs = gridspec.GridSpec(2, 1, height_ratios=self.panel_heights)

        plots=[]
        plots.append(plt.subplot(gs[0]))
        plots.append(plt.subplot(gs[1]))
        if self.point_size == 0:
            plots.append(plt.subplot(gs[2]))
            bottom_plot_index = 2
        else:
            bottom_plot_index = 1

        # ------------------------- TOP PLOT -----------------------
        for h in self.hlines:
            if self.log_y > 0:
                h = math.log(h, self.log_y)
            plots[0].hlines(h, 0, max_x, linestyle='--', linewidth=1, color='black')


        violins = plots[0].violinplot(violin_data, violin_positions, widths=self.violin_width, showmeans=False, showextrema=False, showmedians=False)
        for x, pc in enumerate(violins['bodies']):
            pc.set_facecolor(colours[x])
            pc.set_edgecolor(colours[x])

        if self.point_size == 0:
            plots[0].scatter(scatter_count_x, scatter_count_y, s=scatter_count_sizes, c=scatter_count_colours, linewidth=0)
        else:
            plots[0].scatter(scatter_data_x, scatter_data_y, c=scatter_data_colours, s=self.point_size)

        plots[0].axis([0,max(bottom_scatter_x) + 1,min(scatter_count_y), max(scatter_count_y)])

        y_tick_positions, y_tick_labels = MicPlotter._top_plot_y_ticks(mic_data, self.antibiotic, self.log_y)
        plots[0].yaxis.set_ticks(y_tick_positions)
        plots[0].set_yticklabels(y_tick_labels)
        ylabel = r'$\log_' + str(int(self.log_y)) + '$(MIC) $\mu$g/mL' if self.log_y > 0 else r'MIC $\mu$g/mL'
        plots[0].set_ylabel(ylabel)
        plots[0].set_xticklabels([])
        plots[0].set_title(self.main_title, fontsize=18)

        # ------------------------- BOTTOM PLOT -----------------------
        plots[bottom_plot_index].axis([0,max(bottom_scatter_x) + 1,0,max(bottom_scatter_y) + 1])
        plots[bottom_plot_index].scatter(bottom_scatter_x, bottom_scatter_y, marker='o', s=self.dot_size, color=bottom_colours)
        plots[bottom_plot_index].spines["top"].set_visible(False)
        plots[bottom_plot_index].spines["right"].set_visible(False)
        plots[bottom_plot_index].spines["bottom"].set_visible(False)
        plots[bottom_plot_index].spines["left"].set_visible(False)
        plots[bottom_plot_index].yaxis.set_tick_params(length=0)
        plots[bottom_plot_index].xaxis.set_ticks([])
        plots[bottom_plot_index].set_xticklabels([])
        plots[bottom_plot_index].yaxis.set_ticks([(i+1) for i in range(len(bottom_plot_rows))])
        plots[bottom_plot_index].set_yticklabels(bottom_plot_rows[::-1])

        # ------------------------- RIGHT PLOT -------------------------
        if self.point_size == 0:
            right_x_coord = 0.75
            right_x, right_y, right_sizes = MicPlotter._right_plot_data(scatter_count_sizes, 5, right_x_coord)
            plots[1].scatter(right_x, right_y, s=right_sizes, c="black")
            plots[1].axis('off')
            plots[1].axis([0,4,-2*len(right_y),len(right_y)+1])
            for i, y in enumerate(right_y):
                plots[1].annotate(right_sizes[i], [right_x_coord + 0.75, y-0.2])
            plots[1].annotate("Counts", [right_x_coord - 0.1, len(right_y) + 0.5])

        plt.tight_layout()
        plt.savefig(self.outprefix + '.pdf')


    def run(self):
        mic_data = MicPlotter._load_mic_file(self.mic_file)
        summary_data = MicPlotter._load_summary_file(self.summary_file)
        boxplot_tsv = self.outprefix + '.boxplot.tsv'
        top_plot_data, all_mutations, combinations = MicPlotter._get_top_plot_data(summary_data, mic_data, self.antibiotic, self.use_hets, no_combinations=self.no_combinations, interrupted=self.interrupted, outfile=boxplot_tsv)
        top_plot_data, all_mutations, combinations = MicPlotter._filter_top_plot_data(top_plot_data, all_mutations, combinations, self.min_samples)
        self._make_plot(mic_data, top_plot_data, all_mutations, combinations)
