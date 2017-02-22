import csv
import re
import os
from ariba import common

class Error (Exception): pass

regex_string_to_float = re.compile(r'\s*(?P<lt_or_gt>[<>]?)\s*(?P<equals>=?)\s*(?P<number>[0-9.]+)\s*$')

class MicPlotter:
    def __init__(self,
      antibiotic,
      mic_file,
      summary_file,
      outprefix,
      main_title=None,
      plot_height=15,
      plot_width=15,
      log_y=True,
      plot_types="points,violin",
      jitter_width=0.1,
      jitter_height=0.01,
      no_combinations=False,
      mic_values='0,0.001,0.0025,0.0075,0.015,0.03,0.06,0.125,0.25,0.5,1,2,4,8,16,32,64,128,256,512,1024',
      hlines='0.25,2'
    ):
        self.antibiotic = antibiotic
        self.mic_file = mic_file
        self.summary_file = summary_file
        self.outprefix = outprefix
        self.main_title = self.antibiotic if main_title is None else main_title
        self.plot_height = plot_height
        self.plot_width = plot_width
        self.log_y = log_y
        self.plot_types = set(plot_types.split(','))

        allowed_plot_types = {'point', 'violin', 'boxplot'}
        if not self.plot_types.issubset(allowed_plot_types):
            raise Error('Error in plot_types option. Allowed types are: ' + str(allowed_plot_types) + '. Got: ' +  str(self.plot_types))

        self.jitter_width = jitter_width
        self.jitter_height = jitter_height
        self.no_combinations = no_combinations

        try:
            self.mic_values = [float(x) for x in mic_values.split(',')]
        except:
            raise Error('Error in mic_values option. Needs to be a list of numbers separated by commas. Got this:\n' + mic_values)

        try:
            if len(hlines) == 0:
                self.hlines = []
            else:
                self.hlines = [float(x) for x in hlines.split(',')]
        except:
            raise Error('Error in hlines option. Needs to be a list of numbers separated by commas, or empty. Got this:\n' + hlines)


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
    def _to_boxplot_tsv(cls, summary_data, mic_data, antibiotic, outfile, no_combinations=False):
        ignore_columns = {'assembled', 'match', 'ref_seq', 'pct_id', 'known_var', 'novel_var'}
        all_mutations = set()
        all_mutations_seen_combinations = set()

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

                if len(mutations) == 0:
                    mutations.add('without_mutation')

                all_mutations.update(mutations)
                mutations = list(mutations)
                mutations.sort()
                if no_combinations:
                    for mutation in mutations:
                        all_mutations_seen_combinations.add((mutation,))
                        print(sample, mic_data[sample][antibiotic], mutation, sep='\t', file=f)
                else:
                    all_mutations_seen_combinations.add(tuple(mutations))
                    mutations = '.'.join(mutations)
                    print(sample, mic_data[sample][antibiotic], mutations, sep='\t', file=f)

        return all_mutations, all_mutations_seen_combinations


    @classmethod
    def _to_dots_tsv(cls, all_mutations, combinations, outfile):
        if 'without_mutation' in all_mutations:
            all_mutations.remove('without_mutation')
            combinations.remove(('without_mutation',))
            has_without_mutation = True
        else:
            has_without_mutation = False

        all_mutations = list(all_mutations)
        all_mutations.sort()
        combinations = list(combinations)
        combinations.sort()

        if has_without_mutation:
            all_mutations.append('without_mutation')
            combinations.append(('without_mutation',))

        output_columns = {}
        for combination in combinations:
            output_columns[combination] = [(1 if x in combination else 0) for x in all_mutations]

        with open(outfile, 'w') as f:
            print('Mutation', end='', file=f)
            for x in combinations:
                print('\t', '.'.join(x), sep='', end='', file=f)
            print('', file=f)

            for i in range(len(all_mutations)):
                row = [all_mutations[i]] + [output_columns[x][i] for x in combinations]
                print(*row, sep='\t', file=f)


    def _make_plot(self,
      samples_file,
      dots_file,
    ):
        r_script = self.outprefix + '.R'

        try:
            f = open(r_script, 'w')
        except:
            raise Error('Error opening R script for writing "' + r_script + '"')

        libraries = ['ggplot2', 'RColorBrewer', 'reshape2', 'cowplot', 'latex2exp']
        for lib in libraries:
            print('library(', lib, ')', sep='', file=f)

        print('samples = read.csv(file="', samples_file, r'''", header=TRUE, sep="\t")''', sep='', file=f)
        print('dots = read.csv(file="', dots_file, r'''", header=TRUE, sep="\t", check.names=FALSE)''', sep='', file=f)

        if self.log_y:
            print('use.log = TRUE', file=f)
        else:
            print('use.log = FALSE', file=f)

        print(r'''
dots.melt = melt(dots)
colnames(dots.melt) <- c("var1", "var2", "value")

accent <- brewer.pal(8, 'Accent')
        accentPalette <- colorRampPalette(accent)
        ncols <- length(as.vector(unique(samples$Mutations)))
        cols <- accentPalette(ncols)


names(cols) <- sort(as.vector(unique(samples$Mutations)))
setcols <- c()

for (i in 1:nrow(dots.melt)){
    if (dots.melt[i,3]==1){
        setcols <- c(setcols, cols[as.vector(dots.melt[i,2])])
    }
    else{
      setcols <- c(setcols, NA)
    }
}
dots.melt <- cbind(dots.melt, setcols)

mutations <- levels(dots.melt$var1)
i <- match("without_mutation", mutations)
if (!is.na(i)) {
    mutations <- c(mutations[-i], "without_mutation")
}


dotplot <- ggplot(dots.melt, aes(x=var2, y=var1)) +
  geom_point(aes(fill=setcols, colour=setcols), size=8) +
          scale_fill_identity()+
          scale_colour_identity()+
          ylim(rev(mutations)) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=18),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position="none")

range.mics <- c(''' + ','.join([str(x) for x in self.mic_values]) + r''')
if (use.log){ final.mics <- log(range.mics) }else{ final.mics <- range.mics }
''', file=f)

        if self.log_y:
            print(r'''violinplot <- ggplot(data=samples, aes(x=Mutations, y=log(MIC))) +''', file=f)
        else:
            print(r'''violinplot <- ggplot(data=samples, aes(x=Mutations, y=MIC)) +''', file=f)

        if 'point' in self.plot_types:
            print(r'''    geom_point(aes(color=Mutations), position = position_jitter(width=''', self.jitter_width, ', height=', self.jitter_height, '), size=4, alpha=.5) +', sep='', file=f)

        if 'violin' in self.plot_types:
            print(r'''    geom_violin(aes(color=Mutations),alpha=.10, show.legend = FALSE) +''', file=f)

        if 'boxplot' in self.plot_types:
            print(r'''    geom_boxplot(aes(color=Mutations),alpha=.10, show.legend = FALSE) +''', file=f)

        if self.no_combinations:
            axis_text_x = 'element_text(size=24, angle=45, hjust=1)'
        else:
            axis_text_x = 'element_blank()'

        for x in self.hlines:
            if self.log_y:
                print('    geom_hline(yintercept=log(', x, '), lty=2) +', sep='', file=f)
            else:
                print('    geom_hline(yintercept=', x, ', lty=2) +', sep='', file=f)

        print(r'''    ylab(expression(paste("MIC ", mu, "g/mL"))) +
    scale_colour_manual(values = cols) +
    ggtitle("''' + self.main_title + r'''") +
    scale_y_continuous(breaks=final.mics, labels=range.mics) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            axis.line.x = element_line(color="black"),
            axis.line.y = element_line(color="black"),
            axis.title.x = element_blank(),
            axis.title.y = element_text(size=22),
            axis.text.x = ''' + axis_text_x + r''',
            axis.text.y = element_text(size=24),
            axis.title = element_text(size=20),
            plot.title = element_text(lineheight=.6, size = 24, hjust=.5, face="bold"),
            legend.position="none")
''', file=f)

        if self.no_combinations:
            print('violinplot', file=f)
        else:
            print('plot_grid(violinplot, dotplot, ncol=1, align="v", rel_heights=c(3,1))', file=f)

        print('ggsave("', self.outprefix, '.pdf", useDingbats=FALSE, height=', self.plot_height, ', width=', self.plot_width, ')', sep='', file=f)
        f.close()
        common.syscall('R CMD BATCH ' + r_script)


    def run(self):
        mic_data = MicPlotter._load_mic_file(self.mic_file)
        summary_data = MicPlotter._load_summary_file(self.summary_file)
        boxplot_tsv = self.outprefix + '.boxplot.tsv'
        all_mutations, combinations = MicPlotter._to_boxplot_tsv(summary_data, mic_data, self.antibiotic, boxplot_tsv, no_combinations=self.no_combinations)
        dots_tsv = self.outprefix + '.dots.tsv'
        MicPlotter._to_dots_tsv(all_mutations, combinations, dots_tsv)
        self._make_plot(boxplot_tsv, dots_tsv)
        
