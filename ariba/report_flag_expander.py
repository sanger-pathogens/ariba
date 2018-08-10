import pyfastaq

from ariba import flag

class Error (Exception): pass

class ReportFlagExpander:
    def __init__(self, infile, outfile):
        self.infile = infile
        self.outfile = outfile


    def run(self):
        f_in = pyfastaq.utils.open_file_read(self.infile)
        f_out = pyfastaq.utils.open_file_write(self.outfile)
        flag_index = None

        for line in f_in:
            fields = line.rstrip().split()

            if flag_index is None:
                try:
                    flag_index = fields.index('flag')
                except:
                    raise Error('"flag" column not found in first line of file ' + self.infile +'. Cannot continue')
            else:
                f = flag.Flag(int(fields[flag_index]))
                fields[flag_index] = f.to_comma_separated_string()

            print(*fields, sep='\t', file=f_out)

        f_in.close()
        f_out.close()

