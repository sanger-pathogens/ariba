import re
import pyfastaq

class Error (Exception): pass

name_regex = re.compile(r'^(?P<vfdb_id>V\S+\)) \((?P<name>\S+)\) (?P<description>.*\]) \[(?P<genus_etc>.*)\]$')

class VfdbParser:
    def __init__(self, infile, outprefix):
        self.infile = infile
        self.outprefix = outprefix


    @classmethod
    def _fa_header_to_name_pieces(cls, fa_header):
        m = name_regex.search(fa_header)
        if m is None:
            return None
        else:
            return tuple([m.group(x) for x in ['vfdb_id', 'name', 'description', 'genus_etc']])


    @staticmethod
    def _fa_header_to_name_and_metadata(fa_header):
        name_data = VfdbParser._fa_header_to_name_pieces(fa_header)
        if name_data is None:
            return fa_header, '.'
        else:
            vfdb_id, name, description, genus_etc = name_data
            return name + '.' + vfdb_id + '.' + genus_etc.replace(' ', '_'), description


    def run(self):
        file_reader = pyfastaq.sequences.file_reader(self.infile)
        fa_out = pyfastaq.utils.open_file_write(self.outprefix + '.fa')
        tsv_out = pyfastaq.utils.open_file_write(self.outprefix + '.tsv')

        for seq in file_reader:
            original_id = seq.id
            seq.id, description = self._fa_header_to_name_and_metadata(seq.id)
            if description == '.':
                seq.id = original_id.split()[0]
            print(seq.id, '1', '0', '.', '.', 'Original name: ' + original_id, sep='\t', file=tsv_out)
            print(seq, file=fa_out)

        pyfastaq.utils.close(fa_out)
        pyfastaq.utils.close(tsv_out)

