import os
import pyfastaq

class Error (Exception): pass


class Checker:
    def __init__(self, infile, min_length=1, max_length=10000):
        self.infile = os.path.abspath(infile)
        if not os.path.exists(self.infile):
            raise Error('File not found: "' + self.infile + '". Cannot continue')

        self.min_length = min_length
        self.max_length = max_length


    def check(self, error_code_on_exit=None):
        file_reader = pyfastaq.sequences.file_reader(self.infile)

        for seq in file_reader:
            if not seq.looks_like_gene():
                return False, 'Not a gene', seq
            elif len(seq) < self.min_length:
                return False, 'Too short', seq
            elif len(seq) > self.max_length:
                return False, 'Too long', seq
           
        return True, None, None


    def fix(self, outprefix):
        file_reader = pyfastaq.sequences.file_reader(self.infile)
        old2new_out = outprefix + '.rename'
        fasta_out = outprefix + '.fa'
        bad_seqs_out = outprefix + '.removed.fa'
        log_out = outprefix + '.log'
        names = {}
        old2new_out_fh = pyfastaq.utils.open_file_write(old2new_out)
        fasta_out_fh = pyfastaq.utils.open_file_write(fasta_out)
        bad_seqs_out_fh = pyfastaq.utils.open_file_write(bad_seqs_out)
        log_out_fh = pyfastaq.utils.open_file_write(log_out)

        for seq in file_reader:
            seq.seq = seq.seq.upper()
            if len(seq) < self.min_length:
                print(seq.id, 'Too short. Skipping', sep='\t', file=log_out_fh)
                print(seq, file=bad_seqs_out_fh)
                continue
            elif len(seq) > self.max_length:
                print(seq.id, 'Too long. Skipping', sep='\t', file=log_out_fh)
                print(seq, file=bad_seqs_out_fh)
                continue


            if not seq.looks_like_gene():
                seq.revcomp()
                if seq.looks_like_gene():
                    print(seq.id, 'Reverse complemented', sep='\t', file=log_out_fh)
                else:
                    print(seq.id, 'Does not look like a gene. Skipping', sep='\t', file=log_out_fh)
                    seq.revcomp()
                    print(seq, file=bad_seqs_out_fh)
                    continue
                
            original_id = seq.id
            # replace unwanted characters with underscores
            to_replace = ' '
            seq.id = seq.id.translate(str.maketrans(to_replace, '_' *  len(to_replace)))

            if seq.id in names:
                names[seq.id] += 1
                seq.id += '.' + str(names[seq.id])
            else:
                names[seq.id] = 1

            print(original_id, seq.id, sep='\t', file=old2new_out_fh)
            print(seq, file=fasta_out_fh)

        pyfastaq.utils.close(fasta_out_fh)
        pyfastaq.utils.close(bad_seqs_out_fh)
        pyfastaq.utils.close(log_out_fh)
        pyfastaq.utils.close(old2new_out_fh)
