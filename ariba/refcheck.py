import os
import pyfastaq

class Error (Exception): pass


class Checker:
    def __init__(self, infile, min_length=1, max_length=10000, outprefix=None):
        self.infile = os.path.abspath(infile)
        if not os.path.exists(self.infile):
            raise Error('File not found: "' + self.infile + '". Cannot continue')

        self.min_length = min_length
        self.max_length = max_length
        self.outprefix = outprefix


    def run(self):
        file_reader = pyfastaq.sequences.file_reader(self.infile)
        names = {}

        if self.outprefix is not None:
            old2new_out = self.outprefix + '.rename'
            fasta_out = self.outprefix + '.fa'
            bad_seqs_out = self.outprefix + '.removed.fa'
            log_out = self.outprefix + '.log'
            old2new_out_fh = pyfastaq.utils.open_file_write(old2new_out)
            fasta_out_fh = pyfastaq.utils.open_file_write(fasta_out)
            bad_seqs_out_fh = pyfastaq.utils.open_file_write(bad_seqs_out)
            log_out_fh = pyfastaq.utils.open_file_write(log_out)

        for seq in file_reader:
            seq.seq = seq.seq.upper()
            if len(seq) < self.min_length:
                if self.outprefix is None:
                    return False, 'Too short', seq
                else:
                    print(seq.id, 'Too short. Skipping', sep='\t', file=log_out_fh)
                    print(seq, file=bad_seqs_out_fh)
                    continue
            elif len(seq) > self.max_length:
                if self.outprefix is None:
                    return False, 'Too long', seq
                else:
                    print(seq.id, 'Too long. Skipping', sep='\t', file=log_out_fh)
                    print(seq, file=bad_seqs_out_fh)
                    continue

            if not seq.looks_like_gene():
                if self.outprefix is None:
                    return False, 'Not a gene', seq
                else:
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

            if self.outprefix is None and original_id != seq.id:
                seq.id = original_id
                return False, 'Name has spaces', seq

            if seq.id in names:
                if self.outprefix is None:
                    return False, 'Duplicate name', seq
                else:
                    names[seq.id] += 1
                    seq.id += '.' + str(names[seq.id])
            else:
                names[seq.id] = 1

            if self.outprefix is not None:
                print(original_id, seq.id, sep='\t', file=old2new_out_fh)
                print(seq, file=fasta_out_fh)

        if self.outprefix is not None:
            pyfastaq.utils.close(fasta_out_fh)
            pyfastaq.utils.close(bad_seqs_out_fh)
            pyfastaq.utils.close(log_out_fh)
            pyfastaq.utils.close(old2new_out_fh)

        return True, None, None
