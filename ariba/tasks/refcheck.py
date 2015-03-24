import argparse
import sys
import pyfastaq
import ariba

def run():
    parser = argparse.ArgumentParser(
        description = 'Check or fix resistance genes FASTA file',
        usage = 'ariba refcheck [options] <infile>')
    parser.add_argument('--genetic_code', type=int, help='Number of genetic code to use. Currently supported 1,4,11 [%(default)s]', choices=[1,4,11], default=11, metavar='INT')
    parser.add_argument('-m', '--min_length', type=int, help='Minimum length in nucleotides of gene [%(default)s]', metavar='INT', default=6)
    parser.add_argument('-n', '--max_length', type=int, help='Maximum length in nucleotides of gene [%(default)s]', metavar='INT', default=10000)
    parser.add_argument('-o', '--outprefix', help='Prefix of output files. If this option is used, a fixed file will be output, together with information on what was changed in the input file. If this option is not used, the script dies if any input sequence is not OK')
    parser.add_argument('infile', help='Input file containing genes to be checked', metavar='Filename')
    options = parser.parse_args()

    pyfastaq.sequences.genetic_code = options.genetic_code
    checker = ariba.refcheck.Checker(
        options.infile,
        min_length=options.min_length,
        max_length=options.max_length
    )
 
    if options.outprefix:
        checker.fix(options.outprefix)
    else:
        ok, reason, seq = checker.check()
        if not ok:
            print('The following sequence not OK, for the reason:', reason)
            print(seq)
            sys.exit(1)
