import argparse
import ariba

def run(options):
    f = ariba.flag.Flag(options.flag_in)
    print('Meaning of flag', f)
    print(f.to_long_string())
