import argparse
import ariba

def run():
    parser = argparse.ArgumentParser(
        description = 'Translate the meaning of a flag output by ARIBA',
        usage = 'ariba flag <flag>')
    parser.add_argument('flag_in', type=int, help='Flag to be translated (an integer)', metavar='flag')
    options = parser.parse_args()

    f = ariba.flag.Flag(options.flag_in)
    print('Meaning of flag', f)
    print(f.to_long_string())
