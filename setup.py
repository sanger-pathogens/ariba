import os
import shutil
import sys
import glob
from setuptools import setup, find_packages


required_progs = [
    'bcftools',
    'GapFiller.pl',
    'SSPACE_Basic_v2.0.pl',
    'samtools',
    'smalt',
]

print('Checking that dependencies are found in path...')
found_all_progs = True
for program in required_progs:
    if shutil.which(program) is None:
        found_all_progs = False
        found = ' NOT FOUND'
    else:
        found = ' OK'
    print(found, program, sep='\t')
if not found_all_progs:
    print('Cannot install because some required programs not found.', file=sys.stderr)
    sys.exit(1)


if not(  (shutil.which('velveth') and shutil.which('velvetg')) or shutil.which('spades.py')  ):
    print('Must have velvet (velveth and velvetg) or spades.py installed. Cannot continue', file=sys.stderr)
    sys.exit(1)


print('... all dependencies found')


setup(
    name='ariba',
    version='0.0.1',
    description='ARIBA: Antibiotic Resistance Identification By Assembly',
    packages = find_packages(),
    author='Martin Hunt',
    author_email='path-help@sanger.ac.uk',
    url='https://github.com/sanger-pathogens/ariba',
    scripts=glob.glob('scripts/*'),
    test_suite='nose.collector',
    install_requires=[
        'nose >= 1.3',
        'openpyxl',
        'pyfastaq >= 3.0.1',
        'pysam >= 0.8.1',
        'pymummer>=0.0.2'
    ],
    license='GPLv3',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Programming Language :: Python :: 3 :: Only',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
    ],
)
