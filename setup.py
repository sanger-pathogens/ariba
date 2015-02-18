import os
import shutil
import sys
import glob
from setuptools import setup, find_packages


setup(
    name='ariba',
    version='0.1.2',
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
