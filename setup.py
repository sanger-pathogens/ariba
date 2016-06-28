import os
import shutil
import sys
import glob
from setuptools import setup, find_packages


setup(
    name='ariba',
    version='1.0.1',
    description='ARIBA: Antibiotic Resistance Identification By Assembly',
    packages = find_packages(),
    package_data={'ariba': ['test_run_data/*']},
    author='Martin Hunt',
    author_email='path-help@sanger.ac.uk',
    url='https://github.com/sanger-pathogens/ariba',
    scripts=glob.glob('scripts/*'),
    test_suite='nose.collector',
    tests_require=['nose >= 1.3'],
    install_requires=[
        'openpyxl >= 1.6.2',
        'pyfastaq >= 3.12.0',
        'pysam >= 0.8.1, <= 0.8.3',
        'pymummer>=0.6.1',
    ],
    license='GPLv3',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Programming Language :: Python :: 3 :: Only',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
    ],
)
