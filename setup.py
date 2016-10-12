import os
import shutil
import sys
import glob
from setuptools import setup, find_packages, Extension

minimap_c_files = [
    'bseq.c',
    'index.c',
    'kthread.c',
    'map.c',
    'misc.c',
    'sdust.c',
    'sketch.c',
]

minimap_c_files = [os.path.join('third_party', 'minimap-0.2', x) for x in minimap_c_files]
minimap_c_files.append(os.path.join('ariba', 'ext', 'minimap_ariba.cpp'))
minimap_mod = Extension(
    "minimap_ariba",
    minimap_c_files,
    extra_link_args=['-lz'],
    include_dirs=[os.path.join('third_party', 'minimap-0.2')],
)

fermilite_c_files = [
    'bfc.c',
    'bseq.c',
    'bubble.c',
    'htab.c',
    'ksw.c',
    'kthread.c',
    'mag.c',
    'misc.c',
    'mrope.c',
    'rld0.c',
    'rle.c',
    'rope.c',
    'unitig.c'
]
fermilite_c_files = [os.path.join('third_party', 'fermi-lite-0.1', x) for x in fermilite_c_files]
fermilite_c_files.append(os.path.join('ariba', 'ext', 'fml-asm_ariba.cpp'))
fermilite_mod = Extension(
    "fermilite_ariba",
    fermilite_c_files,
    extra_link_args=['-lz'],
    include_dirs=[os.path.join('third_party', 'fermi-lite-0.1')],
)

vcfcall_mod = Extension(
    "vcfcall_ariba",
    [os.path.join('ariba', 'ext', 'vcfcall_ariba.cpp')],
)

setup(
    ext_modules=[minimap_mod, fermilite_mod, vcfcall_mod],
    name='ariba',
    version='2.3.0',
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
        'dendropy >= 4.1.0',
        'pyfastaq >= 3.12.0',
        'pysam >= 0.9.1',
        'pymummer>=0.8.1',
    ],
    license='GPLv3',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Programming Language :: Python :: 3 :: Only',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
    ],
)
