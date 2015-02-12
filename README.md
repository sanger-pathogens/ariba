ARIBA
=====

Antibiotic Resistance Identification By Assembly


Installation
------------

Dependencies:
  * [samtools and bcftools] [samtools]  version >= 1.2
  * [SSPACE-basic scaffolder] [sspace]
  * [GapFiller] [gapfiller]
  * [SMALT] [smalt]

Install with pip:

    pip3 install pyfastaq

Alternatively, you can download the latest release from this github repository,
or clone the repository. Then run the tests:

    python3 setup.py test

If the tests all pass, install:

    python3 setup.py install



Usage
-----

The installation will put a single script called `ariba` in your path.
The usage is:

    ariba <command> [options]

Run just `ariba` to get a list of tasks. Use `-h` or `--help`
with a task to get the help, for example

    ariba run --help

To run the pipeline to identify the genes present in a pair of FASTQ files:

    ariba run reference_genes.fa reads_1.fastq reads_2.fastq Output_directory



  [smalt]: https://www.sanger.ac.uk/resources/software/smalt/
  [sspace]: http://www.baseclear.com/genomics/bioinformatics/basetools/SSPACE
  [gapfiller]: http://www.baseclear.com/genomics/bioinformatics/basetools/gapfiller
  [samtools]: http://www.htslib.org/
