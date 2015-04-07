ARIBA
=====

Antibiotic Resistance Identification By Assembly

For how to use ARIBA, please see the [ARIBA wiki page] [ARIBA wiki].



Installation
------------

ARIBA has the following dependencies, which need to be installed:
  * [bowtie2] [Bowtie2] version >= 2.1.0
  * [cd-hit] [cdhit] version >= 4.6
  * [samtools and bcftools] [samtools]  version >= 1.2
  * [MUMmer] [mummer] version >= 3.23
  * Either [SPAdes] [spades] version >= 3.5.0 or [Velvet] [velvet] version >= 1.2.07
    (SPAdes is recommended)

ARIBA has the following optional dependencies. If they are installed,
they will be used. Otherwise scaffolding and gap filling will be
skipped.
  * [SSPACE-basic scaffolder] [sspace]
  * [GapFiller] [gapfiller]

Once the dependencies are installed, install ARIBA using pip:

    pip3 install ariba

Alternatively, you can download the latest release from this github repository,
or clone the repository. Then run the tests:

    python3 setup.py test

If the tests all pass, install:

    python3 setup.py install


Usage
-----

Please read the [ARIBA wiki page] [ARIBA wiki] for usage instructions.


  [bowtie2]: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
  [cdhit]: http://weizhongli-lab.org/cd-hit/
  [ARIBA wiki]: https://github.com/sanger-pathogens/ariba/wiki
  [gapfiller]: http://www.baseclear.com/genomics/bioinformatics/basetools/gapfiller
  [mummer]: http://mummer.sourceforge.net/
  [samtools]: http://www.htslib.org/
  [spades]: http://bioinf.spbau.ru/spades
  [sspace]: http://www.baseclear.com/genomics/bioinformatics/basetools/SSPACE
  [velvet]: http://www.ebi.ac.uk/~zerbino/velvet/


