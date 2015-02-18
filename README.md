ARIBA
=====

Antibiotic Resistance Identification By Assembly

For how to use ARIBA, please see the [ARIBA wiki page] [ARIBA wiki].



Installation
------------

ARIBA has the following dependencies, which need to be installed:
  * [samtools and bcftools] [samtools]  version >= 1.2
  * [SSPACE-basic scaffolder] [sspace]
  * [GapFiller] [gapfiller]
  * [MUMmer] [mummer] version >= 3.23
  * [SMALT] [smalt] version >= 0.7.4
  * Either [SPAdes] [spades] version >= 3.5.0 or [Velvet] [velvet] version >= 1.2.07
    (SPAdes is recommended)

Once the dependencies are installed, install ARIBA using pip:

    pip3 install pyfastaq

Alternatively, you can download the latest release from this github repository,
or clone the repository. Then run the tests:

    python3 setup.py test

If the tests all pass, install:

    python3 setup.py install


Usage
-----

Please read the [ARIBA wiki page] [ARIBA wiki] for usage instructions.


  [ARIBA wiki]: https://github.com/sanger-pathogens/ariba/wiki
  [gapfiller]: http://www.baseclear.com/genomics/bioinformatics/basetools/gapfiller
  [mummer]: http://mummer.sourceforge.net/
  [samtools]: http://www.htslib.org/
  [smalt]: https://www.sanger.ac.uk/resources/software/smalt/
  [spades]: http://bioinf.spbau.ru/spades
  [sspace]: http://www.baseclear.com/genomics/bioinformatics/basetools/SSPACE
  [velvet]: http://www.ebi.ac.uk/~zerbino/velvet/


