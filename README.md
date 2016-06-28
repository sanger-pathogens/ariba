Disclaimer
==========

This is pre-publication software that is currently under active development.
Use it at your own risk. Bug reports are welcome, but
user support is not provided at this time.


ARIBA
=====

Antibiotic Resistance Identification By Assembly


For how to use ARIBA, please see the [ARIBA wiki page][ARIBA wiki].



Installation
------------

ARIBA has the following dependencies, which need to be installed:
  * [Python3][python] version >= 3.4
  * [R][r] version >= 2.14.0
  * The R package [ape][ape] version >= 3.1
  * [Bowtie2][bowtie2] version >= 2.1.0
  * [CD-HIT][cdhit] version >= 4.6
  * [Samtools and BCFtools][samtools]  version >= 1.2
  * [MUMmer][mummer] version >= 3.23
  * [SPAdes][spades] version >= 3.5.0
  * [Python2][python] version >= 2.7 (SPAdes needs python2)


ARIBA has the following optional dependencies. If they are installed,
they will be used. Otherwise scaffolding and gap filling will be
skipped.
  * [SSPACE-basic scaffolder][sspace]
  * [GapFiller][gapfiller]

Once the dependencies are installed, install ARIBA using pip:

    pip3 install ariba

Alternatively, you can download the latest release from this github repository,
or clone the repository. Then run the tests:

    python3 setup.py test

If the tests all pass, install:

    python3 setup.py install


### Dependencies and environment variables

By default, ARIBA will look for the dependencies in your `$PATH`, using
the names in the table below. This behaviour can be overridden and
point ARIBA to a specific program using environment variables.
The environment variable is checked first and is used if it is set.
Otherwise ARIBA looks in your `$PATH` for the default name. This applies
to the following dependencies.

| Dependency     |  Default executable    | Environment variable name |
|----------------|------------------------|---------------------------|
| BCFtools       | `bcftools`             | `$ARIBA_BCFTOOLS`         |
| Bowtie2        | `bowtie2`              | `$ARIBA_BOWTIE2`          |
| CD-HIT         | `cd-hit-est`           | `$ARIBA_CDHIT`            |
| GapFiller      | `GapFiller.pl`         | `$ARIBA_GAPFILLER`        |
| R              | `Rscript`              | `$ARIBA_R`                |
| Samtools       | `samtools`             | `$ARIBA_SAMTOOLS`         |
| SPAdes         | `spades.py`            | `$ARIBA_SPADES`           |
| SSPACE         | `SSPACE_Basic_v2.0.pl` | `$ARIBA_SSPACE`           |


For example, you could specify an exact version of a Samtools executable
that you compiled and downloaded in your home directory (assuming BASH):

    export ARIBA_SAMTOOLS=$HOME/samtools-1.2/samtools

### Temporary files

ARIBA can temporarily make a large number of files whilst running, which
are put in a temporary directory made by ARIBA.  The total size of these
files is small, but there can be a many of them. This can be a
problem when running large numbers (100s or 1000s) of jobs simultaneously
on the same file system.
By default, ARIBA creates a temporary directory for these files
inside the output directory of each run.

Each temporary directory
is unique to one run of ARIBA, and is automatically deleted at the end
of the run (even if ARIBA was killed by the user or crashed).
The parent directory of the temporary
directory can be changed using the environment variable
`$ARIBA_TMPDIR`. The temporary directory for each run will be made
inside `$ARIBA_TMPDIR`. For example,

    export $ARIBA_TMPDIR=/tmp

will result in the creation of a new directory inside `/tmp`, which
will have a name of the form

    /tmp/ariba.tmp.abcdef

where the suffix `abcdef` is a random string of characters, chosen
such that `/tmp/ariba.tmp.abcdef` does not already exist.

The temporary directory can also be changed using the option
`--tmp_dir` when running `ariba run`. Using this option takes precedence
over the environment variable `$ARIBA_TMPDIR`. If neither are
set, then ARIBA creates the temporary directory inside
the output directory given to `ariba run`.

The exception to the above is if the option `--noclean` is used.
This forces the temporary directory to be placed in the output
directory, and temporary files are kept. It is intended for
debugging.



Usage
-----

Please read the [ARIBA wiki page][ARIBA wiki] for usage instructions.



Build status: [![Build Status](https://travis-ci.org/sanger-pathogens/ariba.svg?branch=master)](https://travis-ci.org/sanger-pathogens/ariba)


  [bowtie2]: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
  [cdhit]: http://weizhongli-lab.org/cd-hit/
  [ARIBA wiki]: https://github.com/sanger-pathogens/ariba/wiki
  [gapfiller]: http://www.baseclear.com/genomics/bioinformatics/basetools/gapfiller
  [mummer]: http://mummer.sourceforge.net/
  [samtools]: http://www.htslib.org/
  [spades]: http://bioinf.spbau.ru/spades
  [sspace]: http://www.baseclear.com/genomics/bioinformatics/basetools/SSPACE
  [ape]: https://cran.r-project.org/web/packages/ape/index.html
  [r]: https://www.r-project.org/
  [python]: https://www.python.org/


