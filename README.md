ARIBA
=====

Antimicrobial Resistance Identification By Assembly

For methods and benchmarking, please see the [preprint on biorxiv][ariba biorxiv].


For how to use ARIBA, please see the [ARIBA wiki page][ARIBA wiki].



Installation
------------

ARIBA has the following dependencies, which need to be installed:
  * [Python3][python] version >= 3.3.2
  * [Bowtie2][bowtie2] version >= 2.1.0
  * [CD-HIT][cdhit] version >= 4.6
  * [MUMmer][mummer] version >= 3.23


Once the dependencies are installed, install ARIBA using pip:

    pip3 install ariba

ARIBA also depends on several Python packages, all of which are available
via pip, so the above command will get those automatically if they
are not installed. The packages are dendropy >= 4.2.0, matplotlib (no
minimum version required, but only tested on 2.0.0),
pyfastaq >= 3.12.0, pysam >= 0.9.1, and pymummer >= 0.10.1.

Alternatively, you can download the latest release from this github repository,
or clone the repository. Then run the tests:

    python3 setup.py test

If the tests all pass, install:

    python3 setup.py install

### Docker
ARIBA can be run in a Docker container. First of all install Docker, then to install ARIBA run:

    docker pull sangerpathogens/ariba

To use ARIBA you would use a command such as this (substituting in your directories), where your files are assumed to be stored in /home/ubuntu/data:

    docker run --rm -it -v /home/ubuntu/data:/data sangerpathogens/ariba ariba -h


### Debian (testing)
ARIBA is available in the latest version of Debian, and over time will progressively filter through to Ubuntu and other distributions which use Debian. To install it as root:

    sudo apt-get install ariba


### Ubuntu

You can use `apt-get` (see above), or to ensure you get the latest version of ARIBA, the following commands can be
used to install ARIBA and its dependencies. This was tested on a new instance of Ubuntu 16.04.

    sudo  apt-get update
    sudo apt-get install -y python3-dev python3-pip python3-tk zlib1g-dev bowtie2 mummer cd-hit
    export ARIBA_CDHIT=cdhit-est
    sudo pip3 install ariba

### Dependencies and environment variables

By default, ARIBA will look for the dependencies in your `$PATH`, using
the names in the table below. This behaviour can be overridden and
point ARIBA to a specific program using environment variables.
The environment variable is checked first and is used if it is set.
Otherwise ARIBA looks in your `$PATH` for the default name. This applies
to the following dependencies.

| Dependency     |  Default executable    | Environment variable name |
|----------------|------------------------|---------------------------|
| Bowtie2        | `bowtie2`              | `$ARIBA_BOWTIE2`          |
| CD-HIT (est)   | `cd-hit-est`           | `$ARIBA_CDHIT`            |


For example, you could specify an exact version of a bowtie2 executable
that you compiled and downloaded in your home directory (assuming BASH):

    export ARIBA_BOWTIE2=$HOME/bowtie2-2.1.0/bowtie2

Note that ARIBA also runs `bowtie2-build`, for which it uses the
`bowtie2` executable with `-build` appended. So in this case
it would try to use

    $HOME/bowtie2-2.1.0/bowtie2-build

 
### Temporary files
 

ARIBA can temporarily make a large number of files whilst running, which
are put in a temporary directory made by ARIBA.  The total size of these
files is small, but there can be a many of them. This can be a
problem when running large numbers (100s or 1000s) of jobs simultaneously
on the same file system.
The parent directory of the temporary directory is determined in the
following order of precedence:

1. The value of the option `--tmp_dir` (if that option was used)
2. The environment variable `$ARIBA_TMPDIR` (if it is set)
3. The environment variable `$TMPDIR` (if it is set)
4. If none of the above is found, then use the run's output directory.

Each temporary directory
is unique to one run of ARIBA, and is automatically deleted at the end
of the run (even if ARIBA was killed by the user or crashed).
For example,

    export $ARIBA_TMPDIR=/tmp

will result in the creation of a new directory inside `/tmp`, which
will have a name of the form

    /tmp/ariba.tmp.abcdef

where the suffix `abcdef` is a random string of characters, chosen
such that `/tmp/ariba.tmp.abcdef` does not already exist.

The exception to the above is if the option `--noclean` is used.
This forces the temporary directory to be placed in the output
directory, and temporary files are kept. It is intended for
debugging.



Usage
-----

Please read the [ARIBA wiki page][ARIBA wiki] for usage instructions.



Build status: [![Build Status](https://travis-ci.org/sanger-pathogens/ariba.svg?branch=master)](https://travis-ci.org/sanger-pathogens/ariba)

  [ariba biorxiv]: http://biorxiv.org/content/early/2017/04/07/118000
  [bowtie2]: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
  [cdhit]: http://weizhongli-lab.org/cd-hit/
  [ARIBA wiki]: https://github.com/sanger-pathogens/ariba/wiki
  [mummer]: http://mummer.sourceforge.net/
  [python]: https://www.python.org/


