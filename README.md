# ARIBA

Antimicrobial Resistance Identification By Assembly

For how to use ARIBA, please see the [ARIBA wiki page][ARIBA wiki].

PLEASE NOTE: we currently do not have the resources to provide support for Ariba - see the [Feedback/Issues](#feedbackissues) section.

[![Unmaintained](http://unmaintained.tech/badge.svg)](http://unmaintained.tech/)  
[![Build Status](https://travis-ci.org/sanger-pathogens/ariba.svg?branch=master)](https://travis-ci.org/sanger-pathogens/ariba)   
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-brightgreen.svg)](https://github.com/sanger-pathogens/ariba/blob/master/LICENSE)   
[![status](https://img.shields.io/badge/MGEN-10.1099%2Fmgen.0.000131-brightgreen.svg)](http://mgen.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000131)   
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/recipes/ariba/README.html)  
[![Container ready](https://img.shields.io/badge/container-ready-brightgreen.svg)](https://quay.io/repository/biocontainers/ariba)  
[![Docker Build Status](https://img.shields.io/docker/build/sangerpathogens/ariba.svg)](https://hub.docker.com/r/sangerpathogens/ariba)  
[![Docker Pulls](https://img.shields.io/docker/pulls/sangerpathogens/ariba.svg)](https://hub.docker.com/r/sangerpathogens/ariba)  
[![codecov](https://codecov.io/gh/sanger-pathogens/ariba/branch/master/graph/badge.svg)](https://codecov.io/gh/sanger-pathogens/ariba)

## Contents
* [Introduction](#introduction)
* [Quick Start](#quick-start)
* [Installation](#installation)
  * [Required dependencies](#required-dependencies)
  * [Using pip3](#using-pip3)
  * [From Source](#from-source)
  * [Docker](#docker)
  * [Debian (testing)](#debian-testing)
  * [Ubuntu](#ubuntu)
  * [Dependencies and environment variables](#dependencies-and-environment-variables)
* [Temporary files](#temporary-files)
* [Usage](#usage)
* [License](#license)
* [Feedback/Issues](#feedbackissues)
* [Citation](#citation)

## Introduction
ARIBA is a tool that identifies antibiotic resistance genes by running local assemblies.
It can also be used for [MLST calling](https://github.com/sanger-pathogens/ariba/wiki/MLST-calling-with-ARIBA).

The input is a FASTA file of reference sequences (can be a mix of genes and noncoding sequences) and paired sequencing reads. ARIBA reports which of the reference sequences were found, plus detailed information on the quality of the assemblies and any variants between the sequencing reads and the reference sequences.

## Quick Start
Get reference data, for instance from [CARD](https://card.mcmaster.ca/). See [getref](https://github.com/sanger-pathogens/ariba/wiki/Task%3A-getref) for a full list.

    ariba getref ncbi out.ncbi

Prepare reference data for ARIBA:

    ariba prepareref -f out.ncbi.fa -m out.ncbi.tsv out.ncbi.prepareref

Run local assemblies and call variants:

    ariba run out.ncbi.prepareref reads1.fastq reads2.fastq out.run

Summarise data from several runs:

    ariba summary out.summary out.run1/report1.tsv out.run2/report2.tsv out.run3/report3.tsv

Please read the [ARIBA wiki page][ARIBA wiki] for full usage instructions.

## Tutorials
[The Jupyter notebook tutorial](https://github.com/sanger-pathogens/pathogen-informatics-training)

## Installation

If you encounter an issue when installing ARIBA please contact your local system administrator. If you encounter a bug you can log it [here](https://github.com/sanger-pathogens/ariba/issues).

### Required dependencies
  * [Python3][python] version >= 3.6.0
  * [Bowtie2][bowtie2] version >= 2.1.0
  * [CD-HIT][cdhit] version >= 4.6
  * [MUMmer][mummer] version >= 3.23

ARIBA also depends on several Python packages, all of which are available
via pip. Installing ARIBA with pip3 will get these automatically if they
are not already installed:
  * dendropy >= 4.2.0
  * matplotlib>=3.1.0
  * pyfastaq >= 3.12.0
  * pysam >= 0.9.1
  * pymummer >= 0.10.1
  * biopython

### Using pip3
Install ARIBA using pip:

    pip3 install ariba

### From Source
Download the latest release from this github repository or clone it. Run the tests:

    python3 setup.py test

**Note for OS X:** The tests require gawk which will need to be installed separately, e.g. via Homebrew.

If the tests all pass, install:

    python3 setup.py install

Alternatively, install directly from github using:

    pip3 install git+https://github.com/sanger-pathogens/ariba.git #--user

### Docker
ARIBA can be run in a Docker container. First install Docker, then install ARIBA:

    docker pull sangerpathogens/ariba

To use ARIBA use a command like this (substituting in your directories), where your files are assumed to be stored in /home/ubuntu/data:

    docker run --rm -it -v /home/ubuntu/data:/data sangerpathogens/ariba ariba -h

When calling Ariba via Docker (as above) you'll also need to add **/data/** in front of all the passed in file or directory names (e.g. /data/my_output_folder).


### Debian (Ariba version may not be the latest)
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

## Temporary files

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

## Usage
    usage: ariba <command> <options>

    optional arguments:
      -h, --help      show this help message and exit

    Available commands:

	aln2meta      Converts multi-aln fasta and SNPs to metadata
	expandflag    Expands flag column of report file
	flag          Translate the meaning of a flag
	getref        Download reference data
	micplot       Make violin/dot plots using MIC data
	prepareref    Prepare reference data for input to "run"
	pubmlstget    Download species from PubMLST and make db
	pubmlstspecies
		      Get list of available species from PubMLST
	refquery      Get cluster or sequence info from prepareref output
	run           Run the local assembly pipeline
	summary       Summarise multiple reports made by "run"
	test          Run small built-in test dataset
	version       Get versions and exit

Please read the [ARIBA wiki page][ARIBA wiki] for full usage instructions.

## License
ARIBA is free software, licensed under [GPLv3](https://github.com/sanger-pathogens/ariba/blob/master/LICENSE).

## Feedback/Issues
We currently do not have the resources to provide support for Ariba. However, the community might be able to help you out if you report any issues about usage of the software to the [issues page](https://github.com/sanger-pathogens/ariba/issues).

## Citation
If you use this software please cite:

ARIBA: rapid antimicrobial resistance genotyping directly from sequencing reads
Hunt M, Mather AE, Sánchez-Busó L, Page AJ, Parkhill J , Keane JA, Harris SR.
Microbial Genomics 2017. doi: [110.1099/mgen.0.000131](http://mgen.microbiologyresearch.org/content/journal/mgen/10.1099/mgen.0.000131)


  [ariba biorxiv]: http://biorxiv.org/content/early/2017/04/07/118000
  [bowtie2]: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
  [cdhit]: http://weizhongli-lab.org/cd-hit/
  [ARIBA wiki]: https://github.com/sanger-pathogens/ariba/wiki
  [mummer]: http://mummer.sourceforge.net/
  [python]: https://www.python.org/
