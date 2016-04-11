#!/bin/bash
set -e
set -x

start_dir=$(pwd)

BCFTOOLS_VERSION=1.3
BOWTIE2_VERSION=2.2.8
CDHIT_VERSION=4.6.5
SAMTOOLS_VERSION=1.3
MUMMER_VERSION=3.23
SPADES_VERSION=3.6.0

BCFTOOLS_DOWNLOAD_URL="https://github.com/samtools/bcftools/releases/download/1.3/bcftools-${BCFTOOLS_VERSION}.tar.bz2"
BOWTIE2_DOWNLOAD_URL="http://downloads.sourceforge.net/project/bowtie-bio/bowtie2/${BOWTIE2_VERSION}/bowtie2-${BOWTIE2_VERSION}-linux-x86_64.zip"
CDHIT_DOWNLOAD_URL="https://github.com/weizhongli/cdhit/archive/V${CDHIT_VERSION}.tar.gz"
SAMTOOLS_DOWNLOAD_URL="https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2"
MUMMER_DOWNLOAD_URL="http://downloads.sourceforge.net/project/mummer/mummer/${MUMMER_VERSION}/MUMmer${MUMMER_VERSION}.tar.gz"
SPADES_DOWNLOAD_URL="http://spades.bioinf.spbau.ru/release${SPADES_VERSION}/SPAdes-${SPADES_VERSION}-Linux.tar.gz"


# Make an install location
if [ ! -d 'build' ]; then
  mkdir build
fi
cd build
build_dir=$(pwd)

# DOWNLOAD ALL THE THINGS
download () {
  url=$1
  download_location=$2

  if [ -e $download_location ]; then
    echo "Skipping download of $url, $download_location already exists"
  else
    echo "Downloading $url to $download_location"
    wget $url -O $download_location
  fi
}


# --------------- bcftools -----------------
cd $build_dir
download $BCFTOOLS_DOWNLOAD_URL "bcftools-${BCFTOOLS_VERSION}.tar.bz2"
bcftools_dir="$build_dir/bcftools-${BCFTOOLS_VERSION}"
tar -xjf bcftools-${BCFTOOLS_VERSION}.tar.bz2
cd $bcftools_dir
make


# --------------- bowtie2 ------------------
cd $build_dir
download $BOWTIE2_DOWNLOAD_URL "bowtie2-${BOWTIE2_VERSION}.zip"
bowtie2_dir="$build_dir/bowtie2-${BOWTIE2_VERSION}"
unzip -n bowtie2-${BOWTIE2_VERSION}.zip


# --------------- cdhit --------------------
cd $build_dir
download $CDHIT_DOWNLOAD_URL "cdhit-${CDHIT_VERSION}.tar.gz"
tar -zxf cdhit-${CDHIT_VERSION}.tar.gz
cdhit_dir="$build_dir/cdhit-${CDHIT_VERSION}"
cd $cdhit_dir
make


# --------------- samtools -----------------
cd $build_dir
download $SAMTOOLS_DOWNLOAD_URL "samtools-${SAMTOOLS_VERSION}.tar.bz2"
samtools_dir="$build_dir/samtools-${SAMTOOLS_VERSION}"
tar -xjf samtools-${SAMTOOLS_VERSION}.tar.bz2
cd $samtools_dir
make


# --------------- mummer ------------------
cd $build_dir
download $MUMMER_DOWNLOAD_URL "MUMmer${MUMMER_VERSION}.tar.gz"
mummer_dir="$build_dir/MUMmer${MUMMER_VERSION}"
tar -zxf MUMmer${MUMMER_VERSION}.tar.gz
cd $mummer_dir
make


# --------------- spades -----------------
cd $build_dir
download $SPADES_DOWNLOAD_URL "SPAdes-${SPADES_VERSION}-Linux.tar.gz"
spades_dir="$build_dir/SPAdes-${SPADES_VERSION}-Linux/bin"
tar -zxf SPAdes-${SPADES_VERSION}-Linux.tar.gz


cd $start_dir

update_path () {
  new_dir=$1
  if [[ ! "$PATH" =~ (^|:)"${new_dir}"(:|$) ]]; then
    export PATH=${new_dir}:${PATH}
  fi
}

update_path ${bcftools_dir}
update_path ${bowtie2_dir}
update_path ${cdhit_dir}
update_path ${mummer_dir}
update_path ${samtools_dir}
update_path ${spades_dir}


# -------------- R packages ---------------
mkdir -p ~/R/libs
echo "R_LIBS=~/R/libs" > ~/.Renviron
wget https://cran.r-project.org/src/contrib/Archive/ape/ape_3.1.tar.gz
R CMD INSTALL ape_3.1.tar.gz

