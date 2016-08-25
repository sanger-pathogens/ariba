#!/bin/bash
set -e
set -x

start_dir=$(pwd)

BOWTIE2_VERSION=2.2.8
CDHIT_VERSION=4.6.5
MASH_VERSION=1.1
SAMTOOLS_VERSION=1.3
MUMMER_VERSION=3.23

BOWTIE2_DOWNLOAD_URL="http://downloads.sourceforge.net/project/bowtie-bio/bowtie2/${BOWTIE2_VERSION}/bowtie2-${BOWTIE2_VERSION}-linux-x86_64.zip"
CDHIT_DOWNLOAD_URL="https://github.com/weizhongli/cdhit/archive/V${CDHIT_VERSION}.tar.gz"
MASH_DOWNLOAD_URL="https://github.com/marbl/Mash/releases/download/v${MASH_VERSION}/mash-Linux64-v${MASH_VERSION}.tar.gz"
SAMTOOLS_DOWNLOAD_URL="https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2"
MUMMER_DOWNLOAD_URL="http://downloads.sourceforge.net/project/mummer/mummer/${MUMMER_VERSION}/MUMmer${MUMMER_VERSION}.tar.gz"


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


# --------------- mash --------------------
cd $build_dir
download $MASH_DOWNLOAD_URL "mash-Linux64-v${MASH_VERSION}.tar.gz"
mash_dir="$build_dir/mash-Linux64-v${MASH_VERSION}"
tar zxf mash-Linux64-v${MASH_VERSION}.tar.gz


# --------------- mummer ------------------
cd $build_dir
download $MUMMER_DOWNLOAD_URL "MUMmer${MUMMER_VERSION}.tar.gz"
mummer_dir="$build_dir/MUMmer${MUMMER_VERSION}"
tar -zxf MUMmer${MUMMER_VERSION}.tar.gz
cd $mummer_dir
make


# --------------- samtools -----------------
cd $build_dir
download $SAMTOOLS_DOWNLOAD_URL "samtools-${SAMTOOLS_VERSION}.tar.bz2"
samtools_dir="$build_dir/samtools-${SAMTOOLS_VERSION}"
tar -xjf samtools-${SAMTOOLS_VERSION}.tar.bz2
cd $samtools_dir
make



cd $start_dir

update_path () {
  new_dir=$1
  if [[ ! "$PATH" =~ (^|:)"${new_dir}"(:|$) ]]; then
    export PATH=${new_dir}:${PATH}
  fi
}

update_path ${bowtie2_dir}
update_path ${cdhit_dir}
update_path ${mash_dir}
update_path ${mummer_dir}
update_path ${samtools_dir}
