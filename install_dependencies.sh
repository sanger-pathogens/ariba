#!/bin/bash
set -e
set -x

start_dir=$(pwd)

BOWTIE2_VERSION=2.3.1
CDHIT_VERSION=4.6.5
MUMMER_VERSION=3.23

BOWTIE2_DOWNLOAD_URL="http://downloads.sourceforge.net/project/bowtie-bio/bowtie2/${BOWTIE2_VERSION}/bowtie2-${BOWTIE2_VERSION}-legacy-linux-x86_64.zip"
CDHIT_DOWNLOAD_URL="https://github.com/weizhongli/cdhit/archive/V${CDHIT_VERSION}.tar.gz"
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
download $BOWTIE2_DOWNLOAD_URL "bowtie2-${BOWTIE2_VERSION}-legacy.zip"
bowtie2_dir="$build_dir/bowtie2-${BOWTIE2_VERSION}-legacy"
unzip -n bowtie2-${BOWTIE2_VERSION}-legacy.zip


# --------------- cdhit --------------------
cd $build_dir
download $CDHIT_DOWNLOAD_URL "cdhit-${CDHIT_VERSION}.tar.gz"
tar -zxf cdhit-${CDHIT_VERSION}.tar.gz
cdhit_dir="$build_dir/cdhit-${CDHIT_VERSION}"
cd $cdhit_dir
make


# --------------- mummer ------------------
cd $build_dir
download $MUMMER_DOWNLOAD_URL "MUMmer${MUMMER_VERSION}.tar.gz"
mummer_dir="$build_dir/MUMmer${MUMMER_VERSION}"
tar -zxf MUMmer${MUMMER_VERSION}.tar.gz
cd $mummer_dir
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
update_path ${mummer_dir}
