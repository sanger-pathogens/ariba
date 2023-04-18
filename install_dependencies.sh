#!/usr/bin/env bash
set -vexu

install_root=$1

BOWTIE2_VERSION=2.3.1
SPADES_VERSION=3.13.1

apt-get update -qq
apt-get install -y software-properties-common
apt-add-repository universe
apt-get update -qq

apt-get install --no-install-recommends -y \
  build-essential \
  cd-hit \
  curl \
  git \
  libcurl4-gnutls-dev \
  libssl-dev \
  libbz2-dev \
  liblzma-dev \
  mummer \
  python3-dev \
  python3-setuptools \
  python3-pip \
  python3-pysam \
  python3-tk \
  python3-matplotlib \
  unzip \
  wget \
  zlib1g-dev

ln -s -f /usr/bin/python3 /usr/local/bin/python

if [ ! -d $install_root ]; then
  mkdir $install_root
fi
cd $install_root


# --------------- bowtie2 ------------------
cd $install_root
wget -q http://downloads.sourceforge.net/project/bowtie-bio/bowtie2/${BOWTIE2_VERSION}/bowtie2-${BOWTIE2_VERSION}-legacy-linux-x86_64.zip
unzip -n bowtie2-${BOWTIE2_VERSION}-legacy-linux-x86_64.zip
rm bowtie2-${BOWTIE2_VERSION}-legacy-linux-x86_64.zip
mv bowtie2-${BOWTIE2_VERSION}-legacy bowtie2


# --------------- spades -------------------
cd $install_root
wget -q https://github.com/ablab/spades/releases/download/v${SPADES_VERSION}/SPAdes-${SPADES_VERSION}-Linux.tar.gz
tar xf SPAdes-${SPADES_VERSION}-Linux.tar.gz
rm SPAdes-${SPADES_VERSION}-Linux.tar.gz
mv SPAdes-${SPADES_VERSION}-Linux SPAdes
