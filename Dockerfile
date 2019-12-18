FROM ubuntu:18.04

ENV DEBIAN_FRONTEND=noninteractive

MAINTAINER ariba-help@sanger.ac.uk

# Software version numbers
ARG BOWTIE2_VERSION=2.2.9
ARG SPADES_VERSION=3.13.1
ARG ARIBA_TAG=master
ARG ARIBA_BUILD_DIR=/ariba

ARG LOCALE_COUNTRY=en_GB

RUN apt-get -qq update && \
    apt-get install --no-install-recommends -y \
  build-essential \
  cd-hit \
  curl \
  git \
  libbz2-dev \
  liblzma-dev \
  mummer \
  python3-dev \
  python3-setuptools \
  python3-pip \
  python3-tk \
  python3-matplotlib \
  unzip \
  wget \
  zlib1g-dev

# Install locales
RUN apt-get update && apt-get install -y locales-all && rm -rf /var/lib/apt/lists/* 
# Set a default locale.
ENV LANG=${LOCALE_COUNTRY}.UTF-8 \
    LANGUAGE=${LOCALE_COUNTRY}:en 

# Install bowtie
RUN wget -q http://downloads.sourceforge.net/project/bowtie-bio/bowtie2/${BOWTIE2_VERSION}/bowtie2-${BOWTIE2_VERSION}-linux-x86_64.zip \
  && unzip bowtie2-${BOWTIE2_VERSION}-linux-x86_64.zip \
  && rm -f bowtie2-${BOWTIE2_VERSION}-linux-x86_64.zip

# Install SPAdes
RUN wget -q https://github.com/ablab/spades/releases/download/v${SPADES_VERSION}/SPAdes-${SPADES_VERSION}-Linux.tar.gz \
  && tar -zxf SPAdes-${SPADES_VERSION}-Linux.tar.gz \
  && rm -f SPAdes-${SPADES_VERSION}-Linux.tar.gz

# Need MPLBACKEND="agg" to make matplotlib work without X11, otherwise get the error
# _tkinter.TclError: no display name and no $DISPLAY environment variable
ENV ARIBA_BOWTIE2=$PWD/bowtie2-${BOWTIE2_VERSION}/bowtie2 ARIBA_CDHIT=cdhit-est MPLBACKEND="agg"
ENV PATH=$PATH:$PWD/SPAdes-${SPADES_VERSION}-Linux/bin

RUN ln -s -f /usr/bin/python3 /usr/local/bin/python

# Install Ariba
RUN mkdir -p $ARIBA_BUILD_DIR
COPY . $ARIBA_BUILD_DIR
RUN cd $ARIBA_BUILD_DIR \
  && python3 setup.py clean --all \
  && python3 setup.py test \
  && python3 setup.py install \
  && rm -rf $ARIBA_BUILD_DIR

CMD ariba
