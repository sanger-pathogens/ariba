FROM ubuntu:18.04

ENV DEBIAN_FRONTEND=noninteractive

MAINTAINER ariba-help@sanger.ac.uk

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

RUN wget -q http://downloads.sourceforge.net/project/bowtie-bio/bowtie2/2.2.9/bowtie2-2.2.9-linux-x86_64.zip \
  && unzip bowtie2-2.2.9-linux-x86_64.zip \
  && rm -f bowtie2-2.2.9-linux-x86_64.zip

RUN wget -q https://github.com/ablab/spades/releases/download/v3.13.1/SPAdes-3.13.1-Linux.tar.gz \
  && tar -zxf SPAdes-3.13.1-Linux.tar.gz \
  && rm -f SPAdes-3.13.1-Linux.tar.gz

# Need MPLBACKEND="agg" to make matplotlib work without X11, otherwise get the error
# _tkinter.TclError: no display name and no $DISPLAY environment variable
ENV ARIBA_BOWTIE2=$PWD/bowtie2-2.2.9/bowtie2 ARIBA_CDHIT=cdhit-est MPLBACKEND="agg"
ENV PATH=$PATH:$PWD/SPAdes-3.13.1-Linux/bin

RUN cd /usr/local/bin && ln -s /usr/bin/python3 python && cd

RUN git clone https://github.com/sanger-pathogens/ariba.git \
  && cd ariba \
  && git checkout v2.14.2 \
  && python3 setup.py test \
  && python3 setup.py install

CMD ariba
