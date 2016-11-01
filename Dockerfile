#
# This container will install ARIBA from master
#
FROM debian:testing

#
#  Authorship
#
MAINTAINER ap13@sanger.ac.uk

#
# Install the dependancies
#
RUN apt-get update -qq && apt-get install -y git bowtie2 cd-hit fastaq libc6 libfml0 libgcc1 libminimap0 libstdc++6 mummer python3 python3-setuptools python3-dev python3-pysam python3-pymummer python3-dendropy gcc g++ zlib1g-dev

#
# Get the latest code from github and install
#
RUN git clone https://github.com/sanger-pathogens/ariba.git && cd ariba && python3 setup.py install
