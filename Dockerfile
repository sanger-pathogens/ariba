FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

ARG ARIBA_BUILD_DIR=/ariba
ARG DEPS_DIR=/bioinf-tools
ARG LOCALE_COUNTRY=en_GB

# Install locales
RUN apt-get update && apt-get install -y locales-all && rm -rf /var/lib/apt/lists/*
# Set a default locale.
ENV LANG=${LOCALE_COUNTRY}.UTF-8 \
    LANGUAGE=${LOCALE_COUNTRY}:en

RUN mkdir -p $ARIBA_BUILD_DIR
COPY . $ARIBA_BUILD_DIR
RUN $ARIBA_BUILD_DIR/install_dependencies.sh $DEPS_DIR

# Need MPLBACKEND="agg" to make matplotlib work without X11, otherwise get the error
# _tkinter.TclError: no display name and no $DISPLAY environment variable
ENV ARIBA_BOWTIE2=$DEPS_DIR/bowtie2/bowtie2 \
    ARIBA_CDHIT=cdhit-est \
    MPLBACKEND="agg" \
    PATH=$PATH:$DEPS_DIR/SPAdes/bin

# Install Ariba
RUN cd $ARIBA_BUILD_DIR \
  && python3 -m pip install -r requirements.txt \
  && python3 setup.py clean --all \
  && python3 setup.py test \
  && python3 setup.py install \
  && rm -rf $ARIBA_BUILD_DIR

CMD ariba
