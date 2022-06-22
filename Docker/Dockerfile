FROM ubuntu:18.04


MAINTAINER bhaas@broadinstitute.org


RUN apt-get update && apt-get install -y gcc g++ perl python3 automake make \
                                       wget curl libdb-dev \
                                       bzip2 zlibc zlib1g zlib1g-dev  default-jre \
                       python3-setuptools python3-dev build-essential python3-distutils \
                                       unzip libbz2-dev  liblzma-dev git && \
    apt-get clean


# make python3 be the default python
RUN ln -sf /usr/bin/python3 /usr/bin/python

RUN curl https://bootstrap.pypa.io/pip/3.6/get-pip.py -o get-pip.py && \
  python get-pip.py



RUN pip install igv-reports==1.0.1 requests pysam pandas

RUN curl -L https://cpanmin.us | perl - App::cpanminus
RUN cpanm install URI::Escape


ENV SRC /usr/local/src
ENV BIN  /usr/local/bin

## Samtools
ENV SAMTOOLS_VERSION=1.7

RUN SAMTOOLS_URL="https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2" && \
   cd $SRC && \
   wget $SAMTOOLS_URL && \
   tar xvf samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
   cd samtools-${SAMTOOLS_VERSION}/htslib-${SAMTOOLS_VERSION} && ./configure && make && make install && \
   cd ../ && ./configure --without-curses && make && make install
               


COPY intron_counter.py $BIN/

WORKDIR $SRC

ENV CTAT_SPLICING_VERSION=0.0.2
ENV CO=c3e18ca1f0d61f6e9a93ab2c68731d93ac1c0b6d
RUN git clone https://github.com/NCIP/CTAT-SPLICING.git && \
  cd CTAT-SPLICING && \
  git checkout $CO && \
  git submodule init && git submodule update && \
  git submodule foreach --recursive git submodule init && \
  git submodule foreach --recursive git submodule update && \
  make


