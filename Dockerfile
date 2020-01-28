FROM ubuntu

# metadata
LABEL base.image="ubuntu"
LABEL version="1"
LABEL software="nextPARSParer"
LABEL software.version="0.10"
LABEL description=""
LABEL website="https://github.com/Gabaldonlab/nextPARSParer"
LABEL license="GNU General Public License"

MAINTAINER Uciel Chorostecki "uciel.chorostecki@bsc.es"

ARG SOURCE_DIR=/root/src/nextparsPARSER/dependencies
RUN mkdir -p SOURCE_DIR

# ---------------------------------------

RUN echo "Installing python and relevant tools..."
RUN apt-get update && apt-get install -y \ 
    build-essential \
    libxml2-dev \
    libxslt-dev \
    libssl-dev \
    libyaml-dev \
    libffi-dev \
    python \
	python3 \
    python-dev \
    python-pip \
    vim \
	zlib1g-dev 

# General dev tools
RUN echo "Installing general dev tools..."

RUN apt-get install -y git

RUN echo "Latest versions of python tools via pip..."
RUN pip install --upgrade pip \
                          virtualenv \
                          requests

# WORKDIR /root/src/
# RUN git clone https://github.com/Gabaldonlab/nextPARSParser.git

# ---------------------------------------

WORKDIR $SOURCE_DIR

RUN echo "Installing STAR"
RUN git clone https://github.com/alexdobin/STAR.git

RUN echo "Compiling STAR" 
WORKDIR $SOURCE_DIR/STAR/source
RUN make STAR

RUN PATH=$PATH:/root/src/nextparsPARSER/dependencies/STAR/source
RUN export PATH



