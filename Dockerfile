FROM ubuntu:latest
RUN rm /bin/sh && ln -s /bin/bash /bin/sh

RUN apt-get update && apt-get install --yes \
 build-essential \
 gcc-multilib \
 apt-utils \
 zlib1g-dev \
 liblzo2-dev \
 libbz2-dev \
 liblzma-dev \
 libz-dev \
 libxml2-dev \
 libncurses5-dev \
 autotools-dev \
 git \
 perl \
 python \
 python-pip \
 python-dev \
 dos2unix

 WORKDIR /opt
 RUN mkdir remap-gff3
 COPY . /opt/remap-gff3/
 WORKDIR /opt/remap-gff3
 RUN find . -type f -print0 | xargs -0 dos2unix && pip install .
