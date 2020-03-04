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
 perl \
 python3 \
 python3-pip \
 python3-dev \
 dos2unix \
 && pip3 install --upgrade pip

 WORKDIR /opt
 RUN mkdir remap-gff3
 COPY . /opt/remap-gff3/
 WORKDIR /opt/remap-gff3
 RUN find . -type f -print0 | xargs -0 dos2unix && pip3 install .
