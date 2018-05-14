FROM ubuntu:latest
RUN rm /bin/sh && ln -s /bin/bash /bin/sh

RUN apt-get update && apt-get install --yes \
 build-essential \
 zlib1g-dev \
 liblzo2-dev \
 git \
 perl \
 python \
 python-pip \
 python-dev

 WORKDIR /opt
 RUN git clone https://github.com/NAL-i5K/remap-gff3.git
 WORKDIR /opt/remap-gff3
 RUN pip install --upgrade pip && pip install \
  gff3 \
  gff3tool \
  CrossMap