FROM ubuntu:latest
RUN rm /bin/sh && ln -s /bin/bash /bin/sh

RUN apt-get update && apt-get install --yes \
 build-essential \
 gcc-multilib \
 apt-utils \
 zlib1g-dev \
 libncurses5-dev \
 autotools-dev \
 git \
 perl \
 python \
 libbz2-dev \
 liblzma-dev \
 libz-dev \
 libxml2-dev \
 python-pip \
 python-dev

 WORKDIR /opt
 RUN git clone https://github.com/NAL-i5K/remap-gff3.git
 WORKDIR /opt/remap-gff3
 RUN wget https://pypi.python.org/packages/55/db/fa76af59a03c88ad80494fc0df2948740bbd58cd3b3ed5c31319624687cc/bx-python-0.7.3.tar.gz
 RUN pip install --upgrade pip && pip install \
  numpy \
  bx-python-0.7.3.tar.gz \
  CrossMap \
  gff3tool \