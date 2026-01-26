FROM debian:trixie

ENV DEBIAN_FRONTEND=noninteractive
ENV LANG=C.UTF-8
ENV LC_ALL=C.UTF-8
ENV TERM=xterm

RUN  echo "deb http://deb.debian.org/debian bookworm main contrib non-free non-free-firmware" > /etc/apt/sources.list.d/my-sources.list && \
  echo "deb-src http://deb.debian.org/debian bookworm main contrib non-free non-free-firmware" >> /etc/apt/sources.list.d/my-sources.list && \
  apt-get update -y && apt-get upgrade -y \
  && apt-get install -y --no-install-recommends git samtools bcftools \
  build-essential wget ca-certificates zlib1g zlib1g-dev \
  libssl-dev libffi-dev libbz2-dev && \
  apt build-dep -y python3 python3-pip

WORKDIR /
COPY --chmod=755 wakhan /usr/local/bin
RUN mkdir build && cd build && \
  wget https://www.python.org/ftp/python/3.8.20/Python-3.8.20.tgz && \
  tar -xf Python-3.8.20.tgz && cd Python-3.8.20 && \
  ./configure --with-pip && make -j 10 && make install && \
  python3.8 -m pip install --upgrade pip && \
  pip3.8 install setuptools && cd / && rm -fr /build 
RUN cd /opt/ && git clone https://github.com/crankycrank/Wakhan.git && \
  cd Wakhan && git checkout bc48900 && pip3.8 install -r requirements.txt
#RUN pip install . --break-system-packages

# default interactive shell
CMD ["/bin/bash"]

# Add metadata labels
LABEL org.opencontainers.image.source="https://github.com/KolmogorovLab/Wakhan" \
      org.opencontainers.image.description="Wakhan: Copy number aberration analysis tool"
