FROM debian:trixie

RUN DEBIAN_FRONTEND=noninteractive apt-get update -y && apt-get upgrade -y && apt-get install -y --no-install-recommends python3 python3-setuptools python3-pysam libhts3 samtools gawk sed minimap2 mafft exonerate git ca-certificates bwa python3-pip build-essential zlib1g-dev python3-dev

ENV LANG=C.UTF-8                                                                              
ENV LC_ALL=C.UTF-8                                                                            
ENV TERM=xterm                                                                                

WORKDIR /                                                                                         
RUN git clone https://github.com/adamewing/tldr.git && cd tldr \
  git checkout 7c6dfda &&  python3 setup.py install

CMD ["/bin/bash"]
