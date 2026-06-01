# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Dockerfile for annotating VCFs with VEP

FROM ubuntu:focal-20220531
MAINTAINER "Ryan Collins <Ryan_Collins@dfci.harvard.edu>"

# Install essential system tools and Perl environment
ENV DEBIAN_FRONTEND=noninteractive
ENV DEBIAN_PRIORITY=critical
RUN apt-get clean && \
    apt-get -qqy update && \
    apt-get -qqy install --fix-missing \
      apt-utils \
      build-essential \
      zlib1g-dev \
      libbz2-dev \
      liblzma-dev \
      libperl-dev \
      libcurl4-openssl-dev \
      libssl-dev \
      libdbi-perl \
      libdbd-mysql-perl \
      default-libmysqlclient-dev \
      bioperl \
      curl \
      wget \
      file \
      git \
      perl \
      gzip \
      less \
      libmodule-build-perl \
      libarchive-zip-perl \
      libarchive-extract-perl \
      libtry-tiny-perl \
      libwww-perl \
      libscalar-util-numeric-perl \
      libset-intervaltree-perl \
      libdevel-cycle-perl \
      libpadwalker-perl \
      libtest-warn-perl \
      libipc-run-perl \
      cpanminus \
      libpng-dev \
      libfreetype6-dev \
      uuid-dev

# Install additional Perl modules needed for VEP + plugins
RUN cpanm --quiet --notest \
    Bio::Root::Version \
    Test::Warnings \
    PadWalker \
    Devel::Cycle \
    Set::IntervalTree \
    Scalar::Util::Numeric \
    Bio::EnsEMBL::Utils::IO

# Make VEP cache and plugin directories
RUN mkdir /opt/vep_cache /opt/vep_plugins
ENV VEP_CACHE=/opt/vep_cache
ENV VEP_PLUGINS=/opt/vep_plugins

# Install htslib
ARG HTSLIB_VERSION=1.18
RUN cd /opt/ && \
    wget https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 && \
    tar xjf htslib-${HTSLIB_VERSION}.tar.bz2 && \
    rm htslib-${HTSLIB_VERSION}.tar.bz2 && \
    cd htslib-${HTSLIB_VERSION} && \
    ./configure --prefix=/opt/htslib-1.18 --enable-shared && \
    make && \
    make install && \
    cd -
ENV CFLAGS="-I/opt/htslib-1.18/include"
ENV LDFLAGS="-L/opt/htslib-1.18/lib"
ENV LD_LIBRARY_PATH="/opt/htslib-1.18/lib:$LD_LIBRARY_PATH"

# Clone Ensembl repos
RUN cd /opt && \
    git clone https://github.com/Ensembl/ensembl.git && \
    git clone https://github.com/Ensembl/ensembl-test.git && \
    git clone https://github.com/Ensembl/Bio-DB-HTS.git && \
    git clone https://github.com/Ensembl/ensembl-vep.git && \
    cd -
ENV PERL5LIB=/opt/ensembl/modules:/opt/ensembl-test/modules:/opt/Bio-DB-HTS/lib

# Install Bio::DB::HTS
ENV HTSLIB_DIR=/opt/htslib-1.18
RUN cd /opt/Bio-DB-HTS && \
    perl Build.PL && \
    ./Build && \
    rm -f t/10leak.t t/11housekeeping_apache2.t && \
    ./Build test && \
    ./Build install && \
    perl -e "use Bio::DB::HTS" && \
    cd -

# Install libBigWig and UCSC Kent source
RUN cd /opt && \
    git clone https://github.com/dpryan79/libBigWig.git && \
    cd libBigWig && \
    make && \
    make install && \
    cd -
RUN cd /opt && \
    wget http://hgdownload.cse.ucsc.edu/admin/jksrc.zip && \
    unzip jksrc.zip && \
    rm jksrc.zip && \
    cd kent/src/lib && \
    CXXFLAGS=-fPIC CFLAGS=-fPIC CPPFLAGS=-fPIC make && \
    cd - 

# Ensure the Kent library is located where expected by Bio::DB::BigWig
RUN set -eux; \
    MACH="$(uname -m)"; \
    echo "machine arch: $MACH"; \
    mkdir -p /opt/kent/src/lib/"$MACH"; \
    JK="$(find /opt/kent/src/lib -type f -name 'jkweb.a' | head -n1)"; \
    if [ -z "$JK" ]; then echo "ERROR: jkweb.a not found under /opt/kent/src/lib"; exit 1; fi; \
    cp -n "$JK" /opt/kent/src/lib/"$MACH"/
ENV KENT_SRC=/opt/kent/src


# Install Bio::DB::Big with cpanm
RUN cpanm --quiet --notest Bio::DB::Big

# Install Bio::DB::BigWig via manual CPAN tarball build (non-interactive)
RUN set -eux; \
    MACH="$(uname -m)"; \
    HTSLIB_DIR=/opt/htslib-1.18; \
    cd /opt && \
    wget http://www.cpan.org/authors/id/L/LD/LDS/Bio-BigFile-1.07.tar.gz && \
    tar xzf Bio-BigFile-1.07.tar.gz && \
    cd Bio-BigFile-1.07 && \
    perl Build.PL --extra_linker_flags="/opt/kent/src/lib/${MACH}/jkweb.a -L${HTSLIB_DIR} -lhts -lz -lssl -lcrypto" && \
    ./Build && \
    ./Build install && \
    cd - && \
    rm -rf /opt/Bio-BigFile-1.07 /opt/Bio-BigFile-1.07.tar.gz

# Install VEP
ARG VEP_VERSION=115
RUN cd /opt/ensembl-vep && \
    perl INSTALL.pl \
      --NO_HTSLIB \
      --NO_TEST \
      --CACHE_VERSION $VEP_VERSION \
      --AUTO ap \
      --SPECIES homo_sapiens \
      --PLUGINS UTRAnnotator,dbNSFP,SpliceAI,AlphaMissense,PrimateAI \
      --CONVERT \
      --CACHEDIR $VEP_CACHE/ \
      --PLUGINSDIR $VEP_PLUGINS/ && \
    cd -

# Install LOFTEE (not available directly through VEP install)
RUN cd /opt && \
    git clone -b grch38 https://github.com/konradjk/loftee.git && \
    cp -rf loftee/* $VEP_PLUGINS/ && \
    cd -

# Alias VEP and HTSLib tools on path
ENV PATH=/opt/ensembl-vep/:/opt/htslib-1.18/bin:$PATH

# Launch bash as entrypoint
CMD ["/bin/bash"]
