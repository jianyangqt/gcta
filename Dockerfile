# SPDX-License-Identifier: GPL-3.0
FROM ubuntu:24.04 as base

RUN apt -y update -qq &&  apt -y upgrade && \
	DEBIAN_FRONTEND=noninteractive apt -y install \
	--no-install-recommends --no-install-suggests \
		ca-certificates \
		libeigen3-dev \
		libgsl27 \
		libmkl-core \
		libmkl-intel-lp64 \
		libmkl-gnu-thread \
		libopenblas64-0-openmp \
		libopenblas64-0-pthread \
		libspectra-dev \
		libsqlite3-0 \
		libzstd1 \
		zlib1g-dev \
	&& \
	apt -y clean && rm -rf /var/lib/apt/lists/* /tmp/*

ENV MKLROOT=/usr
ENV OPENBLAS=/usr/include/x86_64-linux-gnu/openblas-openmp
ENV EIGEN3_INCLUDE_DIR=/usr/include/eigen3
ENV BOOST_LIB=/usr/include/boost/
ENV SPECTRA_LIB=/usr/include/Spectra

FROM base as builder

RUN apt -y update -qq && DEBIAN_FRONTEND=noninteractive apt -y install \
	--no-install-recommends --no-install-suggests \
		libboost1.74-dev \
		libgsl-dev \
		libmkl-dev \
		libopenblas-openmp-dev \
		libopenblas-pthread-dev \
		libsqlite3-dev \
		libz-dev \
		libzstd-dev \
		autoconf \
		automake \
		cmake \
		g++ \
		gcc \
		git \
		make \
		pkg-config \
        && \
        apt -y clean && rm -rf /var/lib/apt/lists/* /tmp/*

WORKDIR /usr/src

COPY . gcta
RUN cd gcta && git submodule update --init && mkdir build && \
	cd build && cmake .. && make

FROM base as runtime
ARG BASE
ARG RUN_CMD
ARG BUILD_REPO
ARG BUILD_TIME

COPY --chmod=0555 --from=builder /usr/src/gcta/build/ /usr/local/bin/

ENTRYPOINT [ "gcta" ]
LABEL org.opencontainers.image.description="GCTA (Genome-wide Complex Trait Analysis)"
LABEL org.opencontainers.image.authors="thebigcorporation@gmail.com"
LABEL org.opencontainers.image.source="https://github.com/jianyangqt/gcta.git"
