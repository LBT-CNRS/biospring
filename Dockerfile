FROM ubuntu:latest AS base

RUN apt-get update && apt-get install -y \
    build-essential \
    autoconf \
    libtool \
    pkg-config \
    git \
    cmake \
    libnetcdf-c++4-dev \
    python3 \ 
    python3-pip \
    gdb

ARG CMAKE_BUILD_TYPE=Release
ENV CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
RUN echo "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}"

# Martini tools
FROM python:3 AS cg-tools
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir vermouth

# Copy FreeSASA
FROM base AS build_freesasa
COPY freesasa_src /freesasa/freesasa_src
WORKDIR /freesasa/freesasa_src
RUN  autoreconf -i && ./configure --prefix /freesasa --disable-json --disable-xml && \
     make -j1 && make install


FROM base AS build_mddriver
COPY mddriver_src /mddriver/mddriver_src
WORKDIR /mddriver/mddriver_src
# Ensure that we are in the "customdata" branch in the mddriver_src repo
RUN mkdir build ; cd build && \
    cmake .. -DCMAKE_INSTALL_PREFIX=/mddriver && \
    make -j1 && make install


FROM base AS build_biospring
# Copy freesasa lib and header
COPY --from=build_freesasa /freesasa/include/freesasa.h /usr/local/include
COPY --from=build_freesasa /freesasa/lib/libfreesasa.a /usr/local/lib
# Copy mddriver lib, headers and cmake file and set MDDriver_DIR
COPY --from=build_mddriver /mddriver/share/cmake/MDDriverConfig.cmake /mddriver/share/cmake/MDDriverConfig.cmake
COPY --from=build_mddriver /mddriver/include/* /mddriver/include/
COPY --from=build_mddriver /mddriver/lib/libmddriver.a /mddriver/lib/
ARG MDDriver_DIR=/mddriver/share/cmake/
RUN echo "The ARG variable value is $MDDriver_DIR"

# Remind to remove build dir in biospring_src and CMakeCache.txt
COPY ./ /biospring/biospring_src
WORKDIR /biospring/biospring_src
# Ensure that we are in the "develop" branch in the biospring_src repo
RUN mkdir build ; cd build && \
    cmake .. -DCMAKE_INSTALL_PREFIX:PATH=/biospring \
            -DBUILD_TESTS=OFF \
            -DMDDRIVER_SUPPORT=ON \
            -DOPENMP_SUPPORT=ON \
            -DFREESASA_SUPPORT=ON \
            -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} && \
    # Try to set -j1 if compilation error, restart deamon process
    make -j1 && make install

COPY --from=cg-tools /usr/local /usr/local

# Set environment variables as per BioSpring-config.sh
ENV BioSpring_ROOT="/biospring" \
    BioSpring_PATH="/biospring/bin" \
    BioSpring_DATA_PATH="/biospring/share/biospring/data" \
    BioSpring_SCRIPTS_PATH="/biospring/share/biospring/scripts" \
    BioSpring_DIR="/biospring/share/biospring/cmake" \
    BioSpring_LIBRARY_DIR="/biospring/lib" \
    BioSpring_INCLUDE_DIR="/biospring/include" \
    PATH="/biospring/bin:/usr/local/bin/:${PATH}" 

# Create the /data directory for external mounts
RUN mkdir -p /data

# Set /data as the working directory
WORKDIR /data
