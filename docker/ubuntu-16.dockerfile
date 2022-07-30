# syntax=docker/dockerfile:1
FROM ubuntu:16.04

# Disable tzdata questions
ARG DEBIAN_FRONTEND="noninteractive"
ENV TZ="Etc/UTC"

# Install dependencies
RUN <<EOT
    apt-get update
    apt-get upgrade -y
    # C++20 compiler
    apt-get install -y build-essential gcc g++ cmake
    # VTK
    apt-get install -y vtk6 libvtk6-dev
    # Vectorization libs
    apt-get install -y libomp-dev libblas-dev liblapack-dev
    # Boost
    apt-get install -y libboost-dev libboost-program-options-dev
    # other deps
    apt-get install -y git wget zlib1g-dev
EOT

# Copy source files
COPY src /breastPhantom/src
COPY cfg /breastPhantom/cfg
COPY CMakeLists.txt /breastPhantom

# Build binary
WORKDIR /breastPhantom
RUN <<EOT
    mkdir build && cd build
    cmake ..
    make
EOT

# Run the binary
ENTRYPOINT ["build/breastPhantom"]
