# syntax=docker/dockerfile:1
FROM ubuntu:22.04

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
    apt-get install -y vtk9 libvtk9-dev
    # Vectorization libs
    apt-get install -y libomp-dev libblas-dev liblapack-dev
    # Boost
    apt-get install -y libboost-dev libboost-program-options-dev libboost-iostreams-dev
    # other deps
    apt-get install -y git wget
EOT

# Copy source files
COPY src /breastPhantom/src
COPY cfg /breastPhantom/cfg
COPY CMakeLists.txt /breastPhantom

# Build binary
WORKDIR /breastPhantom
RUN <<EOT
    cmake -B build .
    make -C build
EOT

# Run the binary
ENTRYPOINT ["build/breastPhantom"]
