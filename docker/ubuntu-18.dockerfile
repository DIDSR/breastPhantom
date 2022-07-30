# syntax=docker/dockerfile:1
FROM ubuntu:18.04

# Disable tzdata questions
ARG DEBIAN_FRONTEND="noninteractive"
ENV TZ="Etc/UTC"

# Enable extra repositories for GCC and CMake
RUN <<EOT
    apt-get update
    apt-get install -y wget gpg software-properties-common lsb-release

    # Ubuntu Test channel for GCC
    add-apt-repository ppa:ubuntu-toolchain-r/test -y

    # PPA for newer Boost versions
    add-apt-repository ppa:mhier/libboost-latest

    # CMake official binaries
    wget -O - https://apt.kitware.com/keys/kitware-archive-latest.asc | gpg --dearmor - \
        > /usr/share/keyrings/kitware-archive-keyring.gpg
    echo "deb [signed-by=/usr/share/keyrings/kitware-archive-keyring.gpg] https://apt.kitware.com/ubuntu/ $(lsb_release -cs) main" \
        > /etc/apt/sources.list.d/kitware.list

    # update packages
    apt-get update
    apt-get upgrade -y

    # install kitware keyring to avoid broken keys
    rm /usr/share/keyrings/kitware-archive-keyring.gpg
    apt-get install -y kitware-archive-keyring
EOT

# Install dependencies
RUN <<EOT
    # Boost
    apt-get install -y libboost1.68-dev
    # C++20 compiler
    apt-get install -y build-essential gcc-11 g++-11 cmake
    # VTK
    apt-get install -y vtk7 libvtk7-dev
    # Vectorization libs
    apt-get install -y libomp-dev libblas-dev liblapack-dev
    # other deps
    apt-get install -y git wget
EOT

# Copy source files
COPY src /breastPhantom/src
COPY cfg /breastPhantom/cfg
COPY CMakeLists.txt /breastPhantom

# Build binary
WORKDIR /breastPhantom
# RUN <<EOT
#     CC=gcc-11 CXX=g++-11 cmake -B build .
#     make -C build
# EOT

# # Run the binary
# ENTRYPOINT ["build/breastPhantom"]
