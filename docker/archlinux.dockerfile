# syntax=docker/dockerfile:1
FROM archlinux:base-devel

RUN <<EOT
    pacman -Suy --noconfirm

    # build tools
    pacman -S --noconfirm gcc git wget cmake make
    # dependencies
    pacman -S --noconfirm vtk fmt openmp openmpi lapack boost
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
