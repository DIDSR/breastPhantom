Installation
============

The VICTRE Breast Phantom source code is available on `Github <https://github.com/DIDSR/breastPhantom>`_. It is written in C++ and therefore must be compiled for your system before use.

Operating System
----------------

This software has been developed and tested on Debian and Red Hat based Linux distributions. Compilation is via CMake and should be possible on any UNIX-like operationg system (including MacOS) with minimal modification to the build configuration in the provided CMakeLists.txt file, however the following distributions are recommended:

- Ubuntu 22.04
- Arch Linux

Compilation for Ubuntu 22.04
----------------------------

First, install the packages necessary to compile the software::

    sudo apt-get install build-essential git gcc g++ cmake vtk9 libvtk9-dev \
        libblas-dev liblapack-dev libboost-dev libboost-program-options-dev libboost-iostreams-dev

Next, clone the git repository::

    git clone https://github.com/DIDSR/breastPhantom

.. code-block::

    cd breastPhantom

Then build the binary in a new ``build`` directory::

    cmake . -B build

.. code-block::

    make -C build

If the previous steps were successful, the executable ``breastPhantom`` should be present in the ``build`` directory. This file can be copied to a convenient location. No additional installation steps are required. Consult the next section :ref:`running` which describes how to run the executable.

Compilation for Ubuntu 20.04
----------------------------

First, install the packages necessary to compile the software::

    sudo apt-get install build-essential gcc-10 g++-10 cmake vtk7 libvtk7-dev \
        libblas-dev liblapack-dev libboost-dev libboost-program-options-dev libboost-iostreams-dev

Next, clone the git repository::

    git clone https://github.com/DIDSR/breastPhantom

.. code-block::

    cd breastPhantom

Then build the binary in a new ``build`` directory with GCC 10::

    mkdir build && cd build

.. code-block::

    CC=gcc-10 CXX=g++-10 cmake ..

.. code-block::

    make

If the previous steps were successful, the executable ``breastPhantom`` should be present in the ``build`` directory. This file can be copied to a convenient location. No additional installation steps are required. Consult the next section :ref:`running` which describes how to run the executable.

.. note:: When running CMake you may recieve an error message related to references to files ``/usr/lib/x86_64-linux-gnu/libvtkRenderingPythonTkWidgets.so`` and ``/usr/bin/vtk``. This is a known VTK bug but since we are not linking to these files, the error can safely be ignored.

Compilation for Arch Linux
--------------------------

First, install the packages necessary to compile the software::

    sudo pacman -Suy gcc git wget cmake make vtk fmt openmp openmpi lapack boost

Next, clone the git repository::

    git clone https://github.com/DIDSR/breastPhantom

.. code-block::

    cd breastPhantom

Then build the binary in a new ``build`` directory::

    cmake . -B build

.. code-block::

    make -C build

If the previous steps were successful, the executable ``breastPhantom`` should be present in the ``build`` directory. This file can be copied to a convenient location. No additional installation steps are required. Consult the next section which describes how to run the executable.
