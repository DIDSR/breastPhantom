Installation
============

The VICTRE Breast Phantom source code is available `here <https://github.com/DIDSR/breastPhantom>`_.  It is written in C++ and therefore must be compiled for your system before use.

Operating System
----------------

This software has been developed and tested on Debian and Red Hat based Linux distributions.  Compilation is via cmake and should be possible on any UNIX-like operationg system (including MacOS) with minimal
modification to the build configuration in the provided CMakeLists.txt file, however the following distributions are recommended: 

- Ubuntu 16.04
- Centos 7.x

The following build instructions are based on Ubuntu 16.04.

Compilation
-----------

First, install the packages necessary to compile the software (reminder: instructions are specific to Ubuntu,
consult your operating system documentation for analagous package names and installation instructions)  From a terminal command prompt::

    > sudo apt-get install cmake vtk6 libvtk6-dev \
			   libblas-dev liblapack-dev \
			   libboost-dev libboost-program-options-dev \
			   libproj-dev zlib1g-dev

Next, clone the git repository::

    > git clone https://github.com/DIDSR/breastPhantom

From a suitable build directory generate a Makefile using cmake and run make::

    > cmake <path to repository>
    > make

where <path to repository> is the directory containing the cloned repository.

.. note:: when running cmake you may recieve an error message related to references to files "/usr/lib/x86_64-linux-gnu/libvtkRenderingPythonTkWidgets.so" and "/usr/bin/vtk".
	  This is a known VTK bug but since we are not linking to these files, the error can safely be ignored.

Installation
------------

If the previous steps were successful, the executable breastPhantom should be present in the build directory.  This file can be copied to a convenient location.
No additional installation steps are required.  Consult the next section which describes how to run the executable.
