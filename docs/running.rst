Running the Software
====================

The executable breastPhantom is run from the command line and takes a single required command line argument::

    > breastPhantom -c [FILE]

where ``[FILE]`` is a configuration file which specifies all the breast model parameters.

Configuration file
------------------

Several sample configuration files are provided, all ending with the extension ``.cfg``. Configuration files can be modified with your favorite text editor. A configuration file has many parameters, please see [#f1]_ for a complete list with brief explanations.

A few of the parameters that most users will probably want to modify are the following:

*outputDir*
    The directory where the generated phantom files will be written.

*imgRes*
    The phantom voxel size (in mm).  Typically this varies between 0.05 and 0.25 depending on application.

*targetFatFrac*
    The desired breast density in terms of volume fractional fat content.  Should be between 0 and 1, though for typical breasts should vary between approximately 0.4 and 0.95.

.. [#f1] All configuration parameters are explained in :ref:`config params` and can also be listed with ``breastPhantom --help``.

Generating a phantom
--------------------

Once the desired configuration file is created, the software is run from the command line as above.  Several status messages will be generated during the phantom creation process. The time required to generate a phantom can vary significantly and is highly dependent on the voxel size (imgRes parameter) and the processing capabilities of your computer. A low resolution phantom on a fast computer may require 30 mins of processing time, while a high-resolution phantom on a slow computer may take more than a day.  Sections of this software are multi-threaded and should take advantage of multiple cores if available. The number of cores used can be controlled with the environment variable ``OMP_NUM_THREADS`` if desired. All cores will be used by default.

All phantom output files will be placed in the directory specified in the configuration file

Randomization
-------------

The phantom generation process is stochastic and relies upon a random number generator seed.  This seed can be specified in the configuration file using the parameter *seed*, though the default behavior is to randomly generate a seed at runtime by reading from ``/dev/urandom``. In this mode each time the breast phantom generation code is run with a fixed configuration file, a different random breast will be generated that is consistent with the configuration parameters. This allows for the generation of random breast textures at fixed glandularity, for example.
