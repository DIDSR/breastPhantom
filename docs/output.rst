Output Formats
==============

Voxel Values
------------

The main output of the software is a 3D voxelized breast phantom where each voxel is an unsigned 8 bit integer tissue label (possible values 0-255).
Here is a lookup table specifying the labels for each tissue type.

============================  =====
Tissue                        Label
============================  =====
air                           0
fat                           1
skin                          2
glandular                     29
nipple                        33
muscle                        40
ligament                      88
TDLU\ :sup:`1`                95
duct                          125
artery                        150
vein                          225
cancerous mass\ :sup:`2`      200
calcification\ :sup:`2`       250
compression paddle\ :sup:`2`  50
============================  =====

:sup:`1`\ Terminal Duct Lobular Unit

:sup:`2`\ These tissues are not created by the phantom generation software, but can be added by other VICTRE pipeline software.  They are included here for completeness.

Output Files
------------

For each run of the breast phantom generation software several files are created in the output directory specified in the configuration file.
The naming convention is that all output filenames begin with ``p_`` followed by the random number generator seed that was used and an extention specifying the type of output file.
Here is a list of output files:

p\_\ *nnnnnnnn*\ .cfg
    A copy of the configuration file used to generate the phantom

p\_\ *nnnnnnnn*\ .vti
    The phantom volume stored in VTK image format (XML with zlib-compressed data element).  This is a file format native to the Visualization Toolkit.  See https://www.vtk.org for more information.
    This file can be opened in Paraview, an open-source, multi-platform data analysis and visualization application available for download at https://www.paraview.org/ 

p\_\ *nnnnnnnn*\ .raw.gz
    The raw phantom volume stored as 8-bit unsigned integers in a gzip archive.  There is no file header.

p\_\ *nnnnnnnn*\ .mhd
    A metaImage header file containing information about the phantom data stored in the file p\_\ *nnnnnnnn*\ .raw.gz  Parsing this header file will allow you to read and manipulate raw.gz phantom files.
    See https://itk.org/Wiki/MetaIO/Documentation for more information on the MetaImage format
    
p\_\ *nnnnnnnn*\ .loc
    This file lists candidate locations for insertion of masses and calcifications.  Each location is where a TDLU has been randomly generated.  Terminal duct lobular units are a common site for cancer formation.
    Locations are stored in plain text, one location per line with comma separated coordinates in millimeters.  See the phantom geometry section for an explanation of the coordinate system.
