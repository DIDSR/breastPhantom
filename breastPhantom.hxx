/*! \file breastPhantom.hxx
 *  \brief breastPhantom main header file
 *  \author Christian G. Graff
 *  \version 1.0
 *  \date 2018
 *
 *  \copyright To the extent possible under law, the author(s) have
 *  dedicated all copyright and related and neighboring rights to this
 *  software to the public domain worldwide. This software is
 *  distributed without any warranty.  You should have received a copy
 *  of the CC0 Public Domain Dedication along with this software.
 *  If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
 *
 */

#ifndef __BREASTPHANTOM_HXX__
#define __BREASTPHANTOM_HXX__

#include <iostream>
#include <cstdint>

#include <math.h>
#include <unistd.h>
#include <sys/stat.h>
// debug
#include <time.h>
#include <zlib.h>
#include <omp.h>

#include "perlinNoise.hxx"
#include "createDuct.hxx"
#include "createArtery.hxx"
#include "createVein.hxx"

// vtk stuff
#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkMath.h>
#include <vtkPolyData.h>
#include <vtkCleanPolyData.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkPolyDataNormals.h>
#include <vtkGeometryFilter.h>
#include <vtkDecimatePro.h>
#include <vtkTIFFWriter.h>
#include <vtkFloatArray.h>
#include <vtkCellArray.h>
#include <vtkPointData.h>
#include <vtkSortDataArray.h>
#include <vtkParametricSpline.h>
#include <vtkProperty.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkSmoothPolyDataFilter.h>
#include <vtkSurfaceReconstructionFilter.h>
#include <vtkContourFilter.h>
#include <vtkReverseSense.h>
#include <vtkImageData.h>
#include <vtkFillHolesFilter.h>
#include <vtkCellLocator.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkXMLImageDataReader.h>
#include <vtkMetaImageWriter.h>
#include <vtkPointLocator.h>
#include <vtkOctreePointLocator.h>
#include <vtkMinimalStandardRandomSequence.h>
#include <vtkVector.h>

// number of fat lobule Fourier perturbation coefficients
#define NUMCOEFF 3

#include "tissueStruct.hxx"

#endif /* __BREASTPHANTOM_HXX__ */
