/*
 * breastPhantom.hxx
 *
 *  Created on: Dec 22, 2014
 *      Author: cgg
 */

#ifndef BREASTPHANTOM_HXX_
#define BREASTPHANTOM_HXX_

#ifndef __IOS__
	#define __IOS__
	#include <iostream>
#endif

#ifndef __CSDTINT__
	#define __CSTDINT__
	#include <cstdint>
#endif

#include <math.h>
#include <unistd.h>
#include <sys/stat.h>
// debug
#include <time.h>

#ifndef __OMP__
	#define __OMP__
	#include <omp.h>
#endif

#include "perlinNoise.hxx"
#include "createDuct.hxx"
#include "createArtery.hxx"
#include "createVein.hxx"

// vtk stuff

#include <vtk/vtkVersion.h>
#ifndef __VTKSMARTPOINTER__
	#define __VTKSMARTPOINTER__
	#include <vtk/vtkSmartPointer.h>
#endif
#ifndef __VTKMATH__
	#define __VTKMATH__
	#include <vtk/vtkMath.h>
#endif
#include <vtk/vtkPolyData.h>
#include <vtk/vtkCleanPolyData.h>
#include <vtk/vtkVertexGlyphFilter.h>
//#include <vtk/vtkDelaunay3D.h>
#include <vtk/vtkPolyDataNormals.h>
#include <vtk/vtkGeometryFilter.h>
#include <vtk/vtkDecimatePro.h>
#include <vtk/vtkTIFFWriter.h>
#include <vtk/vtkFloatArray.h>
#include <vtk/vtkCellArray.h>
#include <vtk/vtkPointData.h>
#include <vtk/vtkParametricSpline.h>
#include <vtk/vtkProperty.h>
#include <vtk/vtkTransform.h>
#include <vtk/vtkTransformPolyDataFilter.h>
#include <vtk/vtkSmoothPolyDataFilter.h>
#include <vtk/vtkSurfaceReconstructionFilter.h>
#include <vtk/vtkContourFilter.h>
#include <vtk/vtkReverseSense.h>
#ifndef __VTKIMAGEDATA__
	#define __VTKIMAGEDATA__
	#include <vtk/vtkImageData.h>
#endif
#include <vtk/vtkFillHolesFilter.h>
#include <vtk/vtkCellLocator.h>
#include <vtk/vtkXMLPolyDataWriter.h>
#ifndef __VTKIMAGEDATAWRITER__
	#define __VTKIMAGEDATAWRITER__
	#include <vtk/vtkXMLImageDataWriter.h>
#endif

#ifndef __VTKIMAGEDATAREADER__
	#define __VTKIMAGEDATAREADER__
	#include <vtk/vtkXMLImageDataReader.h>
#endif
#include <vtk/vtkPointLocator.h>
#include <vtk/vtkOctreePointLocator.h>
#include <vtk/vtkMinimalStandardRandomSequence.h>
#ifndef __VTKVECTOR__
	#define __VTKVECTOR__
	#include <vtk/vtkVector.h>
#endif

// number of fat lobule Fourier perturbation coefficients
#define NUMCOEFF 3

#ifndef __TISSUESTRUCT__
	#define __TISSUESTRUCT__
	#include "tissueStruct.hxx"
#endif

// status bar
inline void statusBar(unsigned int current, unsigned int total, unsigned int width = 40, unsigned int numUpdate = 50); 


#endif /* BREASTPHANTOM_HXX_ */
