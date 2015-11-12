/*
 * create_vein.hxx
 *
 *  Created on: Dec 22, 2014
 *      Author: cgg
 */

#ifndef CREATEVEIN_HXX_
#define CREATEVEIN_HXX_

#ifndef __VTKSMARTPOINTER__
	#define __VTKSMARTPOINTER__
	#include <vtk/vtkSmartPointer.h>
#endif

#ifndef __VTKIMAGEDATA__
	#define __VTKIMAGEDATA__
	#include <vtk/vtkImageData.h>
#endif

#ifndef __VTKIMAGEDATAWRITER__
	#define __VTKIMAGEDATAWRITER__
	#include <vtk/vtkXMLImageDataWriter.h>
#endif

#ifndef __VTKIMAGEDATAREADER__
	#define __VTKIMAGEDATAREADER__
	#include <vtk/vtkXMLImageDataReader.h>
#endif

#ifndef __TISSUESTRUCT__
	#define __TISSUESTRUCT__
	#include "tissueStruct.hxx"
#endif

#ifndef __VEIN__
	#define __VEIN__
	#include "vein.hxx"
#endif

void generate_vein(vtkImageData* breast, boost::program_options::variables_map vm, int* boundBox,
	tissueStruct* tissue, double* sposPtr, double* sdirPtr, int seed, bool firstTree);


#endif /* CREATEVEIN_HXX_ */
