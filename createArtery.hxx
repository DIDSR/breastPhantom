/*! \file createArtery.hxx
 *  \brief breastPhantom artery creation header file
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

#ifndef CREATEARTERY_HXX_
#define CREATEARTERY_HXX_

#ifndef __VTKSMARTPOINTER__
	#define __VTKSMARTPOINTER__
	#include <vtkSmartPointer.h>
#endif

#ifndef __VTKIMAGEDATA__
	#define __VTKIMAGEDATA__
	#include <vtkImageData.h>
#endif

#ifndef __VTKIMAGEDATAWRITER__
	#define __VTKIMAGEDATAWRITER__
	#include <vtkXMLImageDataWriter.h>
#endif

#ifndef __VTKIMAGEDATAREADER__
	#define __VTKIMAGEDATAREADER__
	#include <vtkXMLImageDataReader.h>
#endif

#ifndef __TISSUESTRUCT__
	#define __TISSUESTRUCT__
	#include "tissueStruct.hxx"
#endif

#ifndef __ARTERY__
	#define __ARTERY__
	#include "artery.hxx"
#endif

void generate_artery(vtkImageData* breast, boost::program_options::variables_map vm, int* boundBox,
		     tissueStruct* tissue, double* sposPtr, double* sdirPtr, double* nipplePos, int seed, int mainSeed, bool firstTree);


#endif /* CREATEARTERY_HXX_ */
