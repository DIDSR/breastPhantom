/*! \file createDuct.hxx
 *  \brief breastPhantom duct creation header file
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

#ifndef CREATEDUCT_HXX_
#define CREATEDUCT_HXX_

#include <vtkVersion.h>

#ifndef __VTKIMAGEDATA__
    #define __VTKIMAGEDATA__
    #include <vtkImageData.h>
#endif

#ifndef __VTKPOINTS__
    #define __VTKPOINTS__
    #include <vtkPoints.h>
#endif

#ifndef __VTKDOUBLEARRAY__
    #define __VTKDOUBLEARRAY__
    #include <vtkDoubleArray.h>
#endif

#ifndef __TISSUESTRUCT__
    #define __TISSUESTRUCT__
    #include "tissueStruct.hxx"
#endif

#ifndef __DUCT__
    #define __DUCT__
    #include "duct.hxx"
#endif

void generate_duct(vtkImageData* breast, boost::program_options::variables_map vm, vtkPoints* TDLUloc, vtkDoubleArray* TDLUattr,
    unsigned char compartmentId, int* boundBox, tissueStruct* tissue, double* sposPtr, double* sdirPtr, int seed);

#endif /* CREATEDUCT_HXX_ */
