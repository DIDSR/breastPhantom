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

#ifndef __CREATEDUCT_HXX__
#define __CREATEDUCT_HXX__

#include <vtkVersion.h>
#include <vtkImageData.h>
#include <vtkPoints.h>
#include <vtkDoubleArray.h>

#include "tissueStruct.hxx"
#include "duct.hxx"

void generate_duct(vtkImageData* breast, boost::program_options::variables_map vm, vtkPoints* TDLUloc, vtkDoubleArray* TDLUattr,
    unsigned char compartmentId, int* boundBox, tissueStruct* tissue, double* sposPtr, double* sdirPtr, int seed);

#endif /* __CREATEDUCT_HXX__ */
