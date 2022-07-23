/*! \file createVein.hxx
 *  \brief breastPhantom vein creation header file
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

#ifndef __CREATEVEIN_HXX__
#define __CREATEVEIN_HXX__

#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkImageData.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkXMLImageDataReader.h>

#include "tissueStruct.hxx"
#include "vein.hxx"

void generate_vein(vtkImageData* breast, boost::program_options::variables_map vm, int* boundBox,
           tissueStruct* tissue, double* sposPtr, double* sdirPtr, double* nipplePos, int seed, int mainSeed, bool firstTree);


#endif /* __CREATEVEIN_HXX__ */
