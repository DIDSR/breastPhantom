/*
 * breastMass.hxx
 *
 *  Created on: Oct 19, 2015
 *      Author: Christian Graff
 */

#ifndef BREASTMASS_HXX_
#define BREASTMASS_HXX_

// boost
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <boost/program_options.hpp>

#include <iostream>
#include <math.h>
#include <unistd.h>
#include <sys/stat.h>

#include <omp.h>

// vtk
#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkMath.h>
#include <vtkPNGWriter.h>
#include <vtkImageData.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkMinimalStandardRandomSequence.h>
#include <vtkBoxMuellerRandomSequence.h>
#include <vtkImageInterpolator.h>

// create spiculation segments
void createBranch(double, double, double, double, double, double, double, double, 
	double, vtkImageData*, vtkMinimalStandardRandomSequence*, vtkBoxMuellerRandomSequence*, boost::program_options::variables_map);


#endif /* BREASTMASS_HXX_ */
