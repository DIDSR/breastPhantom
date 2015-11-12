/*
 * breastSphereMass.cxx
 *
 *      Author: cgg
 *
 *		Input:  voxelized breast model, mass x,y,z,diameter
 *		Output: voxelized breast model with inserted lesion
 */
 
#include "breastSphereMass.hxx"

namespace po = boost::program_options;

int main(int argc, char* argv[]){

	double pi = vtkMath::Pi();
	
	tissueStruct tissue;
        
	tissue.bg = 0;
	tissue.skin = 2;
	tissue.nipple = 33;
	tissue.fat = 1;
	tissue.cooper = 88;
	tissue.gland = 29;
	tissue.TDLU = 95;
	tissue.duct = 125;
	tissue.artery = 150;
	tissue.vein = 225 ;
	tissue.muscle = 40;
	tissue.mass = 200;
	tissue.calc = 250;
	tissue.paddle = 50;
	
	// command line parameters
	po::options_description opts("All options");
	opts.add_options()
		("infile,i", po::value<std::string>(), "input voxelized breast")
		("xloc,x", po::value<double>(), "mass x location (mm)")
		("yloc,y", po::value<double>(), "mass y location (mm)")
		("zloc,z", po::value<double>(), "mass z location (mm)")
		("diam,d", po::value<double>(), "mass diameter (mm)")
	;

	po::variables_map vm;

	// get command line arguments
	po::store(parse_command_line(argc,argv,opts), vm);
	std::string inputFile;
	if(vm.count("infile")){
		inputFile = vm["infile"].as<std::string>();
	} else {
		cerr << "No input breast VTI file specified\n";
		return(EXIT_FAILURE);
	}
	
	double centerPos[3];
	double diameter;
	
	if(vm.count("xloc") && vm.count("yloc") && vm.count("zloc") && vm.count("diam")){
		centerPos[0] = vm["xloc"].as<double>();
		centerPos[1] = vm["yloc"].as<double>();
		centerPos[2] = vm["zloc"].as<double>();
		diameter = vm["diam"].as<double>();
	} else {
		cerr << "Must specify mass x,y,z coordinates and diameter\n";
		return(EXIT_FAILURE);
	}
	
	// read input phantom file
	vtkSmartPointer<vtkXMLImageDataReader> reader = 
		vtkSmartPointer<vtkXMLImageDataReader>::New();
	reader->SetFileName(inputFile.c_str());
	reader->Update();
	
	if(!reader->CanReadFile(inputFile.c_str())){
		cerr << "Cannot read input breast VTI file\n";
		return(EXIT_FAILURE);
	}
	
	vtkSmartPointer<vtkImageData> breast =
		vtkSmartPointer<vtkImageData>::New();
		
	breast->DeepCopy(reader->GetOutput());
	
	// check center is within field of view
	int centerIdx[3];	
	double pcoords[3];
	int inVol = breast->ComputeStructuredCoordinates(centerPos, centerIdx, pcoords);
	
	if(!inVol){
		cerr << "Mass center outside phantom field of view\n";
		return(EXIT_FAILURE);
	}
	
	// check center is within breast
	unsigned char *p = static_cast<unsigned char *>(breast->GetScalarPointer(centerIdx));
								
	if(*p == tissue.bg || *p == tissue.paddle){
		cerr << "Mass center outside breast\n";
		return(EXIT_FAILURE);
	}
	
	int breastExtent[6];
	breast->GetExtent(breastExtent);
	
	double breastSpacing[3];
	breast->GetSpacing(breastSpacing);
	
	// always isotropic resolution
	double imgRes = breastSpacing[0];
	
	int searchD = static_cast<int>(ceil(diameter/2.0/imgRes+1.0));
	
	int searchSpace[6] = {centerIdx[0]-searchD, centerIdx[0]+searchD,
		centerIdx[1]-searchD, centerIdx[1]+searchD, centerIdx[2]-searchD, centerIdx[2]+searchD};
		
	searchSpace[0] = (searchSpace[0] < breastExtent[0]) ? breastExtent[0] : searchSpace[0];
	searchSpace[1] = (searchSpace[1] > breastExtent[1]) ? breastExtent[1] : searchSpace[1];
	searchSpace[2] = (searchSpace[2] < breastExtent[2]) ? breastExtent[2] : searchSpace[2];
	searchSpace[3] = (searchSpace[3] > breastExtent[3]) ? breastExtent[3] : searchSpace[3];
	searchSpace[4] = (searchSpace[4] < breastExtent[4]) ? breastExtent[4] : searchSpace[4];
	searchSpace[5] = (searchSpace[5] > breastExtent[5]) ? breastExtent[5] : searchSpace[5];
	
	bool warnOutside = false;
	int massVoxels = 0;
	
	for(int a=searchSpace[0]; a<=searchSpace[1]; a++){
		for(int b=searchSpace[2]; b<=searchSpace[3]; b++){
			for(int c=searchSpace[4]; c<=searchSpace[5]; c++){
				int idx[3] = {a,b,c};
				double loc[3];
				breast->GetPoint(breast->ComputePointId(idx), loc);
				if(vtkMath::Distance2BetweenPoints(loc, centerPos) <= diameter*diameter/4.0){
					// inside mass
					unsigned char *q = static_cast<unsigned char *>(breast->GetScalarPointer(idx));
					if(*q == tissue.bg || *q == tissue.paddle){
						warnOutside = true;
					} else {
						*q = tissue.mass;
						massVoxels += 1;
					}
				}
			}
		}
	}
	
	cout << "Analytic mass vol. = " << 4.0/3.0*pi*(diameter/2)*(diameter/2)*(diameter/2) << 
		" mm^3, Voxelized mass vol. = " << massVoxels*(imgRes*imgRes*imgRes) << " mm^3\n";
		
	// save with inserted mass
	
	vtkSmartPointer<vtkXMLImageDataWriter> writeMass =
		vtkSmartPointer<vtkXMLImageDataWriter>::New();

	writeMass->SetFileName("breastMass.vti");
	writeMass->SetInputData(breast);
	writeMass->Write();
	
	return(EXIT_SUCCESS);
}
	
						
	
	
	
	
