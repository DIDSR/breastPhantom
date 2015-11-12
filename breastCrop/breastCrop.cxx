/*
 * breastCrop.cxx
 *
 *      Author: cgg
 *
 *		Input:  vti file unsigned int
 *		Output: raw binary unsigned int
 */
 
#include "breastCrop.hxx"

namespace po = boost::program_options;

int main(int argc, char* argv[]){

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
		("infile,i", po::value<std::string>()->required(), "input vti file")
	;

	po::variables_map vm;

	// get command line arguments
	po::store(parse_command_line(argc,argv,opts), vm);
	std::string inputFile = vm["infile"].as<std::string>();

	// read input phantom file
	vtkSmartPointer<vtkXMLImageDataReader> reader = 
		vtkSmartPointer<vtkXMLImageDataReader>::New();
	reader->SetFileName(inputFile.c_str());
	reader->Update();
	
	if(!reader->CanReadFile(inputFile.c_str())){
		cerr << "Cannot read input VTI file\n";
		return(1);
	}
	
	// read into memory
	vtkSmartPointer<vtkImageData> input =
		vtkSmartPointer<vtkImageData>::New();
		
	input->DeepCopy(reader->GetOutput());
	
	int dim2[3];
	input->GetDimensions(dim2);

	// find crop extent
	int extent[6];
	input->GetExtent(extent);
	int cropExtent[6];
	
	bool foundTissue = false;
	int idx = extent[0];
	
	while(!foundTissue){
		for(int j=extent[2]; j<=extent[3]; j++){
			for(int k=extent[4]; k<=extent[5]; k++){
				unsigned char* q = static_cast<unsigned char *>(input->GetScalarPointer(idx,j,k));
				if(*q != tissue.bg && *q != tissue.paddle){
					foundTissue = true;
				}
			}
		}
		if(!foundTissue){
			idx += 1;
		}
	}
	
	cropExtent[0] = idx;
	
	foundTissue = false;
	idx = extent[1];
	
	while(!foundTissue){
		for(int j=extent[2]; j<=extent[3]; j++){
			for(int k=extent[4]; k<=extent[5]; k++){
				unsigned char* q = static_cast<unsigned char *>(input->GetScalarPointer(idx,j,k));
				if(*q != tissue.bg && *q != tissue.paddle){
					foundTissue = true;
				}
			}
		}
		if(!foundTissue){
			idx -= 1;
		}
	}
	
	cropExtent[1] = idx;
	
	foundTissue = false;
	idx = extent[2];
	
	while(!foundTissue){
		for(int i=extent[0]; i<=extent[1]; i++){
			for(int k=extent[4]; k<=extent[5]; k++){
				unsigned char* q = static_cast<unsigned char *>(input->GetScalarPointer(i,idx,k));
				if(*q != tissue.bg && *q != tissue.paddle){
					foundTissue = true;
				}
			}
		}
		if(!foundTissue){
			idx += 1;
		}
	}
	
	cropExtent[2] = idx;
	
	foundTissue = false;
	idx = extent[3];
	
	while(!foundTissue){
		for(int i=extent[0]; i<=extent[1]; i++){
			for(int k=extent[4]; k<=extent[5]; k++){
				unsigned char* q = static_cast<unsigned char *>(input->GetScalarPointer(i,idx,k));
				if(*q != tissue.bg && *q != tissue.paddle){
					foundTissue = true;
				}
			}
		}
		if(!foundTissue){
			idx -= 1;
		}
	}
	
	cropExtent[3] = idx;
	
	foundTissue = false;
	idx = extent[4];
	
	while(!foundTissue){
		for(int i=extent[0]; i<=extent[1]; i++){
			for(int j=extent[2]; j<=extent[3]; j++){
				unsigned char* q = static_cast<unsigned char *>(input->GetScalarPointer(i,j,idx));
				if(*q != tissue.bg && *q != tissue.paddle){
					foundTissue = true;
				}
			}
		}
		if(!foundTissue){
			idx += 1;
		}
	}
	
	cropExtent[4] = idx;
	
	foundTissue = false;
	idx = extent[5];
	
	while(!foundTissue){
		for(int i=extent[0]; i<=extent[1]; i++){
			for(int j=extent[2]; j<=extent[3]; j++){
				unsigned char* q = static_cast<unsigned char *>(input->GetScalarPointer(i,j,idx));
				if(*q != tissue.bg && *q != tissue.paddle){
					foundTissue = true;
				}
			}
		}
		if(!foundTissue){
			idx -= 1;
		}
	}
	
	cropExtent[5] = idx;
	
	// crop phantom
	
	input->Crop(cropExtent);
	
	// output file name
	char outfilename[256];
	char infilename[256];
	
	strcpy(infilename, inputFile.c_str());
	
	char *infilenameHead;
	
	infilenameHead = strtok(infilename, ".");
	
	sprintf(outfilename, "%s_cropped.vti", infilenameHead);
	
	// output cropped vti file
	vtkSmartPointer<vtkXMLImageDataWriter> writer = 
		vtkSmartPointer<vtkXMLImageDataWriter>::New();
		
	writer->SetFileName(outfilename);
	writer->SetInputData(input);
	writer->Write();
	
	return(0);
}

									
	
	
	
