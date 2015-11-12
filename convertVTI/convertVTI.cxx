/*
 * convertVTI.cxx
 *
 *      Author: cgg
 *
 *		Input:  vti file unsigned int
 *		Output: raw binary unsigned int
 */
 
#include "convertVTI.hxx"

namespace po = boost::program_options;

int main(int argc, char* argv[]){

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
	int numEl = dim2[1]*dim2[2];
	
	// output file name
	char outfilename[256];
	char infilename[256];
	
	strcpy(infilename, inputFile.c_str());
	
	char *infilenameHead;
	
	infilenameHead = strtok(infilename, ".");
	
	sprintf(outfilename, "%s_%d_%d_%d.raw", infilenameHead, dim2[0], dim2[1], dim2[2]);
	
	// output raw format
	FILE *out;
	
	out = fopen(outfilename, "wb");
	
	unsigned char* p = static_cast<unsigned char*>(input->GetScalarPointer());
	
	for(int i=0; i<dim2[0]; i++){
		fwrite(&(p[i*numEl]), numEl, sizeof(unsigned char), out);
	}
	fclose(out);
	
	return(0);
}

									
	
	
	
