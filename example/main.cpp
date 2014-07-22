#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <stdint.h>
#include <vector>
#include <string.h>
#include <utility>
#include <algorithm>
#include <omp.h>

#include <time.h>
#include <sys/time.h>

#include "SimpleCT.h"

using namespace std;

int main( int argc, const char** argv ) {

	if (argc != 6) {
		std::cout << "Usage:" << std::endl;
		std::cout << "./simpleExample test_128x128.uint8 0.1 32 ct ct_img.raw" << std::endl;
		return false;
	}

	struct timespec t1, t2;
	clock_gettime(CLOCK_MONOTONIC,  &t1);

	unsigned short dim = 2;

	char fileName[512];
	strcpy(fileName,argv[1]);
	
	char *xSize, *ySize, *zSize, *type;
	
	std::cout << std::endl;
	std::cout << "Input file: \t \t \t " << fileName << std::endl;
	std::cout << "Image dimensions: \t " ;

	type  	= strstr(fileName,".");
	*type 	= '\0';
	type++;

	xSize 	= strstr(fileName,"_");
	*xSize 	= '\0';
	xSize++;

	ySize 	= strstr(xSize,"x");
	*ySize 	= '\0';
	ySize++;

	zSize 	= strstr(ySize,"x");
	if (zSize != NULL) {
		*zSize = '\0';
		zSize++;
		dim = 3;
	}

	std::vector<unsigned int> size;
	size.push_back(atoi(xSize));
	size.push_back(atoi(ySize));

	std::cout << "\t x = " << size[0];
	std::cout << "\t y = " << size[1];

	int numberOfElements 	= size[0]*size[1];

	if (dim == 3) {
		size.push_back(atoi(zSize));
		numberOfElements 	= numberOfElements*size[2];
		std::cout << "\t z = " << size[2];
	}

	std::cout << std::endl;

	std::cout << "File type: \t \t \t " << type << std::endl;

	double* img;
	img 	= new double[numberOfElements];

	FILE * inputImageFile;
	int numberOfReadBytes;
	inputImageFile 	= fopen (argv[1],"r");

	if ( strcmp(type,"uint8") == 0 ) {
		unsigned char* uint8_t_type;
		uint8_t_type 			= new unsigned char[numberOfElements];
		numberOfReadBytes  		= fread (uint8_t_type, 1, numberOfElements, inputImageFile);
		if (numberOfReadBytes != numberOfElements) {
			fputs ("Reading error",stderr);
			exit (1);
		}
		for (int i = 0; i < numberOfElements; ++i)
			img[i] = double(uint8_t_type[i]);
		delete[] uint8_t_type;
	} else if ( strcmp(type,"uint16") == 0 ) {
		unsigned short* uint16_t_type;
		uint16_t_type 			= new unsigned short[numberOfElements];
		numberOfReadBytes  		= fread (uint16_t_type, 2, numberOfElements, inputImageFile);
		if (numberOfReadBytes != numberOfElements) {
			fputs ("Reading error",stderr);
			exit (1);
		}
		for (int i = 0; i < numberOfElements; ++i)
			img[i] = double(uint16_t_type[i]);
		delete[] uint16_t_type;
	} else {
		std::cerr << "Unknown file type : " << type << std::endl;
		return false;
	}
	fclose (inputImageFile);

	ContourTree ct(img,size);
	ct.prune(uint(atoi(argv[3])), double(atof(argv[2])));

	std::cout << std::endl;
	ct.print_Short();

	char ctFileName[512];
	strcpy(ctFileName,argv[4]);

	char ctImgFileName[512];
	strcpy(ctImgFileName,argv[5]);

	ct.print_CT(ctFileName);
	ct.print_CT_img(ctImgFileName);

	cout<<"Intensity threshold: \t\t "	<< double(atof(argv[2])) << endl;
	cout<<"Area/volume threshold: \t\t "	<< uint(atoi(argv[3])) << endl << endl;

	delete[] img;

	clock_gettime(CLOCK_MONOTONIC,  &t2);

	cout<<"Total computation time:\t\t "<< (t2.tv_sec - t1.tv_sec) + (double) (t2.tv_nsec - t1.tv_nsec) * 1e-9 << endl;
	std::cout << std::endl;

	return true;

}
