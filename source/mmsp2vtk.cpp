// File:    mmsp2vtk.cpp
// Purpose: reads MMSP grid containing sparse floats, converts to LegacyVTK
// Output:  VTK file
// Depends: MMSP, zlib

// Questions/Comments to kellet@rpi.edu (Trevor Keller)

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <zlib.h>

#include "MMSP.hpp"

int main(int argc, char* argv[]) {
	if ( argc != 3 ) {
		std::cout << "Usage: " << argv[0] << " data.dat output.vtk\n";
		return ( 1 );
	}

	// file open error check
	std::ifstream input(argv[1]);
	if (!input) {
		std::cerr << "File input error: could not open " << argv[1] << ".\n\n";
		exit(-1);
	}

	// read data type
	std::string type;
	getline(input, type, '\n');

	// grid type error check: read line, "grid:sparse:float"
	if (type.substr(0, 4) != "grid") {
		std::cerr << "File input error: file does not contain grid data." << std::endl;
		exit(-1);
	} else if (type.substr(5, 6) != "sparse") {
		std::cerr << "File input error: grid does not contain sparse data." << std::endl;
		exit(-1);
	} else if (type.substr(12, 5) != "float") {
		std::cerr << "File input error: vector data does not contain floats." << std::endl;
		exit(-1);
	}

	// read grid dimension
	int dim;
	input >> dim;

	std::vector<MMSP::vector<int> > points;
	std::vector<float> weights;

	if (dim == 3) {
		// construct grid object
		MMSP::grid<3, MMSP::sparse<float> > grid(argv[1]);
		std::ofstream output(argv[2]);
		if (!output) {
			std::cerr << "File output error: could not create " << argv[2] << ".\n\n";
			exit(-1);
		}
		// Write VTK header
		output<<"# vtk DataFile Version 2.0\n"
					<<"From file "<<argv[1]<<'\n'
		      <<"ASCII\n"
		      <<"DATASET STRUCTURED_POINTS\n"
		      <<"DIMENSIONS";
		for (int d=0; d<dim; ++d)
		  output<<' '<<MMSP::x1(grid,d)-MMSP::x0(grid,d);
		output<<"\nORIGIN 0 0 0\n"
		      <<"SPACING";
		for (int d=0; d<dim; ++d)
		  output<<' '<<MMSP::dx(grid,d);
		output<<"\nPOINT_DATA "<<MMSP::nodes(grid)
		      <<"\nSCALARS grainid int\n"
		      <<"LOOKUP_TABLE default\n";
		// Write grid data
		for (int n=0; n<nodes(grid); ++n)
		  output<<grid(n).grain_id()<<'\n';
		output.close();
	} else {
		std::cerr << "Error: " << dim << "-D data is not supported!" << std::endl;
		exit(1);
	}

	return 0;
}
