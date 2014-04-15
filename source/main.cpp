// main.cpp, modified from MMSP.main.hpp

// The user must supply the following in any source
// code that includes this file:
//
//     #include"..."
//
//     std::string PROGRAM = "...";
//     std::string MESSAGE = "...";
//     typedef ... GRID2D;
//     typedef ... GRID3D;
//
//     #include"MMSP.main.hpp"
//
//
// The first include must provide the functions
//
//     void generate(int dim,
//                   const char* filename);
//
//     void update(GRID2D& grid, int steps);
//
// which the main() function calls to generate
// example grids or to perform computations.


#ifndef MMSP_MAIN
#define MMSP_MAIN
#include<iostream>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<cstdlib>
#include<cctype>
#include"graingrowth.cpp"

int main(int argc, char* argv[]) {
	MMSP::Init(argc, argv);


	// check argument list
	if (argc < 2) {
		std::cout << PROGRAM << ": bad argument list.  Use\n\n";
		std::cout << "    " << PROGRAM << " --help\n\n";
		std::cout << "to generate help message.\n\n";
		exit(-1);
	}



	// print help message and exit
	if (std::string(argv[1]) == std::string("--help")) {
		std::cout << PROGRAM << ": " << MESSAGE << "\n\n";
		std::cout << "Valid command lines have the form:\n";
		std::cout << "    " << PROGRAM << " ";
		std::cout << "[--help] [--init dimension [outfile]] [infile [outfile] steps [increment]]\n\n";
		std::cout << "A few examples of using the command line follow.\n\n";
		std::cout << "The command\n";
		std::cout << "    " << PROGRAM << " --help\n";
		std::cout << "generates this help message and exits.  ";
		std::cout << "The \"--init\" option can be used to initialize the grid with a Voronoi tessellation, e.g.\n";
		std::cout << "    " << PROGRAM << " --init 2\n";
		std::cout << "generates the Voronoi tessellation on a grid of dimension 2 and writes it to the \n";
		std::cout << "file named \"voronoi.dat\".\n\n";
		std::cout << "    " << PROGRAM << " --init 3 voronoi\n";
		std::cout << " or " << PROGRAM << " --init 3 voronoi.dat\n";
		std::cout << "generates the Voronoi tessellation on a grid of dimension 3 and writes it to the \n";
		std::cout << "file named \"voronoi.dat\".\n\n";
		std::cout << "    " << PROGRAM << " polycrystal.dat 1000\n";
		std::cout << "reads the grid contained within \"polycrystal.dat\" and runs a simulation for 1000 time steps.\n";
		std::cout << "The final grid is written to a file named \"polycrystal.1000.dat\".\n\n";
		std::cout << "    " << PROGRAM << " polycrystal.dat 1000 100\n";
		std::cout << "reads the grid contained within \"polycrystal.dat\" and runs a simulation for 1000 time steps.\n";
		std::cout << "The grid is then written to a file every 100 time steps.\n";
		std::cout << "The resulting files are \nnamed \"polycrystal.0100.dat\", \"polycrystal.0200.dat\", ... \"polycrystal.1000.dat\".\n\n";
		std::cout << "    " << PROGRAM << " voronoi.dat polycrystal.dat 1000\n";
		std::cout << "reads the grid contained within \"voronoi.dat\" and runs a simulation for 1000 time steps.\n";
		std::cout << "The final grid is written to a file named \"polycrystal.1000.dat\".\n\n";
		std::cout << "    " << PROGRAM << " polycrystal.0000.dat 1000 100\n";
		std::cout << "reads the grid contained within \"polycrystal.0000.dat\" and runs a simulation for 1000 time steps.\n";
		std::cout << "The grid is then written to a file every 100 time steps.\n";
		std::cout << "The resulting files are \nnamed \"polycrystal.0100.dat\", \"polycrystal.0200.dat\", ... \"polycrystal.1000.dat\".\n\n";
		std::cout << "    " << PROGRAM << " polycrystal.1000.dat 2000 100\n";
		std::cout << "reads the grid contained within \"polycrystal.1000.dat\" and runs a simulation for 1000 additional time steps.\n";
		std::cout << "The grid is then written to a file every 100 time steps.\n";
		std::cout << "The resulting files are \nnamed \"polycrystal.1100.dat\", \"polycrystal.1200.dat\", ... \"polycrystal.2000.dat\".\n\n";
		exit(0);
	}



	// generate initial grid
	else if (std::string(argv[1]) == std::string("--init")) {
		// check argument list
		if (argc<3 or argc>4) {
			std::cout << PROGRAM << ": bad argument list.  Use\n\n";
			std::cout << "    " << PROGRAM << " --help\n\n";
			std::cout << "to generate help message.\n\n";
			exit(-1);
		}

		// check problem dimension
		if (std::string(argv[2]).find_first_not_of("0123456789") != std::string::npos) {
			std::cout << PROGRAM << ": initial grid must have integral dimension.  Use\n\n";
			std::cout << "    " << PROGRAM << " --help\n\n";
			std::cout << "to generate help message.\n\n";
			exit(-1);
		}

		int dim = atoi(argv[2]);

		// dimension must be 2 or 3
		if (dim<2 or dim>3) {
			std::cout<<PROGRAM<<": initial grid must be of dimension 2 or 3.  Use\n\n";
			std::cout<<"    "<<PROGRAM<<" --help\n\n";
			std::cout<<"to generate help message.\n\n";
			exit(-1);
		}

		// set output file name
		std::string outfile;
		if (argc < 4) outfile = "voronoi.dat";
		else outfile = argv[3];

		// tessellate
		char* filename = new char[outfile.length()];
		for (unsigned int i=0; i<outfile.length(); i++)
			filename[i] = outfile[i];
		MMSP::generate(dim, filename);
		delete [] filename;
	}



	// run simulation
	else {
		// bad argument list
		if (argc<3 or argc>5) {
			std::cout << PROGRAM << ": bad argument list.  Use\n\n";
			std::cout << "    " << PROGRAM << " --help\n\n";
			std::cout << "to generate help message.\n\n";
			exit(-1);
		}

		int steps;
		int increment;
		std::string outfile;

		if (std::string(argv[2]).find_first_not_of("0123456789") == std::string::npos) {
			// set output file name
			outfile = argv[1];

			// must have integral number of time steps
			if (std::string(argv[2]).find_first_not_of("0123456789") != std::string::npos) {
				std::cout << PROGRAM << ": number of time steps must have integral value.  Use\n\n";
				std::cout << "    " << PROGRAM << " --help\n\n";
				std::cout << "to generate help message.\n\n";
				exit(-1);
			}

			steps = atoi(argv[2]);
			increment = steps;

			if (argc > 3) {
				// must have integral output increment
				if (std::string(argv[3]).find_first_not_of("0123456789") != std::string::npos) {
					std::cout << PROGRAM << ": output increment must have integral value.  Use\n\n";
					std::cout << "    " << PROGRAM << " --help\n\n";
					std::cout << "to generate help message.\n\n";
					exit(-1);
				}

				increment = atoi(argv[3]);

				// output increment must be smaller than number of steps
				if (increment > steps) {
					std::cout << PROGRAM << ": output increment must be smaller than number of time steps.  Use\n\n";
					std::cout << "    " << PROGRAM << " --help\n\n";
					std::cout << "to generate help message.\n\n";
					exit(-1);
				}
			}
		}

		else {
			// set output file name
			outfile = argv[2];

			// set number of time steps
			if (std::string(argv[3]).find_first_not_of("0123456789") != std::string::npos) {
				// must have integral number of time steps
				std::cout << PROGRAM << ": number of time steps must have integral value.  Use\n\n";
				std::cout << "    " << PROGRAM << " --help\n\n";
				std::cout << "to generate help message.\n\n";
				exit(-1);
			}

			steps = atoi(argv[3]);
			increment = steps;

			if (argc > 4) {
				// must have integral output increment
				if (std::string(argv[4]).find_first_not_of("0123456789") != std::string::npos) {
					std::cout << PROGRAM << ": output increment must have integral value.  Use\n\n";
					std::cout << "    " << PROGRAM << " --help\n\n";
					std::cout << "to generate help message.\n\n";
					exit(-1);
				}

				increment = atoi(argv[4]);

				// output increment must be smaller than number of steps
				if (increment > steps) {
					std::cout << PROGRAM << ": output increment must be smaller than number of time steps.  Use\n\n";
					std::cout << "    " << PROGRAM << " --help\n\n";
					std::cout << "to generate help message.\n\n";
					exit(-1);
				}
			}
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

		// grid type error check
		if (type.substr(0, 4) != "grid") {
			std::cerr << "File input error: file does not contain grid data." << std::endl;
			exit(-1);
		}

		// read grid dimension
		int dim;
		input >> dim;

		// set output file basename
		int iterations_start(0);
		if (outfile.find_first_of(".") != outfile.find_last_of(".")) {
			std::string number = outfile.substr(outfile.find_first_of(".") + 1, outfile.find_last_of(".") - 1);
			iterations_start = atoi(number.c_str());
		}
		std::string base;
		if (outfile.rfind(".", outfile.find_last_of(".") - 1) == -1) // only one dot found
			base = outfile.substr(0, outfile.find_last_of(".")) + ".";
		else {
			int last_dot = outfile.find_last_of(".");
			int prev_dot = outfile.rfind('.', last_dot - 1);
			std::string number = outfile.substr(prev_dot + 1, last_dot - prev_dot - 1);
			bool isNumeric(true);
			for (int i = 0; i < number.size(); ++i) {
				if (!isdigit(number[i])) isNumeric = false;
			}
			if (isNumeric)
				base = outfile.substr(0, outfile.rfind(".", outfile.find_last_of(".") - 1)) + ".";
			else base = outfile.substr(0, outfile.find_last_of(".")) + ".";
		}

		// set output file suffix
		std::string suffix = "";
		if (outfile.find_last_of(".") != std::string::npos)
			suffix = outfile.substr(outfile.find_last_of("."), std::string::npos);

		// set output filename length
		int length = base.length() + suffix.length();
		if (1) {
			std::stringstream slength;
			slength << steps;
			length += slength.str().length();
		}


		if (dim == 2) {
			// construct grid object
			GRID2D grid(argv[1]);

			// perform computation
			for (int i = iterations_start; i < steps; i += increment) {
				MMSP::update(grid, increment);

				// generate output filename
				std::stringstream outstr;
				int n = outstr.str().length();
				for (int j = 0; n < length; j++) {
					outstr.str("");
					outstr << base;
					for (int k = 0; k < j; k++) outstr << 0;
					outstr << i + increment << suffix;
					n = outstr.str().length();
				}

				// write grid output to file
				char* filename = new char[outfile.length()];
				for (unsigned int i=0; i<outfile.length(); i++)
					filename[i] = outfile[i];
				#ifdef BGQ
				MMSP::output_bgq(grid, filename);
				#else
				MMSP::output(grid, filename);
				#endif
				delete [] filename;
			}
		}

		if (dim == 3) {
			// construct grid object
			GRID3D grid(argv[1]);

			// perform computation
			for (int i = iterations_start; i < steps; i += increment) {
				MMSP::update(grid, increment);

				// generate output filename
				std::stringstream outstr;
				int n = outstr.str().length();
				for (int j = 0; n < length; j++) {
					outstr.str("");
					outstr << base;
					for (int k = 0; k < j; k++) outstr << 0;
					outstr << i + increment << suffix;
					n = outstr.str().length();
				}

				// write grid output to file
				char* filename = new char[outfile.length()];
				for (unsigned int i=0; i<outfile.length(); i++)
					filename[i] = outfile[i];
				#ifdef BGQ
				MMSP::output_bgq(grid, filename);
				#else
				MMSP::output(grid, filename);
				#endif
				delete [] filename;
			}
		}
	}

	MMSP::Finalize();
}

#endif
