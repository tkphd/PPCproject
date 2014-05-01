// graingrowth.cpp
// Coarsening algorithms for 2D and 3D sparse phase field (sparsePF) methods
// Questions/comments to kellet@rpi.edu (Trevor Keller)

#ifndef GRAINGROWTH_UPDATE
#define GRAINGROWTH_UPDATE

#include <iomanip>
#include <vector>
#include <cmath>
#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <cassert>

#include"graingrowth.hpp"
#include"MMSP.hpp"
#include"tessellate.hpp"
#include"output.cpp"

void print_progress(const int step, const int steps, const int iterations);

namespace MMSP {

template <int dim>
MMSP::grid<dim,MMSP::sparse<float> >* generate(int seeds, int nthreads)
{
	#if (defined CCNI) && (!defined MPI_VERSION)
	std::cerr<<"Error: MPI is required for CCNI."<<std::endl;
	exit(1);
	#endif
	int rank=0;
	#ifdef MPI_VERSION
	rank = MPI::COMM_WORLD.Get_rank();
	int np = MPI::COMM_WORLD.Get_size();
	#endif
	if (dim == 2) {
		const int edge = 1024;
		int number_of_fields(seeds);
		if (number_of_fields==0) number_of_fields = static_cast<int>(float(edge*edge)/(M_PI*10.*10.)); // average grain is a disk of radius 10
		#ifdef MPI_VERSION
		while (number_of_fields % np) --number_of_fields;
		#endif
		MMSP::grid<dim,MMSP::sparse<float> >* grid = new MMSP::grid<dim,MMSP::sparse<float> >(0, 0, edge, 0, edge);
		if (rank==0) std::cout<<"Grid origin: ("<<g0(*grid,0)<<','<<g0(*grid,1)<<"),"
												<<" dimensions: "<<g1(*grid,0)-g0(*grid,0)<<" × "<<g1(*grid,1)-g0(*grid,1)
												<<" with "<<number_of_fields<<" grains."<<std::endl;
		#ifdef MPI_VERSION
		number_of_fields /= np;
		#endif

		#if (!defined MPI_VERSION) && ((defined CCNI) || (defined BGQ))
		std::cerr<<"Error: CCNI requires MPI."<<std::endl;
		std::exit(1);
		#endif
		tessellate<dim,float>(*grid, number_of_fields);
		if (rank==0) std::cout<<"Tessellation complete."<<std::endl;
		#ifdef MPI_VERSION
		MPI::COMM_WORLD.Barrier();
		#endif
		return grid;
	} else if (dim == 3) {
		const int edge = 512;
		int number_of_fields(seeds);
		if (number_of_fields==0) number_of_fields = static_cast<int>(float(edge*edge*edge)/(4./3*M_PI*10.*10.*10.)); // Average grain is a sphere of radius 10 voxels
		#ifdef MPI_VERSION
		while (number_of_fields % np) --number_of_fields;
		#endif
		MMSP::grid<dim,MMSP::sparse<float> >* grid = new MMSP::grid<dim,MMSP::sparse<float> >(0,0,edge,0,edge,0,edge);
		if (rank==0) std::cout<<"Grid origin: ("<<g0(*grid,0)<<','<<g0(*grid,1)<<','<<g0(*grid,2)<<"),"
												<<" dimensions: "<<g1(*grid,0)-g0(*grid,0)<<" × "<<g1(*grid,1)-g0(*grid,1)<<" × "<<g1(*grid,2)-g0(*grid,2)
												<<" with "<<number_of_fields<<" grains."<<std::endl;
		#ifdef MPI_VERSION
		number_of_fields /= np;
		#endif

		#if (!defined MPI_VERSION) && ((defined CCNI) || (defined BGQ))
		std::cerr<<"Error: CCNI requires MPI."<<std::endl;
		std::exit(1);
		#endif
		tessellate<dim,float>(*grid, number_of_fields);
		#ifdef MPI_VERSION
		MPI::COMM_WORLD.Barrier();
		#endif
		if (rank==0) std::cout<<"Tessellation complete."<<std::endl;
		return grid;
	}
	return NULL;
}


void generate(int dim, char* filename, int seeds, int nthreads) {
	#if (defined CCNI) && (!defined MPI_VERSION)
	std::cerr<<"Error: MPI is required for CCNI."<<std::endl;
	exit(1);
	#endif
	int rank=0;
	#ifdef MPI_VERSION
	rank = MPI::COMM_WORLD.Get_rank();
	#endif
	if (dim == 2) {
		MMSP::grid<2,MMSP::sparse<float> >* grid2=generate<2>(seeds,nthreads);
		assert(grid2!=NULL);
		#ifdef BGQ
		output_bgq(*grid2, filename);
		#else
		output(*grid2, filename);
		#endif
		if (rank==0) std::cout<<"Wrote initial file to "<<filename<<"."<<std::endl;
	}

	if (dim == 3) {
		MMSP::grid<3,MMSP::sparse<float> >* grid3=generate<3>(seeds,nthreads);
		assert(grid3!=NULL);
		#ifdef BGQ
		output_bgq(*grid3, filename);
		#else
		output(*grid3, filename);
		#endif
		if (rank==0) std::cout<<"Wrote initial file to "<<filename<<"."<<std::endl;
	}
}

template <int dim> void update(MMSP::grid<dim, sparse<float> >& grid, int steps, int nthreads) {
	#if (!defined MPI_VERSION) && ((defined CCNI) || (defined BGQ))
	std::cerr<<"Error: MPI is required for CCNI."<<std::endl;
	exit(1);
	#endif
	int rank=0;
	#ifdef MPI_VERSION
 	rank=MPI::COMM_WORLD.Get_rank();
	#endif
	const float dt = 0.01;
	const float width = 10.0;
	const float gamma = 1.0;
	const float eps = 4.0 / acos(-1.0) * sqrt(0.5 * gamma * width);
	const float w = 4.0 * gamma / width;
	const float mu = 1.0;
	const float epsilon = 1.0e-8;

	static int iterations = 1;
	#ifdef DEBUG
	if (iterations==1 && rank==0)
		printf("CFL condition Co=%2.2f.\n", mu*eps*eps*dt/(dx(grid, 0)*dx(grid,0)));
	#endif

 	if (rank==0) print_progress(0, steps, iterations);

	for (int step = 0; step < steps; step++) {
		// update grid must be overwritten each time
		MMSP::grid<dim, sparse<float> > update(grid);
		ghostswap(grid);

		for (int i = 0; i < nodes(grid); i++) {
			vector<int> x = position(grid, i);

			// determine nonzero fields within
			// the neighborhood of this node
			// (2 adjacent voxels along each cardinal direction)
			sparse<int> s;
			for (int j = 0; j < dim; j++)
				for (int k = -1; k <= 1; k++) {
					x[j] += k;
					for (int h = 0; h < length(grid(x)); h++) {
						int index = MMSP::index(grid(x), h);
						set(s, index) = 1;
					}
					x[j] -= k;
				}
			float S = float(length(s));

			// if only one field is nonzero,
			// then copy this node to update
			if (S < 2.0) update(i) = grid(i);
			else {
				// compute laplacian of each field
				sparse<float> lap = laplacian(grid, i);

				// compute variational derivatives
				sparse<float> dFdp;
				for (int h = 0; h < length(s); h++) {
					int hindex = MMSP::index(s, h);
					for (int j = h + 1; j < length(s); j++) {
						int jindex = MMSP::index(s, j);
						// Update dFdp_h and dFdp_j, so the inner loop can be over j>h instead of j≠h
						set(dFdp, hindex) += 0.5 * eps * eps * lap[jindex] + w * grid(i)[jindex];
						set(dFdp, jindex) += 0.5 * eps * eps * lap[hindex] + w * grid(i)[hindex];
					}
				}

				// compute time derivatives
				sparse<float> dpdt;
				for (int h = 0; h < length(s); h++) {
					int hindex = MMSP::index(s, h);
					for (int j = h + 1; j < length(s); j++) {
						int jindex = MMSP::index(s, j);
						set(dpdt, hindex) -= mu * (dFdp[hindex] - dFdp[jindex]);
						set(dpdt, jindex) -= mu * (dFdp[jindex] - dFdp[hindex]);
					}
				}

				// compute update values
				float sum = 0.0;
				for (int h = 0; h < length(s); h++) {
					int index = MMSP::index(s, h);
					float value = grid(i)[index] + dt * (2.0 / S) * dpdt[index]; // Extraneous factor of 2?
					if (value > 1.0) value = 1.0;
					if (value < 0.0) value = 0.0;
					if (value > epsilon) set(update(i), index) = value;
					sum += update(i)[index];
				}

				// project onto Gibbs simplex (enforce Σφ=1)
				float rsum = 0.0;
				if (fabs(sum) > 0.0) rsum = 1.0 / sum;
				for (int h = 0; h < length(update(i)); h++) {
					int index = MMSP::index(update(i), h);
					set(update(i), index) *= rsum;
				}
			}
		} // Loop over nodes(grid)
		if (rank==0) print_progress(step+1, steps, iterations);
		swap(grid, update);
	} // Loop over steps
	ghostswap(grid);
	++iterations;
}

template <class T> std::ostream& operator<<(std::ostream& o, sparse<float>& s) {
	o<<"    Index	Value\n";
	for (int i=0; i<length(s); ++i) {
		int index = MMSP::index(s, i);
		o<<"    "<<std::setw(5)<<std::right<<index<<"  "<<s[index]<<'\n';
	}
	return o;
}

} // namespace MMSP

void print_progress(const int step, const int steps, const int iterations) {
	char* timestring;
	static unsigned long tstart;
	struct tm* timeinfo;

	if (step==0) {
		tstart = time(NULL);
		std::time_t rawtime;
		std::time( &rawtime );
		timeinfo = std::localtime( &rawtime );
		timestring = std::asctime(timeinfo);
		timestring[std::strlen(timestring)-1] = '\0';
		std::cout<<"Pass "<<std::setw(3)<<std::right<<iterations<<": "<<timestring<<" ["<<std::flush;
	} else if (step==steps-1) {
		unsigned long deltat = time(NULL)-tstart;
		std::cout << "•] "
							<<std::setw(2)<<std::right<<deltat/3600<<"h:"
							<<std::setw(2)<<std::right<<(deltat%3600)/60<<"m:"
							<<std::setw(2)<<std::right<<deltat%60<<"s"
							<<" (File "<<std::setw(5)<<std::right<<iterations*steps<<")."<<std::endl;
	} else if ((20 * step) % steps == 0) std::cout<<"• "<<std::flush;
}

#endif

#include"MMSP.main.hpp"
