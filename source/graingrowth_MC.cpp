// graingrowth.hpp
// Algorithms for 2D and 3D isotropic Monte Carlo grain growth
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#ifndef GRAINGROWTH_UPDATE
#define GRAINGROWTH_UPDATE
#include <iomanip>
#include <vector>
#include <cmath>
#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <pthread.h>
#include "rdtsc.h"
#include"graingrowth_MC.hpp"
#include"MMSP.hpp"
#include"tessellate.hpp"
#include"output.cpp"

void print_progress(const int step, const int steps, const int iterations);

namespace MMSP
{
template <int dim>
MMSP::grid<dim,int >* generate(int seeds, int nthreads)
{
	#if (defined CCNI) && (!defined MPI_VERSION)
	std::cerr<<"Error: MPI is required for CCNI."<<std::endl;
	exit(1);
	#endif
	#ifdef MPI_VERSION
	int np = MPI::COMM_WORLD.Get_size();
	#endif

	if (dim == 2) {
		const int edge = 1024;
		int number_of_fields(seeds);
		if (number_of_fields==0) number_of_fields = static_cast<int>(float(edge*edge)/(M_PI*10.*10.)); // average grain is a disk of radius 10
		#ifdef MPI_VERSION
		while (number_of_fields % np) --number_of_fields;
		#endif
		MMSP::grid<dim,int >* grid = new MMSP::grid<dim,int>(0, 0, edge, 0, edge);

		#ifdef MPI_VERSION
		number_of_fields /= np;
		#endif

		#if (!defined MPI_VERSION) && ((defined CCNI) || (defined BGQ))
		std::cerr<<"Error: CCNI requires MPI."<<std::endl;
		std::exit(1);
		#endif
		tessellate<dim,int>(*grid, number_of_fields, nthreads);
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
		MMSP::grid<dim,int >* grid = new MMSP::grid<dim,int>(0, 0, edge, 0, edge, 0, edge);

		#ifdef MPI_VERSION
		number_of_fields /= np;
		#endif

		tessellate<dim,int >(*grid, number_of_fields, nthreads);
		#ifdef MPI_VERSION
		MPI::COMM_WORLD.Barrier();
		#endif
		return grid;
	}
	return NULL;
}

void generate(int dim, char* filename, int seeds, int nthreads)
{
	#if (defined CCNI) && (!defined MPI_VERSION)
	std::cerr<<"Error: MPI is required for CCNI."<<std::endl;
	exit(1);
	#endif
	int rank=0;
	#ifdef MPI_VERSION
	rank = MPI::COMM_WORLD.Get_rank();
	#endif
	if (dim == 2) {
		MMSP::grid<2,int>* grid2=generate<2>(seeds,nthreads);
		assert(grid2!=NULL);
		#ifdef BGQ
		output_bgq(*grid2, filename);
		#else
		output(*grid2, filename);
		#endif
		#ifndef SILENT
		if (rank==0) std::cout<<"Wrote initial file to "<<filename<<"."<<std::endl;
		#endif
	}

	if (dim == 3) {
		MMSP::grid<3,int>* grid3=generate<3>(seeds,nthreads);
		assert(grid3!=NULL);
		#ifdef BGQ
		output_bgq(*grid3, filename);
		#else
		output(*grid3, filename);
		#endif
		#ifndef SILENT
		if (rank==0) std::cout<<"Wrote initial file to "<<filename<<"."<<std::endl;
		#endif
	}

}



template <int dim> struct flip_index {
	MMSP::grid<dim, int>* grid;
	int front_low_left_corner[dim];
	int back_up_right_corner[dim];
	int sublattice;
};


template <int dim> void* flip_index_helper( void* s )
{
	flip_index<dim>* ss = static_cast<flip_index<dim>*>(s);
	int sublattice = ss->sublattice;
	vector<int> x (dim,0);

	const double kT = 0.50;

	// choose a random node
	srand(time(NULL)); /* seed random number generator */
	int range;

	range=(ss->back_up_right_corner[0]-ss->front_low_left_corner[0]);

	int num_of_grids = (ss->back_up_right_corner[0]-ss->front_low_left_corner[0]);
	for (int k=1; k<dim; k++)
		num_of_grids *= (ss->back_up_right_corner[k]-ss->front_low_left_corner[k]);

	for (int h=0; h<num_of_grids; h++) {
		if (range%2==0) { //range is even

			if (sublattice==0) {
				x[0] = ss->front_low_left_corner[0] + 2*(rand()%(range/2+1));  //even
			} else
				x[0] = ss->front_low_left_corner[0] + 2*(rand()%(range/2)) + 1; //odd
		} else { //range is odd, can only happen at the last pthread
			if (sublattice==0)
				x[0] = ss->front_low_left_corner[0]+2*(rand()%((range+1)/2));  //even
			else
				x[0] = ss->front_low_left_corner[0]+2*(rand()%((range+1)/2))+1; //odd
		}

		for (int k=1; k<dim; k++)
			x[k] = ss->front_low_left_corner[k]+rand()%(ss->back_up_right_corner[k]-ss->front_low_left_corner[k]+1);


		int spin1 = (*(ss->grid))(x);

		// determine neighboring spins
		vector<int> r(x);
		sparse<bool> neighbors;
		for (int i=-1; i<=1; i++) {
			for (int j=-1; j<=1; j++) {
				r[0] = x[0] + i;
				r[1] = x[1] + j;
				int spin = (*(ss->grid))(r);
				set(neighbors,spin) = true;
			}
		}

		// choose a random neighbor spin
		int spin2 = index(neighbors,rand()%length(neighbors));

		if (spin1!=spin2) {
			// compute energy change
			double dE = -1.0;
			for (int i=-1; i<=1; i++) {
				for (int j=-1; j<=1; j++) {
					r[0] = x[0] + i;
					r[1] = x[1] + j;
					int spin = (*(ss->grid))(r);
					dE += (spin!=spin2)-(spin!=spin1);
				}
			}

			// attempt a spin flip
			double r = double(rand())/double(RAND_MAX);
			if (dE<=0.0) (*(ss->grid))(x) = spin2;
			else if (r<exp(-dE/kT)) (*(ss->grid))(x) = spin2;
		}
	}
	pthread_exit(0);
	return NULL;
}


template <int dim> void update(MMSP::grid<dim, int>& grid, int steps, int nthreads)
{
	#if (!defined MPI_VERSION) && ((defined CCNI) || (defined BGQ))
	std::cerr<<"Error: MPI is required for CCNI."<<std::endl;
	exit(1);
	#endif
	int rank=0;
	#ifdef MPI_VERSION
	rank=MPI::COMM_WORLD.Get_rank();
	#endif

	MPI::COMM_WORLD.Barrier();

	unsigned long start = rdtsc();

	int front_low_left_corner[dim];
	int back_up_right_corner[dim];

	pthread_t* p_threads = new pthread_t[nthreads];
	flip_index<dim>* mat_para = new flip_index<dim> [nthreads];
	pthread_attr_t attr;
	pthread_attr_init (&attr);

	//check if num of the pthread is too large, if so, reduce it.
	if ((x1(grid, 0)-x0(grid, 0)-nthreads-1)/nthreads<1) {
		std::cerr<<"ERROR: number of pthread is too large, please reduce it to a value <= "<<x1(grid, 0)-x0(grid, 0)-nthreads-1<<std::endl;
		exit(0);
	}

	#ifndef SILENT
	static int iterations = 1;
	if (rank==0) print_progress(0, steps, iterations);
	#endif
	for (int step=0; step<steps; step++) {
		for (int sublattice=0; sublattice!= 2; sublattice++) {
			for (int i=0; i!= nthreads ; i++ ) {

				front_low_left_corner[0] = x0(grid, 0) + ((x1(grid, 0)-x0(grid, 0)-nthreads-1)/nthreads)*i+i;
				for (int k=1; k<dim; k++)
					front_low_left_corner[k] = x0(grid, k);

				if (i==(nthreads-1))
					back_up_right_corner[0] = x1(grid, 0);
				else
					back_up_right_corner[0] = x0(grid, 0) + ((x1(grid, 0)-x0(grid, 0)-nthreads-1)/nthreads)*(i+1)+i;
				for (int k=1; k<dim; k++)
					back_up_right_corner[k] = x1(grid, k);

				mat_para[i].grid = &grid;
				for (int jj=0; jj<dim; jj++) {
					mat_para[i].front_low_left_corner[jj] = front_low_left_corner[jj];
					mat_para[i].back_up_right_corner[jj] = back_up_right_corner[jj];
				}
				mat_para[i].sublattice= sublattice;

				pthread_create(&p_threads[i], &attr, flip_index_helper<dim>, (void*) &mat_para[i] );

			}//loop over threads

			for (int i=0; i!= nthreads ; i++)
				pthread_join(p_threads[i], NULL);

			MPI::COMM_WORLD.Barrier();

			ghostswap(grid); // once looped over a "color", ghostswap.
		}//loop over color
		#ifndef SILENT
		if (rank==0) print_progress(step+1, steps, iterations);
		#endif
	}//loop over step
	#ifndef SILENT
	++iterations;
	#endif

	delete [] p_threads ;
	p_threads=NULL;
	delete [] mat_para ;
	mat_para=NULL;

    unsigned long update_timer = rdtsc()-start;
    unsigned long total_update_time;
    MPI::COMM_WORLD.Allreduce(&update_timer, &total_update_time, 1, MPI_UNSIGNED_LONG, MPI_SUM);
    if(rank==0) std::cout<<"Monte Carlo total update time is "<<total_update_time<<std::endl;
}

}

#ifndef SILENT
void print_progress(const int step, const int steps, const int iterations)
{
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
	} else if (step==steps) {
		unsigned long deltat = time(NULL)-tstart;
		std::cout << "•] "
							<<std::setw(2)<<std::right<<deltat/3600<<"h:"
							<<std::setw(2)<<std::right<<(deltat%3600)/60<<"m:"
							<<std::setw(2)<<std::right<<deltat%60<<"s"
							<<" (File "<<std::setw(5)<<std::right<<iterations*steps<<")."<<std::endl;
	} else if ((20 * step) % steps == 0) std::cout<<"• "<<std::flush;
}
#endif

#endif

#include"MMSP.main.hpp"

// Formatted using astyle:
//  astyle --style=linux --indent-col1-comments --indent=tab --indent-preprocessor --pad-header --align-pointer=type --keep-one-line-blocks --suffix=none
