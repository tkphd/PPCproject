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

#define MPI_VERSION /*-----------commets out when on BGQ------------*/


void print_progress(const int step, const int steps, const int iterations);

namespace MMSP
{
template <int dim>
unsigned long generate(MMSP::grid<dim,int >*& grid, int seeds, int nthreads)
{
	#if (defined CCNI) && (!defined MPI_VERSION)
	std::cerr<<"Error: MPI is required for CCNI."<<std::endl;
	exit(1);
	#endif
	#ifdef MPI_VERSION
	int np = MPI::COMM_WORLD.Get_size();
	#endif

	unsigned long timer=0;
	if (dim == 2) {
		const int edge = 1024;
		int number_of_fields(seeds);
		if (number_of_fields==0) number_of_fields = static_cast<int>(float(edge*edge)/(M_PI*10.*10.)); // average grain is a disk of radius 10
		#ifdef MPI_VERSION
		while (number_of_fields % np) --number_of_fields;
		#endif
		grid = new MMSP::grid<dim,int>(0, 0, edge, 0, edge);

		#ifdef MPI_VERSION
		number_of_fields /= np;
		#endif

		#if (!defined MPI_VERSION) && ((defined CCNI) || (defined BGQ))
		std::cerr<<"Error: CCNI requires MPI."<<std::endl;
		std::exit(1);
		#endif
		timer = tessellate<dim,int>(*grid, number_of_fields, nthreads);
		#ifdef MPI_VERSION
		MPI::COMM_WORLD.Barrier();
		#endif
	} else if (dim == 3) {
		const int edge = 8;
		int number_of_fields(seeds);
		if (number_of_fields==0) number_of_fields = static_cast<int>(float(edge*edge*edge)/(4./3*M_PI*10.*10.*10.)); // Average grain is a sphere of radius 10 voxels
		#ifdef MPI_VERSION
		while (number_of_fields % np) --number_of_fields;
		#endif
		grid = new MMSP::grid<dim,int>(0, 0, edge, 0, edge, 0, edge);

		#ifdef MPI_VERSION
		number_of_fields /= np;
		#endif

		timer = tessellate<dim,int >(*grid, number_of_fields, nthreads);
		#ifdef MPI_VERSION
		MPI::COMM_WORLD.Barrier();
		#endif
	}
	return timer;
}

unsigned long generate(int dim, char* filename, int seeds, int nthreads)
{
	#if (defined CCNI) && (!defined MPI_VERSION)
	std::cerr<<"Error: MPI is required for CCNI."<<std::endl;
	exit(1);
	#endif
	int rank=0;
	#ifdef MPI_VERSION
	rank = MPI::COMM_WORLD.Get_rank();
	#endif

	unsigned long timer = 0;
	if (dim == 2) {
		MMSP::grid<2,int>* grid2=NULL;
		timer = generate<2>(grid2,seeds,nthreads);
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
		MMSP::grid<3,int>* grid3=NULL;
		timer = generate<3>(grid3,seeds,nthreads);
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
	return timer;
}

int LargeNearestInteger(int a, int b){
  if(a%b==0) return a/b;
  else return a/b+1;
}


template <int dim> struct flip_index {
	MMSP::grid<dim, int>* grid;
  int num_of_cells_in_thread;
	int sublattice;
  int cell_coord[dim];
  int lattice_cells_each_dimension[dim];
};


template <int dim> void* flip_index_helper( void* s )
{
  srand(time(NULL)); /* seed random number generator */
	flip_index<dim>* ss = static_cast<flip_index<dim>*>(s);

	vector<int> x (dim,0);
	const double kT = 0.50;
  int first_cell_start_coordinates[dim];
  for(int k=0; k<dim; k++) first_cell_start_coordinates[k] = x0(*(ss->grid), k);
  for(int i=0; i<dim; i++){
    if(x0(*(ss->grid), i)%2!=0) first_cell_start_coordinates[i]--;
  }

	int num_of_grids = (ss->num_of_cells_in_thread)*pow(2,dim);
  int cell_coords_selected[dim];
	for (int h=0; h<num_of_grids; h++) {
	  // choose a random cell to flip
    int cell_numbering_in_thread = rand()%(ss->num_of_cells_in_thread); //choose a cell to flip, from 0 to num_of_cells_in_thread-1
    if(dim==2){
      cell_coords_selected[1]=((ss->cell_coord)[1]+cell_numbering_in_thread)%(ss->lattice_cells_each_dimension)[1];//1-indexed
      cell_coords_selected[0]=(ss->cell_coord)[0]+(((ss->cell_coord)[1]+cell_numbering_in_thread)/(ss->lattice_cells_each_dimension)[0]);
    }else if(dim==3){
      cell_coords_selected[2]=((ss->cell_coord)[2]+cell_numbering_in_thread)%(ss->lattice_cells_each_dimension)[2];//1-indexed
      cell_coords_selected[1]=(ss->cell_coord)[1]+(((ss->cell_coord)[2]+cell_numbering_in_thread)/(ss->lattice_cells_each_dimension)[2])%(ss->lattice_cells_each_dimension)[1];
      cell_coords_selected[0]=(ss->cell_coord)[0]+(((ss->cell_coord)[1]+((ss->cell_coord)[2]+cell_numbering_in_thread)/(ss->lattice_cells_each_dimension)[2]))/(ss->lattice_cells_each_dimension)[1];
    }
 
    for(int i=0; i<dim; i++){
      x[i]=first_cell_start_coordinates[i]+2*cell_coords_selected[i];
    }


/*    #ifdef MPI_VERSION
      int rank=MPI::COMM_WORLD.Get_rank(); 
      if(rank==0 && dim==3)
	std::cout<<"coordinates are (before consider boundary)"<<x[0]<<"  "<<x[1]<<"  "<<x[2]<<std::endl;
     #endif
*/

    if(dim==2){
      switch(ss->sublattice){
        case 1:break;// 0,0
        case 2:x[1]++; break; //0,1
        case 3:x[0]++; break; //1,0
        case 4:x[0]++; x[1]++; break; //1,1
      }
    }
    if(dim==3){
      switch(ss->sublattice){
        case 1:break;// 0,0,0
        case 2:x[2]++; break; //0,0,1
        case 3:x[1]++; break; //0,1,0
        case 4:x[2]++; x[1]++; break; //0,1,1
        case 5:x[0]++; break; //1,0,0
        case 6:x[2]++; x[0]++; break; //1,0,1
        case 7:x[1]++; x[0]++; break; //1,1,0
        case 8:x[2]++; x[1]++; x[0]++; //1,1,1
      }
    }


      //      if(dim==3 && rank==0) std::cout<<"cell_coord is "<<(ss->cell_coord)[0]<<"  "<<(ss->cell_coord)[1]<<"  "<<(ss->cell_coord)[2]<<"\n";
      // if(dim==3 && rank==0) std::cout<<"cell_coords_selected is "<<cell_coords_selected[0]<<"  "<<cell_coords_selected[1]<<"  "<<cell_coords_selected[2]<<"\n";
/*      if(dim==3 && rank==0){
      std::cout<<"sublattice is "<<ss->sublattice<<"\n";
      std::cout<<"cell_coords_selected is "<<cell_coords_selected[0]<<"  "<<cell_coords_selected[1]<<"  "<<cell_coords_selected[2]<<"\n";  
      std::cout<<"coordinates are "<<x[0]<<"  "<<x[1]<<"  "<<x[2]<<std::endl;
      }*/


  
    bool site_out_of_domain = false;
    for(int i=0; i<dim; i++)
      if(x[i]<x0(*(ss->grid), i) || x[i]>x1(*(ss->grid), i)){
        site_out_of_domain = true;
        break;//break from the for int i loop
      }
    if(site_out_of_domain == true)
      continue; //continue the int h loop

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

template <int dim> unsigned long update(MMSP::grid<dim, int>& grid, int steps, int nthreads)
{
	#if (!defined MPI_VERSION) && ((defined CCNI) || (defined BGQ))
	std::cerr<<"Error: MPI is required for CCNI."<<std::endl;
	exit(1);
	#endif
	int rank=0;
	unsigned int np=0;
//	#ifdef MPI_VERSION

	rank=MPI::COMM_WORLD.Get_rank();
	np=MPI::COMM_WORLD.Get_size();

	MPI::COMM_WORLD.Barrier();
//	#endif

	unsigned long update_timer = 0;

	pthread_t* p_threads = new pthread_t[nthreads];
	flip_index<dim>* mat_para = new flip_index<dim> [nthreads];
	pthread_attr_t attr;
	pthread_attr_init (&attr);

	#ifndef SILENT
	static int iterations = 1;
	if (rank==0) print_progress(0, steps, iterations);
	#endif


/*-----------------------------------------------*/
/*---------------generate cells------------------*/
/*-----------------------------------------------*/
  int dimension_length=0, number_of_lattice_cells=1;
  int lattice_cells_each_dimension[dim];
  for(int i=0; i<dim; i++){
    dimension_length = x1(grid, i)-x0(grid, i);
    if(x0(grid, 0)%2==0)
      lattice_cells_each_dimension[i] = dimension_length/2+1;
    else
      lattice_cells_each_dimension[i] = 1+LargeNearestInteger(dimension_length,2);
    number_of_lattice_cells *= lattice_cells_each_dimension[i];
  }
  if(rank==0){ 
    std::cout<<"lattice_cells_each_dimension is ";
    for(int i=0; i<dim; i++)
      std::cout<<lattice_cells_each_dimension[i]<<"  ";
    std::cout<<"\n";
  }

  
//----------assign cells for each pthreads
  int num_of_cells_in_thread = number_of_lattice_cells/nthreads;
  //check if num of the pthread is too large, if so, reduce it.                                                                                                                                         
  if (num_of_cells_in_thread<1) {
    std::cerr<<"ERROR: number of pthread is too large, please reduce it to a value <= "<<number_of_lattice_cells<<std::endl;
    exit(0);
  }


  if(rank==0)
  std::cout<<"number_of_lattice_cells is "<<number_of_lattice_cells<<"\n";

  if(rank==0)
  std::cout<<"nthreads is "<<nthreads<<"\n";

  if(rank==0)
    std::cout<<x0(grid, 0)<<","<<x0(grid,1)<<","<<x0(grid, 2)<<" ->  "<<x1(grid, 0)<<","<<x1(grid,1)<<","<<x1(grid, 2)<<std::endl;

  for (int i=0; i<nthreads; i++) {
    mat_para[i].grid = &grid;
    if(i==(nthreads-1)) mat_para[i].num_of_cells_in_thread = number_of_lattice_cells - num_of_cells_in_thread*(nthreads-1);
    else mat_para[i].num_of_cells_in_thread = num_of_cells_in_thread;
    if(rank==0) std::cout<<"num_of_cells_in_thread is "<<mat_para[i].num_of_cells_in_thread<<" in thread "<<i<<"\n";
    for(int k=0; k<dim; k++) mat_para[i].lattice_cells_each_dimension[k]=lattice_cells_each_dimension[k];
  }


  int cell_coord[dim];//record the start coordinates of each pthread domain.


	for (int step=0; step<steps; step++){
		unsigned long start = rdtsc();
    int num_of_sublattices = 8;
		for (int sublattice=0; sublattice!= num_of_sublattices; sublattice++) {
			for (int i=0; i!= nthreads ; i++) {
        int cell_numbering = num_of_cells_in_thread*i; //0-indexed, celling_numbering is the start cell numbering
        if(rank==0) std::cout<<"cell_numbering is "<<cell_numbering<<"\n";
        if(dim==2){
          cell_coord[1]=cell_numbering%lattice_cells_each_dimension[1];//0-indexed
          cell_coord[0]=(cell_numbering/lattice_cells_each_dimension[1]);
	  if(cell_coord[0]>=lattice_cells_each_dimension[0]){
	    std::cerr<<"the cell coordinates is wrong!"<<std::endl;
	    exit(1);
	  }
        }else if(dim==3){
          cell_coord[2]=cell_numbering%lattice_cells_each_dimension[2];//0-indexed
          cell_coord[1]=(cell_numbering/lattice_cells_each_dimension[2])%lattice_cells_each_dimension[1];
          cell_coord[0]=(cell_numbering/lattice_cells_each_dimension[2])/lattice_cells_each_dimension[1];
	  if(cell_coord[0]>=lattice_cells_each_dimension[0]){
	    std::cerr<<"the cell coordinates is wrong!"<<std::endl;
	    exit(1);
	  }
        }
	// if(dim==3 && rank==0) std::cout<<"cell_coord is "<<cell_coord[0]<<"  "<<cell_coord[1]<<"  "<<cell_coord[2]<<"\n";
				mat_para[i].sublattice=sublattice;
        for(int k=0; k<dim; k++) mat_para[i].cell_coord[k]=cell_coord[k];
				pthread_create(&p_threads[i], &attr, flip_index_helper<dim>, (void*) &mat_para[i] );

			}//loop over threads

			for (int i=0; i!= nthreads ; i++)
				pthread_join(p_threads[i], NULL);

			#ifdef MPI_VERSION
			MPI::COMM_WORLD.Barrier();
			#endif

			ghostswap(grid); // once looped over a "color", ghostswap.
		}//loop over color
		#ifndef SILENT
		if (rank==0) print_progress(step+1, steps, iterations);
		#endif
		update_timer += rdtsc()-start;
	}//loop over step
	#ifndef SILENT
	++iterations;
	#endif

	delete [] p_threads ;
	p_threads=NULL;
	delete [] mat_para ;
	mat_para=NULL;

	unsigned long total_update_time=update_timer;
	#ifdef MPI_VERSION
	MPI::COMM_WORLD.Allreduce(&update_timer, &total_update_time, 1, MPI_UNSIGNED_LONG, MPI_SUM);
	#endif
	//if(rank==0) std::cout<<"Monte Carlo total update time is "<<total_update_time<<std::endl;
	return total_update_time/np; // average update time
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
