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

#include"MMSP.hpp"
#include"tessellate.hpp"
#include"graingrowth.hpp"
#include"output.cpp"


namespace MMSP{

void generate(int dim, char* filename, int seeds=0){
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
		MMSP::grid<2,int > grid(0,0,edge,0,edge);
		if (rank==0) std::cout<<"Grid origin: ("<<g0(grid,0)<<','<<g0(grid,1)<<"),"
												<<" dimensions: "<<g1(grid,0)-g0(grid,0)<<" × "<<g1(grid,1)-g0(grid,1)
												<<" with "<<number_of_fields<<" grains."<<std::endl;
		#ifdef MPI_VERSION
		number_of_fields /= np;
		#endif

		#if (!defined MPI_VERSION) && ((defined CCNI) || (defined BGQ))
		std::cerr<<"Error: CCNI requires MPI."<<std::endl;
		std::exit(1);
		#endif
		tessellate<2,int>(grid, number_of_fields);
		#ifdef MPI_VERSION
		MPI::COMM_WORLD.Barrier();
		#endif
    if (rank==0) std::cout<<"Tessellation complete."<<std::endl;
	  #ifdef BGQ
	  output_bgq(grid, filename);
	  #else
	  output(grid, filename);
	  #endif
	} 
	else if (dim == 3) {
	  //	  const int edge = 512;
	  const int edge = 64;
	  int number_of_fields(seeds);
	  if (number_of_fields==0) number_of_fields = static_cast<int>(float(edge*edge*edge)/(4./3*M_PI*10.*10.*10.)); // Average grain is a sphere of radius 10 voxels
	  #ifdef MPI_VERSION
	  while (number_of_fields % np) --number_of_fields;
	  #endif
	  MMSP::grid<3,int > grid(0,0,edge,0,edge,0,edge);
	  if (rank==0) std::cout<<"Grid origin: ("<<g0(grid,0)<<','<<g0(grid,1)<<','<<g0(grid,2)<<"),"
				<<" dimensions: "<<g1(grid,0)-g0(grid,0)<<" × "<<g1(grid,1)-g0(grid,1)<<" × "<<g1(grid,2)-g0(grid,2)
				<<" with "<<number_of_fields<<" grains."<<std::endl;
	  #ifdef MPI_VERSION
	  number_of_fields /= np;
	  #endif

#if (!defined MPI_VERSION) && ((defined CCNI) || (defined BGQ))
	  std::cerr<<"Error: CCNI requires MPI."<<std::endl;
	  std::exit(1);
	  #endif
	  tessellate<3,int >(grid, number_of_fields);
	  #ifdef MPI_VERSION
	  MPI::COMM_WORLD.Barrier();
	  #endif
	  if (rank==0) std::cout<<"Tessellation complete."<<std::endl;
	  #ifdef BGQ
	  output_bgq(grid, filename);
	  #else
	  output(grid, filename);
	  #endif

	}
}


  template <int dim> struct flip_index{
    MMSP::grid<dim, int>* grid;
    vector<int> front_low_left_corner;
    vector<int> back_up_right_corner;
    int sublattice;
  };


void * flip_index_helper( void * s ){

  flip_index<3> *ss = static_cast<flip_index<3>*>(s);
  int sublattice = ss->sublattice;
  vector<int> x;
	const double kT = 0.50;
  // choose a random node 
  if(sublattice==0 || sublattice==2 || sublattice==4 || sublattice==6)
    x[0] = ss->front_low_left_corner[0]+2*rand()%((ss->back_up_right_corner[0]-ss->front_low_left_corner[0])/2);  //even
  else
    x[0] = ss->front_low_left_corner[0]+2*rand()%((ss->back_up_right_corner[0]-ss->front_low_left_corner[0])/2-1)+1; //odd
    
  if(sublattice==0 || sublattice==1 || sublattice==4 || sublattice==5)
    x[1] = ss->front_low_left_corner[1]+2*rand()%((ss->back_up_right_corner[1]-ss->front_low_left_corner[1])/2);  //even
  else
    x[1] = ss->front_low_left_corner[1]+2*rand()%((ss->back_up_right_corner[1]-ss->front_low_left_corner[1])/2-1)+1; //odd

  if(sublattice==0 || sublattice==1 || sublattice==2 || sublattice==3)
    x[2] = ss->front_low_left_corner[2]+2*rand()%((ss->back_up_right_corner[2]-ss->front_low_left_corner[2])/2);  //even
  else
    x[2] = ss->front_low_left_corner[2]+2*rand()%((ss->back_up_right_corner[2]-ss->front_low_left_corner[2])/2-1)+1; //odd


  int spin1 = (*(ss->grid))(x);
			
	// determine neighboring spins
	vector<int> r(x);
	sparse<bool> neighbors;
	for (int i=-1; i<=1; i++){
		for (int j=-1; j<=1; j++){
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
		for (int i=-1; i<=1; i++)
			for (int j=-1; j<=1; j++) {
			  r[0] = x[0] + i;
			  r[1] = x[1] + j;
			  int spin = (*(ss->grid))(r);
				dE += (spin!=spin2)-(spin!=spin1);
			}

	  // attempt a spin flip
		double r = double(rand())/double(RAND_MAX);
		if (dE<=0.0) (*(ss->grid))(x) = spin2;
		else if (r<exp(-dE/kT)) (*(ss->grid))(x) = spin2;
	}
  pthread_exit(0);
  return NULL;
}


template <int dim> void update(MMSP::grid<dim, int>& grid, int steps) {
#if (!defined MPI_VERSION) && ((defined CCNI) || (defined BGQ))
    std::cerr<<"Error: MPI is required for CCNI."<<std::endl;
    exit(1);
    #endif
    int rank=0;
    #ifdef MPI_VERSION
    rank=MPI::COMM_WORLD.Get_rank();
    #endif

    /*------------------*/
    if(dim == 2)
      std::cout<<"2d case needs to be done"<<std::endl;
    std::exit(1);

  int nthreads = 4;
	pthread_t * p_threads = new pthread_t[nthreads];
	flip_index<dim> * mat_para =	new flip_index<dim> [nthreads];
	pthread_attr_t attr;
	pthread_attr_init (&attr);

  int num_of_grids_on_edge = x1(grid, 0) - x0(grid, 0); // local grid edge length 
  vector<int> front_low_left_corner(3,0);
  vector<int> back_up_right_corner(3,0);

  for(int sublattice=0; sublattice!= 8; sublattice++){
	  for(int i=0; i!= nthreads ; i++ ) {

      front_low_left_corner[0] = x0(grid, 0) + (x1(grid, 0)-x0(grid, 0))/nthreads*i;
      front_low_left_corner[1] = x0(grid, 1);
      front_low_left_corner[2] = x0(grid, 2);
      back_up_right_corner[0] = x1(grid, 0) + (x1(grid, 0)-x0(grid, 0))/nthreads*(i+1)-1;
      back_up_right_corner[1] = x1(grid, 1);
      back_up_right_corner[2] = x1(grid, 2);

      mat_para[i].grid = &grid;
      mat_para[i].front_low_left_corner = front_low_left_corner;
      mat_para[i].back_up_right_corner = back_up_right_corner;
	  	mat_para[i].sublattice= sublattice;
	  	pthread_create(&p_threads[i], &attr, flip_index_helper, (void *) &mat_para[i] );
	  }

  	for(int i=0; i!= nthreads ; i++) {
  		pthread_join(p_threads[i], NULL);
	  }
    ghostswap(grid); // once loopd over a "color", ghostswap.
  }

	delete [] p_threads ;
  delete [] mat_para ;
}


template <class T> std::ostream& operator<<(std::ostream& o, int& s) {
  o<<"    IndexValue\n";
  for (int i=0; i<length(s); ++i) {
    int index = MMSP::index(s, i);
    o<<"    "<<std::setw(5)<<std::right<<index<<"  "<<s[index]<<'\n';
  }
  return o;
}

} 
#endif

#include"MMSP.main.hpp"
