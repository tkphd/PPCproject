// graingrowth.cpp
// Voronoi tessellations for synthetic microstructure test
// Questions/comments to kellet@rpi.edu (Trevor Keller)
#if (defined CCNI) || (defined BGQ)
#include<mpi.h>
#endif
#include<cstdio>
#include"MMSP.hpp"

std::string PROGRAM = "graingrowth";
std::string MESSAGE = "Voronoi tessellation and isotropic grain growth code";

typedef MMSP::grid<2,MMSP::sparse<float> > GRID2D;
typedef MMSP::grid<3,MMSP::sparse<float> > GRID3D;

