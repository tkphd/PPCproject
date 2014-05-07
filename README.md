PPCproject
==========

RPI CSCI-6360 Parallel Programming &amp; Computing 2014: Group Project

If you wish to contribute to this code, please fork the repository into your GitHub account before cloning.

This is a group project intended to implement POSIX threads and to improve MPI-IO performance in MMSP (http://github.com/mesoscale/mmsp).

The code depends on the Mesoscale Microstructure Simulation Project (MMSP) libraries [1].

The source directory contains code for grain growth simulations using a sparse phase-field model [2] and a Potts Monte Carlo model [3].

Before building this code, download the MMSP source code and create symbolic links to ```MMSP/include``` and ```MMSP/algorithms```.

E.g., if you store MMSP in ~/Downloads/mmsp, run
```
ln -s ~/Downloads/mmsp/include source/include
ln -s ~/Downloads/mmsp/algorithms source/algorithms
```

In ```source/```, the Make commands default to the phase-field model:
e.g., ```make``` or ```make parallel``` will build ```graingrowth.out``` or ```parallel_GG.out```, respectively, both of which depend on ```graingrowth.cpp```.

To compile the Monte Carlo version, append "mc" to the commands:
e.g., ```make mc``` or ```make parallelmc``` will build ```graingrowth.out``` or ```parallel_MC.out```, respectively, both of which depend on ```graingrowth_MC.cpp```.

To build this code on an IBM Blue Gene/Q supercomputer, make sure the IBM xl_r compiler and zlib library are available.
On the Rensselaer Polytechnic Institute BG/Q, AMOS, ```module load xl_r experimental/zlib```.
Then, ```make bgq``` or ```make bgqmc```.

The recommended mode for this program is using --nonstop, which will intialize and coarsen a simulation domain.
e.g., use 
```
mpirun -np 6 ./parallel_GG.out --nonstop 3 voronoi.000.dat 200 200 4
```
to initialize and save a startup file (```voronoi.000.dat```) and a checkpoint file after 200 steps (```voronoi.200.dat```), using six MPI ranks with 4 pthreads each.
Make sure you have hardware support for 24 threads before running this command!

Note that AMOS is big-endian. Most consumer PCs are little-endian.
If you have downloaded a binary data file generated using MMSP on AMOS, ```make wrongendian``` and run it on the file.
Afterwards, you can use standard MMSP utilities (e.g., ```mmsp2vti```) to export a file for visualization.


References:

1.  J. Gruber and T. Keller. The Mesoscale Microstructure Simulation Project. http://github.com/mesoscale/mmsp

2.  I. Steinbach and F. Pezzolla. "A Generalized Field Method for Multiphase Transformations Using Interface Fields."
    Physica D 134 (1999) 385-393. DOI: 10.1016/S0167-2789(99)00129-3

3.  S.A. Wright, S.J. Plimpton, T.P. Swiler, R.M. Fye, M.F. Young, and E.A. Holm. "Potts-model grain growth simulations: Parallel algorithms and applications."
    SAND Report (1997) 1925.
