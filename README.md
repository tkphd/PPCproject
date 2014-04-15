PPCproject
==========

RPI CSCI-6360 Parallel Programming &amp; Computing 2014: Group Project

If you wish to contribute to this code, please fork the repository into your GitHub account before cloning.

This is a group project intended to implement POSIX threads and to improve MPI-IO performance in MMSP (http://github.com/mesoscale/mmsp).

The source directory contains code for grain growth simulations using a sparse phase-field method [1].

The code depends on the Mesoscale Microstructure Simulation Project (MMSP) libraries [2].

Before building this code, download the MMSP source code and create symbolic links to MMSP/include and MMSP/algorithms.

E.g., if you store MMSP in ~/Downloads/mmsp, run

	ln -s ~/Downloads/mmsp/include source/include

	ln -s ~/Downloads/mmsp/algorithms source/algorithms

To build the code in source/, make or make parallel.

To build on a Blue Gene/Q in source/, module load xl_r experimental/zlib, then make bgq.
Note that the RPI Blue Gene/Q, AMOS, is configured big-endian. Most consumer PCs are little-endian.
If you have downloaded a binary data file generated using MMSP on AMOS, make wrongendian and run it on the file.
Afterwards, you can use standard MMSP utilities (e.g., mmsp2vti) to export a file for visualization.


References:

1.  I. Steinbach and F. Pezzolla. "A Generalized Field Method for Multiphase Transformations Using Interface Fields."
    Physica D 134 (1999) 385-393. DOI: 10.1016/S0167-2789(99)00129-3

2.  J. Gruber. The Mesoscale Microstructure Simulation Project. http://github.com/mesoscale/mmsp
