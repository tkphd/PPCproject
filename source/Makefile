# Makefile
# GNU Makefile for Voronoi tessellation and sparse phase-field grain growth
# Questions/comments to kellet@rpi.edu (Trevor Keller)

# includes
incdir = include
algodir = algorithms

# IBM XL compiler
BG_XL = /bgsys/drivers/ppcfloor/comm/xl
BG_INC = -I$(BG_PATH)/include
BG_LIB = -L$(BG_PATH)/lib

# compilers/flags
compiler = g++ -O3 -Wall
pcompiler = mpic++ -O3 -Wall
flags = -I$(incdir) -I$(algodir) -I$(algodir)/topology

# RPI CCI AMOS compilers/flags
#qcompiler = mpic++ -g -qarch=qp -qtune=qp -qflag=w -qstrict -qreport
qcompiler = mpic++ -O5 -qarch=qp -qtune=qp -qflag=w -qstrict -qprefetch=aggressive -qsimd=auto -qhot=fastmath -qinline=level=10
#qflags = $(CFLAGS) $(BG_INC) $(BG_LIB) $(LDFLAGS) $(flags)

# ONLY uncomment the following if <module load experimental/zlib> FAILS.
qflags = $(BG_INC) $(BG_LIB) $(flags) -I/bgsys/apps/CCNI/zlib/zlib-1.2.7/include -L/bgsys/apps/CCNI/zlib/zlib-1.2.7/lib

# dependencies
core = $(incdir)/MMSP.utility.hpp \
       $(incdir)/MMSP.grid.hpp \
       $(incdir)/MMSP.sparse.hpp

# the program
graingrowth.out: main.cpp graingrowth.cpp tessellate.hpp $(core)
	$(compiler) -DPHASEFIELD $(flags) $< -o $@ -lz -pthread

parallel: main.cpp graingrowth.cpp tessellate.hpp $(core)
	$(pcompiler) -DPHASEFIELD $(flags) -include mpi.h $< -o parallel_GG.out -lz

bgqmc: main.cpp graingrowth_MC.cpp tessellate.hpp $(core)
	$(qcompiler) $(qflags) -DBGQ -DSILENT $< -o q_MC.out -lz

bgq: main.cpp graingrowth.cpp tessellate.hpp $(core)
	$(qcompiler) $(qflags) -DBGQ -DSILENT -DPHASEFIELD $< -o q_GG.out -lz

wrongendian: wrongendian.cpp
	$(compiler) $< -o $@.out -lz -pthread

mmsp2vtk: mmsp2vtk.cpp $(core)
	$(compiler) $(flags) $< -o $@ -lz

clean:
	rm -rf graingrowth.out parallel_GG.out q_GG.out q_MC.out wrongendian.out
