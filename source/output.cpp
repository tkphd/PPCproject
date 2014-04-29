// output.cpp
// Custom output methods for IBM Blue Gene/Q supercomputers:
// specifically the RPI CCI BG/Q, "AMOS".

#ifdef MPI_VERSION

#include"MMSP.grid.hpp"

namespace MMSP {

template <int dim,typename T>
void output_bgq(const MMSP::grid<dim,T>& GRID, char* filename)
{
	/* MPI-IO to the filesystem with writes aligned to blocks */

	MPI::COMM_WORLD.Barrier();
	const unsigned int rank = MPI::COMM_WORLD.Get_rank();
	const unsigned int np = MPI::COMM_WORLD.Get_size();
	MPI_Request request;
	MPI_Status status;

	// Read filesystem block size (using statvfs). Default to 4096 B.
	struct statvfs buf;
	const unsigned long blocksize = (statvfs(".", &buf) == -1)?4096:buf.f_bsize;
	#ifdef DEBUG
	if (rank==0) std::cout<<"Block size is "<<blocksize<<" B."<<std::endl;
	#endif

	// file open error check
	//MPI::File output = MPI::File::Open(MPI::COMM_WORLD, filename, MPI::MODE_CREATE | MPI::MODE_WRONLY, MPI::INFO_NULL);
	MPI::INFO info;
	MPI_Info_create(&info);
	MPI_Info_set(&info, "IBM_largeblock_io", "true");
	MPI_File output;
	MPI_File_open(MPI::COMM_WORLD, filename, MPI::MODE_WRONLY|MPI::MODE_CREATE|MPI::MODE_EXCL, info, &output);
	if (!output) {
		if (rank==0) std::cerr << "File output error: could not open " << filename << "." << std::endl;
		exit(-1);
	}
	MPI_File_set_size(output, 0);

	// create buffer pointers
  unsigned long* datasizes = NULL;
  unsigned long* offsets = NULL;
  unsigned long* aoffsets = NULL;
  unsigned long* misalignments = NULL;
	char* databuffer=NULL;
  char* headbuffer=NULL;
  char* filebuffer=NULL;
  unsigned int* writeranks=NULL;
	MPI_Request* recvrequests = NULL;
	MPI_Status* recvstatuses = NULL;

	// get grid data to write
	const unsigned long size=write_buffer(GRID, databuffer);
	assert(databuffer!=NULL);
	// Generate MMSP header from rank 0
	unsigned long header_offset=0;
	if (rank==0) {
		std::stringstream outstr;
		// get grid data type
		std::string type = name(GRID);
		outstr << type << '\n';
		outstr << dim << '\n';
		outstr << MMSP::fields(GRID) << '\n';

		for (int i=0; i<dim; i++) outstr << MMSP::g0(GRID,i) << " " << MMSP::g1(GRID,i) << '\n'; // global grid dimensions
		for (int i=0; i<dim; i++) outstr << MMSP::dx(GRID,i) << '\n'; // grid spacing

		// Write MMSP header to file
		header_offset=outstr.str().size();
		headbuffer = new char[header_offset+sizeof(np)];
		memcpy(headbuffer, outstr.str().c_str(), header_offset);
		memcpy(headbuffer+header_offset, reinterpret_cast<const char*>(&np), sizeof(np));
		header_offset+=sizeof(np);
	}
	MPI::COMM_WORLD.Bcast(&header_offset, 1, MPI_UNSIGNED_LONG, 0); // broadcast header size from rank 0
  #ifdef DEBUG
  if (rank==0) std::cout<<"Prepared file header."<<std::endl;
  #endif
	MPI::COMM_WORLD.Barrier();

	// Compute file offsets based on buffer sizes
  datasizes = new unsigned long[np];
  MPI::COMM_WORLD.Allgather(&size, 1, MPI_UNSIGNED_LONG, datasizes, 1, MPI_UNSIGNED_LONG);
  #ifdef DEBUG
  if (rank==0) std::cout<<"Synchronized data sizes."<<std::endl;
  #endif

  // Determine disk space requirement
  unsigned long filesize=header_offset;
  for (unsigned int i=0; i<np; ++i) filesize+=datasizes[i];
	MPI::COMM_WORLD.Barrier();

	offsets = new unsigned long[np];
	offsets[0]=header_offset;
	for (unsigned int n=1; n<np; ++n) {
		assert(datasizes[n] < static_cast<unsigned long>(std::numeric_limits<int>::max()));
		offsets[n]=offsets[n-1]+datasizes[n-1];
	}
	offsets[0]=0;
	#ifdef DEBUG
	assert(datasizes[rank]==size);
	if (rank==0) std::cout<<"  Synchronized data offsets on "<<np<<" ranks. Total size: "<<offsets[np-1]+datasizes[np-1]<<" B."<<std::endl;
	#endif

	// Calculate number of  writers & write size
	unsigned long blocks = filesize/blocksize;
	while (blocks*blocksize<filesize)	++blocks;
	const unsigned int nwriters = (blocks>np)?np:blocks;
	const unsigned long writesize=blocksize*(blocks/nwriters);
	assert(writesize % blocksize==0);
	const unsigned long excessblocks=blocks % nwriters;
	bool isWriter=false;
	#ifdef DEBUG
	if (rank==0) std::cout<<"  Preparing "<<nwriters<<" aggregator/writers; writesize is "<<writesize<<" B, with "<<excessblocks<<" excess blocks."<<std::endl;
	#endif

	// Scan to determine which ranks are writers
	writeranks = new unsigned int[nwriters+1];
	aoffsets = new unsigned long[nwriters];
	writeranks[nwriters]=np-1; // generalization for last writer's benefit
	for (unsigned int w=0; w<nwriters; w++) {
		static unsigned int i=0;
		unsigned long ws=(w<=excessblocks)?writesize+blocksize:writesize;
		// file offset of the w^th writer
		aoffsets[w]=(w>0)?ws+aoffsets[w-1]:0;
		while ((aoffsets[w] > offsets[i]+datasizes[i]) && i<np)
			i++;
		writeranks[w]=i;
		if (rank==i)
			isWriter=true;
		i++;
	}

	// Determine which rank to send data to
	unsigned int prevwriter=nwriters, nextwriter=0;
	if (rank==0) {
		prevwriter=0;
	} else {
		while (writeranks[prevwriter]>=rank)
			--prevwriter;
	}
	if (rank>=writeranks[nwriters-1]) {
		nextwriter=nwriters;
	} else {
		while (writeranks[nextwriter]<=rank)
			++nextwriter;
	}

	unsigned long ws = writesize;
	if (nextwriter<=excessblocks)
		ws+=blocksize;
	if (rank>=writeranks[nwriters-1])
		ws=filesize-aoffsets[nwriters-1]; // last block may be only partially-filled

	/*
	#ifdef DEBUG
	if (rank==0)
		std::cout<<"Filesize is "<<filesize<<" B, or "<<blocks<<" blocks with "<<excessblocks<<" extra."<<std::endl;
	for (unsigned int r=0; r<np; r++) {
		MPI::COMM_WORLD.Barrier();
		if (rank==r) {
			for (unsigned int w=0; w<nwriters; w++)
				if (writeranks[w]==rank)
					printf("Rank %2u is a writer. Offset: %6lu B. Writesize: %6lu B. Datasize: %6lu B.\n", rank, aoffsets[w], ws, datasizes[rank]);
		}
		MPI::COMM_WORLD.Barrier();
	}
	if (rank==0) std::cout<<std::endl;
	MPI::COMM_WORLD.Barrier();
	#endif
	*/

	unsigned long deficiency=0;
	if (rank>0) {
		unsigned long prevws = (prevwriter>=excessblocks)?writesize:writesize+blocksize;
		deficiency = aoffsets[prevwriter]+prevws - offsets[rank];
		if (deficiency>datasizes[rank])
			deficiency=datasizes[rank];
	}
	// Collect block misalignments
   misalignments = new unsigned long[np];
   MPI::COMM_WORLD.Barrier();
   MPI::COMM_WORLD.Allgather(&deficiency, 1, MPI_UNSIGNED_LONG, misalignments, 1, MPI_UNSIGNED_LONG);

	#ifdef DEBUG
	if (datasizes[rank]-deficiency>ws)
		std::fprintf(stderr, "Error on Rank %u, alignment: buffered %lu B > writesize %lu B.\n", rank, datasizes[rank]-deficiency, ws);
	#endif
	/*
	#ifdef DEBUG
	for (unsigned int r=0; r<np; ++r) {
		MPI::COMM_WORLD.Barrier();
		if (r==rank) {
			std::fprintf(stderr, "Rank %2u: lower=%2u, defect=%6lu B, upper=%2u\n", rank, writeranks[prevwriter], misalignments[rank], writeranks[nextwriter]);
		}
		MPI::COMM_WORLD.Barrier();
	}
	if (rank==0) std::cout<<std::endl;
	#endif
	*/

	// Accumulate data
	const unsigned int silentranks=writeranks[nextwriter]-rank; // number of MPI ranks between this rank and the next writer
	MPI_Request sendrequest;
	MPI::COMM_WORLD.Barrier();
	if (isWriter) {
		// This rank is a writer.
		assert(misalignments[rank] < datasizes[rank]);
		// This rank is a writer
		#ifdef DEBUG
		if (rank>0 && writeranks[prevwriter+1]!=rank)
			std::fprintf(stderr, "Error on Rank %u, writer ID: %u != %u\n", rank, writeranks[prevwriter+1], rank);
		#endif

		// Copy local data into filebuffer
		filebuffer = new char[ws];
		char* p = filebuffer;
		if (rank==0) {
			memcpy(p, headbuffer, header_offset);
			p+=header_offset;
		}
		#ifdef DEBUG
		if (datasizes[rank]-misalignments[rank]>ws)
			std::fprintf(stderr, "Error on Rank %u, memcpy: %lu B > %lu B\n", rank, datasizes[rank]-misalignments[rank], ws);
		#endif
		char* q=databuffer+misalignments[rank];
		memcpy(p, q, datasizes[rank]-misalignments[rank]);
		p+=datasizes[rank]-misalignments[rank];

		// Recv remote data into filebuffer
		if (silentranks>0) {
			recvrequests = new MPI_Request[silentranks];
			recvstatuses = new MPI_Status[silentranks];
		}
		for (unsigned int i=0; i<silentranks && rank+i+1<np; i++) {
			unsigned int recv_proc = rank+i+1;
			assert(recv_proc!=rank && recv_proc<np);
			#ifdef DEBUG
			if (recv_proc<rank || recv_proc>np)
				std::fprintf(stderr, "Error on Rank %u, receiving: recv_proc=%i\n", rank, recv_proc);
			#endif
			unsigned long recv_size = misalignments[recv_proc];
			if (recv_size==0) continue;
			#ifdef DEBUG
			if (p+recv_size>filebuffer+ws)
				std::fprintf(stderr, "Error on Rank %u, receiving from %i: %lu B > %lu B\n", rank, recv_proc, p-filebuffer, ws-recv_size);
			#endif
			//recvrequests[i] = MPI::COMM_WORLD.Irecv(p, recv_size, MPI_CHAR, recv_proc, recv_proc);
			MPI_Irecv(p, recv_size, MPI_CHAR, recv_proc, recv_proc, MPI::COMM_WORLD, &recvrequests[i]);
			p+=recv_size;
		}
		#ifdef DEBUG
		if (p-filebuffer!=int(ws))
			std::fprintf(stderr, "Error on Rank %u, total received: %i B != %lu B\n", rank, int(p-filebuffer), ws);
		#endif
		if (rank>0 && misalignments[rank]>0) {
			q=databuffer;
			assert(writeranks[prevwriter]<rank);
			//sendrequest = MPI::COMM_WORLD.Isend(q, misalignments[rank], MPI_CHAR, writeranks[prevwriter], rank);
			MPI_Isend(q, misalignments[rank], MPI_CHAR, writeranks[prevwriter], rank, MPI::COMM_WORLD, &sendrequest);
		}
	}
	//MPI::COMM_WORLD.Barrier();
	if (misalignments[rank] >= datasizes[rank]) {
		assert(writeranks[prevwriter]<rank && writeranks[prevwriter]<np);
		//sendrequest = MPI::COMM_WORLD.Isend(databuffer, datasizes[rank], MPI_CHAR, writeranks[prevwriter], rank);
		MPI_Isend(databuffer, datasizes[rank], MPI_CHAR, writeranks[prevwriter], rank, MPI::COMM_WORLD, &sendrequest);
	}
	if (recvrequests != NULL)
		MPI_Waitall(silentranks, recvrequests, recvstatuses);
		//MPI::Request::Waitall(silentranks, recvrequests);
	if (rank>0) MPI_Wait(&sendrequest, &status);
	MPI::COMM_WORLD.Barrier();

	// Write to disk
	if (filebuffer!=NULL) {
		unsigned int w=0;
		while (writeranks[w]!=rank) ++w;
		assert(w<nwriters);
		if (w==nwriters-1)
			assert(filesize-aoffsets[w]==ws);
		//output.Write_at(aoffsets[w], filebuffer, ws, MPI_CHAR);
		MPI_File_iwrite_at(output, aoffsets[w], filebuffer, ws, MPI_CHAR, &request);
		MPI_Wait(&request, &status);
	}

	MPI::COMM_WORLD.Barrier();
	//output.Close();
	MPI_File_close(&output);
	if (recvrequests!=NULL) {
		delete [] recvrequests; recvrequests=NULL;
		delete [] recvstatuses; recvstatuses=NULL;
	}
	delete [] misalignments; misalignments = NULL;
	delete [] writeranks; writeranks=NULL;
	delete [] offsets; offsets=NULL;
	delete [] aoffsets; aoffsets=NULL;
	delete [] datasizes; datasizes=NULL;
	delete [] databuffer;	databuffer=NULL;
	if (filebuffer!=NULL) {
		delete [] filebuffer; filebuffer=NULL;
	}
}
} // namespace MMSP
#endif
