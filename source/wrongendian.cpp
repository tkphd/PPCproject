// wrongendian.cpp
// Change the endianness of MMSP data files
// Questions/comments to trevor.keller@gmail.com (Trevor Keller)

#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include<zlib.h>
#include<cstring>
#include<cstdlib>
#include<cstdio>
#include<pthread.h>

typedef struct {
	std::ifstream ifile;
	std::ofstream* ofile;
	unsigned long offset;
	int block;
} swap_thread;

#if defined(TIMING) || defined(DEBUG)
// rdtsc: hardware cycle counter for high-precision timing
#if defined(__i386__)
static __inline__ unsigned long long rdtsc(void)
{
	unsigned long long int x;
	__asm__ volatile (".byte 0x0f, 0x31" : "=A" (x));
	return x;
}

#elif defined(__x86_64__)
static __inline__ unsigned long long rdtsc(void)
{
	unsigned hi, lo;
	__asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
	return ( (unsigned long long)lo)|( ((unsigned long long)hi)<<32 );
}

#elif defined(__powerpc__)
static __inline__ unsigned long long rdtsc(void)
{
	unsigned long long int result=0;
	unsigned long int upper, lower,tmp;
	__asm__ volatile(
		"0:                    \n"
		"\tmftbu  %0           \n"
		"\tmftb   %1           \n"
		"\tmftbu  %2           \n"
		"\tcmpw   %2,%0        \n"
		"\tbne    0b           \n"
	: "=r"(upper),"=r"(lower),"=r"(tmp)
  );
  result = upper;
  result = result<<32;
  result = result|lower;

  return(result);
}
#endif

// Initialize timers
unsigned long tsc=0, readtimer=0, writetimer=0, swaptimer=0, pswaptimer=0, zutimer=0, zctimer=0;
#endif

template <typename T>
void swap_endian(T& n)
{
	// http://stackoverflow.com/questions/105252/how-do-i-convert-between-big-endian-and-little-endian-values-in-c
	union {
		T u;
		unsigned char u8[sizeof(T)];
	}
	source, dest;
	source.u = n;
	for (size_t k=0; k<sizeof(T); k++)
		dest.u8[k] = source.u8[sizeof(T)-k-1];
	n = dest.u;
}

template <typename T>
void swap_buffer(T* buffer, const unsigned long& size)
{
	// Swap endianness of buffer in-place
	T value;
	for (T* p=buffer; p<buffer+size; p+=sizeof(value))
		swap_endian<T>(*p);
}

template <typename T>
void* swap_block_kernel(void* x);

pthread_mutex_t write_lock;

// Define grid variables globablly, to ensure pthreads have access
std::string type;
int dim;
int fields;
int x0[3] = {0, 0, 0};
int x1[3] = {0, 0, 0};
float dx[3] = {1.0, 1.0, 1.0};
int blocks;

bool scalar_type, vector_type, sparse_type;
bool bool_type, char_type, unsigned_char_type, int_type, unsigned_int_type, long_type, unsigned_long_type, short_type, unsigned_short_type, float_type, double_type, long_double_type;

int main(int argc, char* argv[])
{
	#if defined(TIMING) || defined(DEBUG)
	unsigned long exectimer=rdtsc();
	#endif
	// command line error check
	if (argc < 3) {
		std::cout<<"Usage: "<<argv[0]<<" [--help] infile outfile [threads]\n\n";
		exit(-1);
	}

	// help diagnostic
	if (std::string(argv[1]) == "--help" || argc<3) {
		std::cout<<argv[0]<<": Change the endianness of MMSP data files.\n";
		std::cout<<"Usage: "<<argv[0]<<" [--help] infile outfile\n\n";
		std::cout<<"Questions/comments to trevor.keller@gmail.com (Trevor Keller).\n\n";
		exit(0);
	}

	// file open error check
	std::ifstream input(argv[1]);
	if (!input) {
		std::cerr<<"File input error: could not open "<<argv[1]<<".\n"<<std::endl;
		exit(-1);
	}

	#if defined(TIMING) || defined(DEBUG)
	#ifdef BGQ
	double clock_rate=1600000000.0;
	#else
	// Read clock rate from GNU/Linux machine
	double clock_rate=2666700000.0;
	std::ifstream fh("/sys/devices/system/cpu/cpu0/cpufreq/scaling_available_frequencies");
	if (!fh)
		clock_rate=2667000000.0;
	else {
		fh>>clock_rate;
		fh.close();
	}
	clock_rate*=1000;
	#endif
	#endif

	// read data type
	#if defined(TIMING) || defined(DEBUG)
	tsc=rdtsc();
	#endif
	getline(input, type, '\n');
	#if defined(TIMING) || defined(DEBUG)
	readtimer+=rdtsc()-tsc;
	#endif

	// grid type error check
	if (type.substr(0, 4) != "grid") {
		std::cerr<<"File input error: file does not contain grid data.\n"<<std::endl;
		exit(-1);
	}

	// parse data type
	//scalar_type = (type.find("scalar") != std::string::npos);
	//vector_type = (type.find("vector") != std::string::npos);
	sparse_type = (type.find("sparse") != std::string::npos);
	bool_type = (type.find("bool") != std::string::npos);
	char_type = (type.find("char") != std::string::npos);
	unsigned_char_type = (type.find("unsigned char") != std::string::npos);
	int_type = (type.find("int") != std::string::npos);
	unsigned_int_type = (type.find("unsigned int") != std::string::npos);
	long_type = (type.find("long") != std::string::npos);
	unsigned_long_type = (type.find("unsigned long") != std::string::npos);
	short_type = (type.find("short") != std::string::npos);
	unsigned_short_type = (type.find("unsigned short") != std::string::npos);
	float_type = (type.find("float") != std::string::npos);
	double_type = (type.find("double") != std::string::npos);
	long_double_type = (type.find("long double") != std::string::npos);

	if (not bool_type    and
  		not char_type    and  not unsigned_char_type   and
			not int_type     and  not unsigned_int_type    and
			not long_type    and  not unsigned_long_type   and
			not short_type   and  not unsigned_short_type  and
			not float_type   and
			not double_type  and  not long_double_type)
	{
		std::cerr<<"File input error: unknown grid data type ("<<type<<").\n"<<std::endl;
		exit(-1);
	} else if (not sparse_type) {
		std::cerr<<"File input error: grid type ("<<type<<") is not implemented.\n"<<std::endl;
		exit(-1);
	}

	// file open error check
	std::ofstream output(argv[2]);
	if (!output) {
		std::cerr<<"File output error: could not open "<<argv[2]<<".\n"<< std::endl;
		exit(-1);
	}

	// Default to 2 threads, unless otherwise specified.
	const unsigned int nthreads=(argc==4)?atoi(argv[3]):2;
	if (nthreads<1) {
		std::cerr<<"POSIX thread error: "<<nthreads<<" threads is too few.\n"<< std::endl;
		exit(-1);
	}

	// write data type
	#if defined(TIMING) || defined(DEBUG)
	tsc=rdtsc();
	#endif
	output<<type<<'\n';
	#if defined(TIMING) || defined(DEBUG)
	writetimer+=rdtsc()-tsc;
	#endif

	// copy grid dimension
	#if defined(TIMING) || defined(DEBUG)
	tsc=rdtsc();
	#endif
	input>>dim;
	#if defined(TIMING) || defined(DEBUG)
	readtimer+=rdtsc()-tsc;
	tsc=rdtsc();
	#endif
	output<<dim<<'\n';
	#if defined(TIMING) || defined(DEBUG)
	writetimer+=rdtsc()-tsc;
	#endif

	// copy number of fields
	#if defined(TIMING) || defined(DEBUG)
	tsc=rdtsc();
	#endif
	input>>fields;
	#if defined(TIMING) || defined(DEBUG)
	readtimer+=rdtsc()-tsc;
	tsc=rdtsc();
	#endif
	output<<fields<<'\n';
	#if defined(TIMING) || defined(DEBUG)
	writetimer+=rdtsc()-tsc;
	#endif

	// copy grid sizes
	for (int i = 0; i < dim; i++){
		#if defined(TIMING) || defined(DEBUG)
		tsc=rdtsc();
		#endif
		input >> x0[i] >> x1[i];
		#if defined(TIMING) || defined(DEBUG)
		readtimer+=rdtsc()-tsc;
		tsc=rdtsc();
		#endif
		output<<x0[i]<<' '<<x1[i]<<'\n';
		#if defined(TIMING) || defined(DEBUG)
		writetimer+=rdtsc()-tsc;
		#endif
	}

	// copy cell spacing
	for (int i = 0; i < dim; i++){
		#if defined(TIMING) || defined(DEBUG)
		tsc=rdtsc();
		#endif
		input >> dx[i];
		#if defined(TIMING) || defined(DEBUG)
		readtimer+=rdtsc()-tsc;
		tsc=rdtsc();
		#endif
		output<<dx[i]<<'\n';
		#if defined(TIMING) || defined(DEBUG)
		writetimer+=rdtsc()-tsc;
		#endif
	}

	// ignore trailing endlines
	#if defined(TIMING) || defined(DEBUG)
	tsc=rdtsc();
	#endif
	input.ignore(10, '\n');
	#if defined(TIMING) || defined(DEBUG)
	readtimer+=rdtsc()-tsc;
	#endif


	// copy number of blocks
	#if defined(TIMING) || defined(DEBUG)
	tsc=rdtsc();
	#endif
	input.read(reinterpret_cast<char*>(&blocks), sizeof(blocks));
	#if defined(TIMING) || defined(DEBUG)
	readtimer+=rdtsc()-tsc;
	tsc=rdtsc();
	#endif
	swap_endian(blocks);
	#if defined(TIMING) || defined(DEBUG)
	swaptimer+=rdtsc()-tsc;
	tsc=rdtsc();
	#endif
	output.write(reinterpret_cast<const char*>(&blocks), sizeof(blocks));
	#if defined(TIMING) || defined(DEBUG)
	writetimer+=rdtsc()-tsc;
	#endif
	#ifdef DEBUG
	std::cout<<blocks<<" blocks"<<std::endl;
	#endif


	unsigned long pos=input.tellg();
	int b=0;
	while (b < blocks) {
		pthread_t* p_threads = new pthread_t[nthreads];
		swap_thread* swap_threads = new swap_thread[nthreads];
		pthread_attr_t attr;
		pthread_attr_init (&attr);
		int spawned=0;
		for (int i=0; i<nthreads; i++) {
			if (b<blocks) {
				//#ifdef DEBUG
				//printf("Block %d of %d: offset %lu.\n", b, blocks, pos);
				//#endif
				input.seekg(pos);
				swap_threads[i].block = b-1;
				swap_threads[i].ifile.open(argv[1]);
				swap_threads[i].offset = pos;
				swap_threads[i].ofile = &output;

				pthread_create(&p_threads[i], &attr, swap_block_kernel<float>, (void *) &swap_threads[i] );
				spawned++;

				b++;

				if (b<blocks) {
					pos+=4*dim*sizeof(int)+sizeof(unsigned long);
					unsigned long datasize;
					input.seekg(pos);
					input.read(reinterpret_cast<char*>(&datasize), sizeof(unsigned long));
					swap_endian(datasize);
					pos+=sizeof(unsigned long) + datasize*sizeof(Bytef);
				}
			}
		}
		// Wait for pthreads to exit
		for (int i=0; i<spawned; i++)
			pthread_join(p_threads[i], NULL);
		#ifdef DEBUG
		std::cout<<std::endl;
		#endif

		pthread_attr_destroy(&attr);
		delete [] p_threads;
		delete [] swap_threads;
	}
	#ifdef DEBUG
	std::cout<<"Finished loop."<<std::endl;
	#endif

	input.close();
	output.close();

	std::cout<<"Endianness of "<<argv[1]<<" successfully inverted."<<std::endl;

	#if defined(TIMING) || defined(DEBUG)
	printf("%2.2e sec. reading data.\n", readtimer/clock_rate);
	printf("%2.2e sec. uncompressing data.\n", zutimer/clock_rate);
	printf("%2.2e sec. swapping data.\n", swaptimer/clock_rate);
	printf("%2.2e sec. swapping data, ×%i threads.\n", pswaptimer/clock_rate, nthreads);
	printf("%2.2e sec. compressing data.\n", zctimer/clock_rate);
	printf("%2.2e sec. writing data.\n", writetimer/clock_rate);
	printf("\n%2.2e sec. total execution time.\n", (rdtsc()-exectimer)/clock_rate);
	#endif
}

template <typename T>
void* swap_block_kernel(void* x)
{
	swap_thread* st = static_cast<swap_thread*>(x);
	st->ifile.seekg(st->offset);
	// copy block limits
	int lmin[3] = {0, 0, 0};
	int lmax[3] = {0, 0, 0};
	for (int j = 0; j < dim; j++) {
		#if defined(TIMING) || defined(DEBUG)
		tsc=rdtsc();
		#endif
		st->ifile.read(reinterpret_cast<char*>(&lmin[j]), sizeof(lmin[j]));
		st->ifile.read(reinterpret_cast<char*>(&lmax[j]), sizeof(lmax[j]));
		#if defined(TIMING) || defined(DEBUG)
		readtimer+=rdtsc()-tsc;
		#endif
		#if defined(TIMING) || defined(DEBUG)
		tsc=rdtsc();
		#endif
		swap_endian<int>(lmin[j]);
		swap_endian<int>(lmax[j]);
		#if defined(TIMING) || defined(DEBUG)
		swaptimer+=rdtsc()-tsc;
		#endif
	}

	// copy boundary conditions
	int blo[dim];
	int bhi[dim];
	for (int j = 0; j < dim; j++) {
		#if defined(TIMING) || defined(DEBUG)
		tsc=rdtsc();
		#endif
		st->ifile.read(reinterpret_cast<char*>(&blo[j]), sizeof(blo[j]));
		st->ifile.read(reinterpret_cast<char*>(&bhi[j]), sizeof(bhi[j]));
		#if defined(TIMING) || defined(DEBUG)
		readtimer+=rdtsc()-tsc;
		#endif
		#if defined(TIMING) || defined(DEBUG)
		tsc=rdtsc();
		#endif
		swap_endian<int>(blo[j]);
		swap_endian<int>(bhi[j]);
		#if defined(TIMING) || defined(DEBUG)
		swaptimer+=rdtsc()-tsc;
		#endif
	}

	// copy grid data
	unsigned long size_in_mem, size_on_disk;
	#if defined(TIMING) || defined(DEBUG)
	tsc=rdtsc();
	#endif
	st->ifile.read(reinterpret_cast<char*>(&size_in_mem), sizeof(size_in_mem)); // read raw size
	st->ifile.read(reinterpret_cast<char*>(&size_on_disk), sizeof(size_on_disk)); // read compressed size
	#if defined(TIMING) || defined(DEBUG)
	readtimer+=rdtsc()-tsc;
	#endif
	#if defined(TIMING) || defined(DEBUG)
	tsc=rdtsc();
	#endif
	swap_endian<unsigned long>(size_in_mem);
	swap_endian<unsigned long>(size_on_disk);
	#if defined(TIMING) || defined(DEBUG)
	swaptimer+=rdtsc()-tsc;
	#endif
	#ifdef DEBUG
	printf("Block %d: Reading %lu B (%lu KB) into %lu B (%lu KB) buffer.\n", st->block+1, size_on_disk, size_on_disk/1024, size_in_mem, size_in_mem/1024);
	#endif

	// Invert buffered values
	Bytef* buffer = new Bytef[size_on_disk];
	#if defined(TIMING) || defined(DEBUG)
	tsc=rdtsc();
	#endif
	st->ifile.read(reinterpret_cast<char*>(buffer), size_on_disk);
	#if defined(TIMING) || defined(DEBUG)
	readtimer+=rdtsc()-tsc;
	#endif
	st->ifile.close();
	if (size_on_disk!=size_in_mem) {
		size_in_mem=1.125*size_in_mem+12;
		char* raw = new char[size_in_mem];
		// Uncompress data
		#if defined(TIMING) || defined(DEBUG)
		tsc=rdtsc();
		#endif
		int status = uncompress(reinterpret_cast<Bytef*>(raw), &size_in_mem, buffer, size_on_disk);
		#if defined(TIMING) || defined(DEBUG)
		zutimer+=rdtsc()-tsc;
		#endif
		if (status<0) std::cerr<<"\nzlib error: "<<status<<std::endl;
		switch( status ) {
			case Z_OK:
				break;
			case Z_MEM_ERROR:
				std::cerr << "Uncompress: out of memory.\n" << std::endl;
				st->ofile->close();
				delete [] raw;
				delete [] buffer;
				exit(-1);
				break;
			case Z_BUF_ERROR:
				std::cerr << "Uncompress: output buffer ("<<size_in_mem<<" B) wasn't large enough for data ("<<size_on_disk<<" B).\n" << std::endl;
				st->ofile->close();
				delete [] raw;
				delete [] buffer;
				exit(-1);
				break;
		}
		//#ifdef DEBUG
		//printf("Block %d: Decompressed %lu B (%lu KB).\n", st->block+1, size_in_mem, size_in_mem/1024);
		//#endif
		// Invert raw data
		char* p = raw;// + sizeof(blockedvals);
		if (sparse_type) {
			if (float_type) {
				for (char* q=p; q<raw+size_in_mem; q+=sizeof(int)) {
					// Swap the number of floats in each sparse object
					int nfloats=0;
					#if defined(TIMING) || defined(DEBUG)
					tsc=rdtsc();
					#endif
					swap_buffer<int>(reinterpret_cast<int*>(q), 1);
					#if defined(TIMING) || defined(DEBUG)
					swaptimer+=rdtsc()-tsc;
					#endif
					memcpy(&nfloats, reinterpret_cast<int*>(q), 1);
					q+=sizeof(int);
					for (int j=0; j<nfloats; j++) {
						// MMSP::sparse::to_buffer copies n * item<T>;
						// each item contains one int and one T
						swap_buffer<int>(reinterpret_cast<int*>(q), 1);
						q+=sizeof(int);
						swap_buffer<float>(reinterpret_cast<float*>(q), 1);
						q+=sizeof(float);
					}
				}
			} else {
				std::cerr<<"ERROR: Grid type ("<<type<<") is not implemented.\n"<<std::endl;
				delete [] raw;
				delete [] buffer;
				exit(-1);
			}
		} else {
			std::cerr<<"ERROR: Grid type ("<<type<<") is not implemented.\n"<<std::endl;
			delete [] raw;
			delete [] buffer;
			exit(-1);
		}
		delete [] buffer; buffer=NULL;
		// Re-compress
		size_on_disk=1.125*size_in_mem+12;
		buffer = new Bytef[size_on_disk];
		#if defined(TIMING) || defined(DEBUG)
		tsc=rdtsc();
		#endif
		status = compress2(buffer, &size_on_disk, reinterpret_cast<Bytef*>(raw), size_in_mem, 9);
		#if defined(TIMING) || defined(DEBUG)
		zctimer+=rdtsc()-tsc;
		#endif
		if (status<0) std::cerr<<"\n\nzlib error: "<<status<<'\n'<<std::endl;
		switch(status) {
			case Z_OK:
				break;
			case Z_MEM_ERROR:
				std::cerr << "Compress: out of memory.\n" << std::endl;
				st->ofile->close();
				delete [] buffer;
				exit(-1);
				break;
			case Z_BUF_ERROR:
				std::cerr << "Uncompress: output buffer ("<<size_on_disk<<" B) wasn't large enough for data ("<<size_in_mem<<" B).\n" << std::endl;
				st->ofile->close();
				delete [] buffer;
				exit(-1);
				break;
		}
		//#ifdef DEBUG
		//printf("Block %d: Compressed into %lu B (%lu KB).\n", st->block+1, size_on_disk, size_on_disk/1024);
		//#endif
		delete [] raw; raw=NULL;
	}	else {
		if (sparse_type) {
			if (float_type) {
				for (char* q=reinterpret_cast<char*>(buffer); q<reinterpret_cast<char*>(buffer)+size_in_mem; q+=sizeof(int)) {
					// Swap the number of floats in each sparse object
					int nfloats=0;
					#if defined(TIMING) || defined(DEBUG)
					tsc=rdtsc();
					#endif
					swap_buffer<int>(reinterpret_cast<int*>(q), 1);
					#if defined(TIMING) || defined(DEBUG)
					swaptimer+=rdtsc()-tsc;
					#endif
					memcpy(&nfloats, reinterpret_cast<int*>(q), 1);
					q+=sizeof(int);
					for (int j=0; j<nfloats; j++) {
						swap_buffer<float>(reinterpret_cast<float*>(q), 1);
						q+=sizeof(float);
					}
				}
				// Compress
				size_on_disk=1.25*size_in_mem+12;
				Bytef* raw=buffer;
				buffer = new Bytef[size_on_disk];
				#if defined(TIMING) || defined(DEBUG)
				tsc=rdtsc();
				#endif
				int status = compress2(buffer, &size_on_disk, raw, size_in_mem, 9);
				#if defined(TIMING) || defined(DEBUG)
				zctimer+=rdtsc()-tsc;
				#endif
				switch(status) {
					case Z_OK:
						break;
					case Z_MEM_ERROR:
						std::cerr << "Compress: out of memory.\n" << std::endl;
						delete [] buffer;
						exit(-1);
						break;
					case Z_BUF_ERROR:
						std::cerr << "Compress: output buffer wasn't large enough.\n" << std::endl;
						delete [] buffer;
						exit(-1);
						break;
				}
				//#ifdef DEBUG
				//printf("Block %d: Compressed %lu B into %lu B (%lu KB).\n", st->block+1, size_in_mem, size_on_disk, size_on_disk/1024);
				//#endif
				delete [] raw; raw=NULL;
			}
		} else {
			std::cerr<<"ERROR: Grid type ("<<type<<") is not implemented.\n"<<std::endl;
			exit(-1);
		}
	}

	#if defined(TIMING) || defined(DEBUG)
	tsc=rdtsc();
	#endif

	pthread_mutex_lock(&write_lock);
	#if defined(TIMING) || defined(DEBUG)
	tsc=rdtsc();
	#endif
	for (int j = 0; j < dim; j++) {
		st->ofile->write(reinterpret_cast<const char*>(&lmin[j]), sizeof(lmin[j]));
		st->ofile->write(reinterpret_cast<const char*>(&lmax[j]), sizeof(lmax[j]));
	}
	for (int j = 0; j < dim; j++) {
		st->ofile->write(reinterpret_cast<const char*>(&blo[j]), sizeof(blo[j]));
		st->ofile->write(reinterpret_cast<const char*>(&bhi[j]), sizeof(bhi[j]));
	}
	st->ofile->write(reinterpret_cast<const char*>(&size_in_mem), sizeof(size_in_mem));
	st->ofile->write(reinterpret_cast<const char*>(&size_on_disk), sizeof(size_on_disk));
	st->ofile->write(reinterpret_cast<const char*>(buffer), size_on_disk);
	#if defined(TIMING) || defined(DEBUG)
	writetimer+=rdtsc()-tsc;
	#endif
	pthread_mutex_unlock(&write_lock);

	delete [] buffer;

	//#ifdef DEBUG
	//printf("Block %d: Finished.\n", st->block+1);
	//#endif

	pthread_exit((void*)0);
	return NULL;
}
