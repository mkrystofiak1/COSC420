#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<string.h>
#include<mpi.h>

//Generate some number and parallel write to a file:
//1. Without views
//2. With views


int main(int argc, char ** argv)
{
	MPI_Init(NULL, NULL);
	MPI_Comm world = MPI_COMM_WORLD;
	int rank, worldSize;
	MPI_Comm_size(world, &worldSize);
	MPI_Comm_rank(world, &rank);
	MPI_File fh;


	const char* fname = "outfile.data";
	int N = 10; //every proc write N ints to the file
	int offset;
	int nums[N]; //local array for now, need to malloc for large N
	
	srand(time(0)+rank);
	for (int i = 0; i < N; i++)
	{
		nums[i] = rand()%100;	
	}
	
	MPI_File_open(
		world, //communicator
		fname, //name of file
		MPI_MODE_CREATE | MPI_MODE_RDWR, //mode
	        MPI_INFO_NULL, //info structure
		&fh); //file handle

	offset= N*rank*sizeof(int); //where to start writing in the file by rank
	
	//to view the file as a human, need to decode the ascii
	//The 'hexdump' utility is made for this!
	// For n=10 we do:
	// NB: the 4 is because each int is 4 bytes
	// So 10/4 says to hexdump "for i = 1 to 10 read next 4 bytes and do format"
	// using %3d because we know each is at most two digits long, so 3 gives us
	// at least one whitespace in the output
	// hexdump -v -e '10/4 "%3d"' -e '"\n"' outfile.data
	//do window manually
	MPI_File_write_at(
		fh,
		offset,
		nums,
		N,
		MPI_INT,
		MPI_STATUS_IGNORE);

	MPI_File_close(&fh);


	MPI_Finalize();
	return 0;
}
