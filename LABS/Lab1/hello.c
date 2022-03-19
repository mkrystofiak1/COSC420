//Author: Mitchell Krystofiak
//Class: COSC420 - Lab1
//Date: September 10, 2021
//Description: A simple hello world program to test the HPCL cluster and print out
//	       the processors name and the nodes rank and worldsize.

#include<stdio.h>
#include<mpi.h>

int main(int argc, char ** argv){
	
	MPI_Init(&argc, &argv);

	char name[MPI_MAX_PROCESSOR_NAME]; //Maximum length of string return by MPI_Get_Processor_Name
	int nameLen;

	MPI_Comm world = MPI_COMM_WORLD;
	int worldSize;
	int myRank;

	MPI_Comm_size(world, &worldSize); //Determines size of the group associated with communicator.
	MPI_Comm_rank(world, &myRank); //Determines the rank of the calling process in the communicator.
	
	MPI_Get_processor_name(name, &nameLen); //Gets the name of the processor.

	printf("Hello, world! from rank %d out of %d. My processor name is %s.\n", myRank, worldSize, name);
	MPI_Finalize();
	return 0;
}
