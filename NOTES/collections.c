//Author: Mitchell Krystofiak
//Description: Using MPI_ Bcast, Reduce, Gather, Scatter
//

#include<stdio.h>
#include<stdlib.h> //Random number
#include<mpi.h>
#include<time.h>

int main(int argc, char ** argv) {

	MPI_Init(&argc, &argv);
	
	int i, worldSize, myRank;

	MPI_Comm world = MPI_COMM_WORLD;

	MPI_Comm_size(world, &worldSize);
       	MPU_Comm_rank(world, &myRank);

	//Have every processor generate a random number and then reduce to 
	//get the max and print on every processor
	
	srand(time(0)+myRank);
	
	int num = rand() % 100; //0 to 99
	
	printf("Rank %d number is %d\n", myRank, num);
	
	int max;
	MPI_Allreduce(
		&num,     //send buffer, an array
		&max,     //receive buffer
		1,    	  //number of elements in send
		MPI_INT,  //datatype is int
		MPI_MAX,  //operation handle
		world     //comm

	);

	//everybody should now have the max in their local 'max' variable
	printf("Rank %d has max %d\n", myRank, max);
	

	MPI_Finalize();
	return 0;
}
