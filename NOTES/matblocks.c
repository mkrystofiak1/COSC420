#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<string.h>
#include<mpi.h>
#include"matrix.h"

//1. Generate a 10x10 matrix on the root
//2. Scatter it in 'windows' to four processors


int main(int argc, char ** argv)
{
	MPI_Init(NULL, NULL);
	MPI_Comm world = MPI_COMM_WORLD;
	int rank, worldSize;
	MPI_Comm_size(world, &worldSize);
	MPI_Comm_rank(world, &rank);

	matrix A;
	A.rows = A.cols = 10;
	A.data = NULL;
	int N = A.rows*A.cols;

	if (rank == 0)
	{
		A.data = malloc(N*sizeof(int));
		srand(time(0)+rank);
		for (int i = 0; i < N; i++)
		{
			A.data[i] = rand()%100;	
		}
		printMatrix(&A);
	}
	
	MPI_Datatype window;
	MPI_Type_vector(
		A.rows/2, //count -- number of contiguous blocks
		A.cols/2, //block length -- number of elements in each block
		A.cols,   //stride -- number of elements between block starts
		MPI_INT,  //oldtype
		&window  //newtype -- handle to the struct for modification
		);
	MPI_Type_commit(&window);
	
	// warmup: just send the first window to rank 1
	
	if (rank == 0)
	{
		MPI_Send(A.data, 1, window, 1, 0, world);
	}
	else
	{
		matrix Awin;
		Awin.rows = Awin.cols = A.rows/2;
		Awin.data = malloc(Awin.rows*Awin.cols*sizeof(int));
		MPI_Recv(Awin.data, Awin.rows*Awin.cols, MPI_INT, 0, 0, world, MPI_STATUS_IGNORE);
		printMatrix(&Awin);
	}	
		
	MPI_Type_free(&window);
	MPI_Finalize();
	return 0;
}
