#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<mpi.h>

int main(int argc, char ** argv)
{
	MPI_Init(&argc, &argv);
	MPI_Comm world = MPI_COMM_WORLD;
	int myRank, worldSize;
	MPI_Comm_size(world, &worldSize);
	MPI_Comm_rank(world, &myRank);

	int iterations = 1999999999;
	
	int inside = 0;
	int outside = 0;
	int TotalInside = 0;
	int TotalOutside = 0;
	
	srand(time(0) + myRank);
	for (int i = myRank; i < iterations; i+=worldSize)
	{
		double x = (double)rand() / (double)RAND_MAX;
		double y = (double)rand() / (double)RAND_MAX;
		if ((x*x + y*y) <= 1)
			inside++;
		outside++;
	}
	
	MPI_Reduce(&inside, &TotalInside, 1, MPI_INT, MPI_SUM, 0, world);
	MPI_Reduce(&outside, &TotalOutside, 1, MPI_INT, MPI_SUM, 0, world);
	
	if (myRank == 0)
	{
		double pi = 4*((double)TotalInside)/((double)TotalOutside);
		printf("Inside/Outside = %f\n", pi);
	}

	MPI_Finalize();
	return 0;
}
