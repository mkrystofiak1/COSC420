//Author: Mitchell Krystofiak
//Class: COSC420 - Lab1
//Date: September 10, 2021
//Description: Brute force primality testing. If composite, a factored output is provided.

#include<stdio.h>
#include<mpi.h>
#include<math.h>
#include<stdbool.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>

int main(int argc, char ** argv){

	MPI_Init(&argc, &argv);
	unsigned long long n;
	if (argc != 2)
	{
		srand(time(0));
		n = rand() % 1000;
		printf("No input to test! Using random num: %llu\n", n);	
	}
	else 
	{
		n = strtoull(argv[1], NULL, 10);
	}
	
	MPI_Comm world = MPI_COMM_WORLD;
	int worldSize;
	int myRank;

	MPI_Comm_size(world, &worldSize);
	MPI_Comm_rank(world, &myRank);
		
	//n = strtoull(argv[1], NULL, 10);
	
	bool result = true;

	if (n <= 1)
                result = false;

        unsigned long long i;
        double nn = n;
        double sqrtn = sqrt(nn);

        for (i = 2 + myRank; i <= (int)(sqrtn); i+=worldSize)
        {
		printf("Rank %d of %d is testing index %llu\n", myRank, worldSize, i);
                if ((n % i) == 0)
                {
			printf("Rank %d found Factorization: %llu = %llu * %llu\n", myRank, n, i, n/i);
                        result = false;
			break;
                }
        }
	
		
	printf("Rank %d of %d returned %d\n\n", myRank, worldSize, result);

	MPI_Finalize();
	return 0;
}
