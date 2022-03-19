#include<stdio.h>
#include<stdlib.h>
#include<mpi.h>
#include<time.h>
#include<string.h>

//want a subroutine to print an array in a buffer manner for thread safety
char* printArr(int* arr, int n, int rank) 
{
	//build a string by appending the ints
	char* buf = malloc(n*(4+1)+1);
	buf[0] = '\0';
	const char* fnt = "%4d ";
	char tmp[4+1+1];
	for(int i = 0; i < n; i++)
	{
		sprintf(tmp, fnt, arr[i]);
		strcat(buf, tmp);	
	}
	printf("Rank %d numbers: %s\n", rank, buf );

	//every malloc needs a free!
	free(buf);
}

int main(int argc, char ** argv) {

	MPI_Init(&argc, &argv);
	
	int worldSize, rank;
	MPI_Comm world = MPI_COMM_WORLD;
	MPI_Comm_size(world, &worldSize);
	MPI_Comm_rank(world, &rank);

	/*
	 * Root node (rank = 0) make an array
	 * - everybody will need an array
	 * - only the root will have a big array
	 * Send 100 elements to each other node
	 * Those nodes will compute the max of their sub array
	 * Gather all the maxes back to the root and then find the final max
	*/
	int N = 10;
	int* nums = NULL;
       	int* sendNums = NULL; //working array, different on each node
	nums = malloc(N*size(int));

	if (rank == 0) //root
	{
		srand(time(0));
		//allocate space for 100 numbers per node
		//malloc needs the size in bytes, sizeof(int) gives the size of an int in bytes
		sendNums = malloc(worldSize*N*sizeof(int));
		for (int i =0; i < N*worldSize; i++)
		{
			sendNums[i] = rand()%100;
		}
		//everybody report what they received
		printArr(sendNums, N*worldSize, rank);
	}
		//printf("Node %d calling Scatter...\n", rank);
		MPI_Scatter(
			sendNums,  //holds things we want to send
			N,         //how many do each process get
			MPI_INT,   //type being sent
			nums,
			N,         //everybodies receiving 100
			MPI_INT,   //type being received
			0,         //root rank
			world      //communication handler
		);
		
		printArr(nums, N, rank);
		
		//find maxes of every node
		int max = 0;
		for(int i = 0; i < N; i++)
		{
			if (nums[i] > max){
				max = nums[i];
			}
		}

		//now gather the maxes back to the root
		int* maxes = NULL;
		if (rank == 0){
			maxes = malloc(worldSize*sizeof(int));
		}

		MPI_Gather(
			&max,
			1,
			MPI_INT,
			maxes,
			1,
			MPI_INT,
			0,
			world
		);

		if (rank == 0) {
			printf("Maxes received:\n");
			printArr(maxes, worldSize, rank);
		}

		free(nums);

		if (rank == 0) {
			free(sendNums);
		}
		if(rank == 0) free(sendNumes);

	
		MPI_Finalize();
		return 0;
}
	
	



