#include<stdio.h>
#include<mpi.h>

int main(int argc, char ** argv)
{
	MPI_Init(NULL,NULL);

	MPI_Comm world = MPI_COMM_WORLD;
	
	int worldSize, myRank;
	MPI_Comm_size(world, &worldSize);
	MPI_Comm_rank(world, &myRank);
	int limit = atoi[argv[1]];

	if (worldSize != 2)
	{
		fprintf(stderr, "Only run this with N=2 processors!\n");
		MPI_Finalize();
		return 1;
	}
	
	if (argc != 2)
	{
		fprintf(stderr, "Usage: ./ring [token limit]\n");
		MPI_Finalize();
		return -1;
	}

	int token = 0;
	
	while( token < limit) 
	{
		if (token % worldSize == rank) // % worldSize)
		{
			//sender block
			token++;
			printf("Rank %d sending token %d to %d\n", rank, token, (myRank+1)%worldSize);
			MPI_Ssend(&token, 1, MPI_INT,
				(myRank+1)%worldSize , //destination rank
				0, //tag - used to associated different send/recv's
				world //comm
				);
			if (token + worldSize -1 > limit) break;
		}
		else
		{
			//receiver block
			MPI_Recv(&token, 1, MPI_INT,
				rank == 0 ? worldSize-1: (rank-1), //source is rank 0
				0, //tag is 0 being sent
				world, //comm
				MPI_STATUS_IGNORE //don't care about status right now
						  //but later will store meta data about the transaction
				);
			printf("rank %d received token %d from %d\n", myRank, token, rank == 0 ? worldSize-1: (rank-1));
		}
	}
	printf("Rank %d done\n", myRank);
	MPI_Finalize();
	return 0;
}

