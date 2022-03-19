#include<stdio.h>
#include<mpi.h>

int main(int argc, char ** argv)
{
	MPI_Init(NULL,NULL);

	MPI_Comm world = MPI_COMM_WORLD;
	
	int worldSize, myRank;
	MPI_Comm_size(world, &worldSize);
	MPI_Comm_rank(world, &myRank);
	
	if (worldSize != 2)
	{
		fprintf(stderr, "Only run this with N=2 processors!\n");
		MPI_Finalize();
		return 1;
	}

	int token = 0;
	
	while( token < 10) 
	{
		if (myRank == 0)
		{
			token = 10;
			MPI_Send(&token, 1, MPI_INT,
				1, //destination rank
				0, //tag - used to associated different send/recv's
				world //comm
				);
		}
		else
		{
			MPI_Recv(&token, 1, MPI_INT,
				0, //source is rank 0
				0, //tag is 0 being sent
				world, //comm
				MPI_STATUS_IGNORE //don't care about status right now
						  //but later will store meta data about the transaction
				);
			printf("ranke %d received token %d\n", myRank, token);	
		}
	}
	printf("Rank %d done\n", myRank);
	MPI_Finalize();
	return 0;
}

