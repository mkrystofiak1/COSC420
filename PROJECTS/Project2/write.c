//Author: Mitchell Krystofiak
//Class: COSC 420 - Project 2 - Arxiv Search Engine
//Date: November 21, 2021
//Description: This program is meant to determine maximum string lengths
//             for storing purposes. It will go through the file, find the
//             largest lengths, and then write the important data to a file.

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<unistd.h>
#include<mpi.h>

int main(int argc, char ** argv)
{
	MPI_Init(&argc, &argv);
	MPI_Comm world = MPI_COMM_WORLD;
	int worldSize, myRank;
	MPI_Comm_size(world, &worldSize);
	MPI_Comm_rank(world, &myRank);


	FILE *fh;
	if ((fh = fopen("./arxiv-metadata.txt","r")) < 0)
	{
		printf("Error: File open error - arxiv-metadata.txt\n");
		MPI_Finalize();
		return 1;
	}
	
	char *buf = malloc(sizeof(char)*200000);
	int maxes[] = {0,0,0,0};
	int lineCount = 0;
	int count = 1;
	int numPapers = 0;
	int length = 0;
	
	//Not used yet
	int linebytes = 0;
	size_t totalbytes = 0;

	//Testing
	int updates[] = {0,0,0,0};

	if (myRank == 0)
	{	
		while (lineCount < 8140940)
		{
			memset(buf, '\0', 200000);

			fscanf(fh, "%[^\n]\n", buf);
			//to count the number of bytes to write to file
			//might be good to parallelize this so each node writes x files
			linebytes = sizeof(char)*length;
			totalbytes += linebytes;
			if (count == 1) //id
			{	
				if (length > maxes[0])
				{
					maxes[0] = length;
					updates[0]++;
				}
			}
			else if (count == 2) //title
			{
				if (length > maxes[1])
				{
					maxes[1] = length;
					updates[1]++;
				}
			}
			else if (count == 3) //author
			{	
				if (length > maxes[2])
				{
					maxes[2] = length;
					updates[2]++;
				}
			}
			else if (count == 4) //abstract
			{
				if (length > maxes[3])
				{
					maxes[3] = length;
					updates[3]++;
				}
			}
			else //plusline
			{
				count = 0;	
				numPapers++;
			}
			
			count++;
			lineCount++;

			//printf("%s\n", buf);
		}

		printf("Max id length: %d\n", maxes[0]);
		printf("Max title length: %d\n", maxes[1]);
		printf("Max author length: %d\n", maxes[2]);
		printf("Max abstract length: %d\n", maxes[3]);
		printf("Total papers: %d\n", numPapers);
		printf("Updates (id,ti,au,ab): (%d, %d, %d, %d)\n", updates[0], updates[1], updates[2], updates[3]);
		printf("Total bytes in file: %zu\n", totalbytes);
	
		/*
		fclose(fh);
		if ((fh = fopen("./arxiv-citations.txt", "r")) < 0)
		{
			printf("Error: File open error - arxiv-citations.txt!\n");
			return 1;
		}
		lineCount = 0;
		maxes[0] = 0;
		while (lineCount < 10913892)
		{
			memset(buf, '\0', 200000);

			fscanf(fh, "%[^\n]\n", buf);
			length = strlen(buf);
			if (length > maxes[0])
			{
				maxes[0] = length;
			}
		}
		printf("Max line length in Citations: %d\n", maxes[0]);
		*/
	}

	fclose(fh);
	free(buf);
	MPI_Finalize();
	return 0;
}

//int parseFile(
