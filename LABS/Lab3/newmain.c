//Author: Mitchell Krystofiak
//Class: COSC 420 - Lab4
//Date: October 31, 2021 - November 7, 2021
//Description: Testing Gauss Jordan method.

#include"matrix.h"

/*
 * Write Matrix Data to a file to read in later.
 */
void writeToFile(matrix *A, int rows, int cols, MPI_Comm world, int worldSize, int myRank)
{
	MPI_File fh;

	int* sendCounts = malloc(sizeof(int)*worldSize);
	int* dsplc = malloc(sizeof(int)*worldSize);

	int tmp = 0;
	for (int i = 0; i < worldSize; i++)
	{
		sendCounts[i] = (rows*cols)/worldSize;
		dsplc[i] = tmp;
		tmp += sendCounts[i];
	}

	if ((rows*cols) % worldSize > 0) 
	{
		sendCounts[worldSize-1] += (rows*cols) % worldSize;
	}

	double* Atemp = malloc(sizeof(double)* sendCounts[myRank]);

	MPI_Scatterv(A->data, sendCounts, dsplc, MPI_DOUBLE, Atemp, sendCounts[myRank], MPI_DOUBLE, 0, world);

	MPI_Offset offset = myRank * sizeof(double) * sendCounts[myRank];

	MPI_File_open(world, "Matrix.data",
		MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

	MPI_File_write_at(fh, offset, Atemp, sendCounts[myRank], MPI_DOUBLE, MPI_STATUS_IGNORE);
	MPI_File_close(&fh);
	free(sendCounts);
	free(dsplc);
	free(Atemp);
}

/*
 * Read in a file to work with.
 */
double* readFromFile(int rows, int cols, MPI_Comm world, int worldSize, int myRank)
{
	MPI_File fh;
	matrix read;
	
	if (myRank == 0)
	{
		zeroMatrix(&read, rows, cols);
	}
	else
	{
		initNullMatrix(&read, rows, cols);
	}

	int* sendCounts = malloc(sizeof(int)*worldSize);
	int* dsplc = malloc(sizeof(int)*worldSize);

	int tmp = 0;
	for (int i = 0; i < worldSize; i++)
	{
		sendCounts[i] = (rows*cols)/worldSize;
		dsplc[i] = tmp;
		tmp += (rows*cols)/worldSize;
	}

	if ((rows*cols) % worldSize > 0)
	{
		sendCounts[worldSize-1] += (rows*cols) % worldSize;
	}

	double* Atemp = malloc(sizeof(double)*sendCounts[myRank]);

	MPI_Offset offset = myRank * sizeof(double) * sendCounts[myRank];

	MPI_File_open(world, "Matrix.data",
		MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);

	MPI_File_read_at(fh, offset, Atemp, sendCounts[myRank], MPI_DOUBLE, MPI_STATUS_IGNORE);
	MPI_File_close(&fh);

	MPI_Gatherv(Atemp, sendCounts[myRank], MPI_DOUBLE, read.data, sendCounts, dsplc, MPI_DOUBLE, 0, world);

	free(sendCounts);
	free(dsplc);
	free(Atemp);

	return read.data;
}

int main(int argc, char ** argv)
{
	MPI_Init(&argc, &argv);
	MPI_Comm world = MPI_COMM_WORLD;
	int worldSize, myRank;

	MPI_Comm_size(world, &worldSize);
	MPI_Comm_rank(world, &myRank);

	if (myRank == 0)
	{
		printf("\nBeginning Gauss Jordan Elimination test cases!\n");
		printf("----------------------------------------------\n\n");
	}

	int testcases[11] = {50, 100, 250, 400, 800, 1400, 2000, 3000, 4000, 8000, 10000};
	//int xtracases[] = {15000, 30000, 50000, 75000, 100000};
	int testworks[4] = {2, 5, 8, 10};
	int tests = 11;
	int tests2 = 4;
	double t1, t2;

	matrix A, B, C, D;
	
	int file = 1;

	if (file)
	{
		matrix X, Y;
		if (myRank == 0)
		{
			printf("Testing Writing/Reading from a file!\n\n");
			randMatrix(&X, 10, 10, 1);
			zeroMatrix(&Y, 10, 10);
		}
		else
		{
			initNullMatrix(&X, 10, 10);
			initNullMatrix(&Y, 10, 10);
		}

		writeToFile(&X, 10, 10, world, worldSize, myRank);
		Y.data = readFromFile(10, 10, world, worldSize, myRank);

		if (myRank == 0)
		{
			printf("Matrix Y from file:\n\n");
			printMatrix(&Y);
		}
	}

	if (argv[1])
	{
		if (myRank == 0)
		{
			printf("Argv[1] is populated, demonstrating read in!\n\n");
		}
		//TODO
	}
	
	if (myRank == 0)
	{
		printf("Testing some 'Make Sure It Works' cases:\n\n");
	}	


	for (int i = 0; i < tests2; i++)
	{
		if (myRank == 0)
		{
			randMatrix(&A, testworks[i], testworks[i], 1);
			idenMatrix(&B, testworks[i], testworks[i]);
			zeroMatrix(&C, testworks[i], testworks[i]);
			zeroMatrix(&D, testworks[i], testworks[i]);
		}
		else
		{
			initNullMatrix(&A, testworks[i], testworks[i]);
			initNullMatrix(&B, testworks[i], testworks[i]);
			initNullMatrix(&C, testworks[i], testworks[i]);
			initNullMatrix(&D, testworks[i], testworks[i]);
		}

		MPI_Barrier(world);
		t1 = MPI_Wtime();
		C.data = Gauss_Jordan(&A, &B, &world, worldSize, myRank);
		MPI_Barrier(world);
		t2 = MPI_Wtime();
		MPI_Barrier(world);
		
			
		if (myRank == 0)
		{
			printf("Matrix A:\n\n");
			printMatrix(&A);
			printf("Matrix B:\n\n");
			printMatrix(&B);
			printf("x Matrix C (post Gauss):\n\n");
			printMatrix(&C);
		}
		D.data = multiply(&C, &A, world, worldSize, myRank);

		if (myRank == 0)
		{
			printf("A*C = D:\n\n");
			printMatrix(&D);

			free(A.data);
			free(B.data);
			free(C.data);
			free(D.data);
		}
	}
	

	for (int i = 0; i < tests; i++)
	{
		
		/*
		 * Matrix solution
		 */
		if (myRank == 0)
		{
			randMatrix(&A, testcases[i], testcases[i], 1);
			randMatrix(&B, testcases[i], testcases[i], 1);
			zeroMatrix(&C, testcases[i], testcases[i]);
		}
		else
		{
			initNullMatrix(&A, testcases[i], testcases[i]);
			initNullMatrix(&B, testcases[i], testcases[i]);	
			initNullMatrix(&C, testcases[i], testcases[i]);
		}

		MPI_Barrier(world);
		t1 = MPI_Wtime();
		C.data = Gauss_Jordan(&A, &B, &world, worldSize, myRank);
		MPI_Barrier(world);
		t2 = MPI_Wtime();
		MPI_Barrier(world);
		if (myRank == 0)
		{
			free(A.data);
			free(B.data);
			free(C.data);
		}

		if (myRank == 0)
		{
			printf("Matrix Size: A = %dx%d, b = %dx%d, time taken: %f seconds.\n", testcases[i], testcases[i], testcases[i], testcases[i], t2-t1);
		}

		/*
		 * Matrix solution - Identity
		 */
		if (myRank == 0)
		{
			randMatrix(&A, testcases[i], testcases[i], 1);
			idenMatrix(&B, testcases[i], testcases[i]);
			zeroMatrix(&C, testcases[i], testcases[i]);
		}
		else
		{
			initNullMatrix(&A, testcases[i], testcases[i]);
			initNullMatrix(&B, testcases[i], testcases[i]);	
			initNullMatrix(&C, testcases[i], testcases[i]);
		}
		
		MPI_Barrier(world);
		t1 = MPI_Wtime();
		C.data = Gauss_Jordan(&A, &B, &world, worldSize, myRank);
		MPI_Barrier(world);
		t2 = MPI_Wtime();
		MPI_Barrier(world);
		if (myRank == 0)
		{
			free(A.data);
			free(B.data);
			free(C.data);
		}

		if (myRank == 0)
		{
			printf("Matrix Size: A = %dx%d, Identity b = %dx%d, time taken: %f seconds.\n", testcases[i], testcases[i], testcases[i], testcases[i], t2-t1);
		}
		
		/*
		 * vector solution
		 */
		if (myRank == 0)
		{
			randMatrix(&A, testcases[i], testcases[i], 1);
			randMatrix(&B, testcases[i], 1, 1);
			zeroMatrix(&C, testcases[i], 1);
		}
		else
		{
			initNullMatrix(&A, testcases[i], testcases[i]);
			initNullMatrix(&B, testcases[i], 1);
			initNullMatrix(&C, testcases[i], 1);
		}

		MPI_Barrier(world);
		t1 = MPI_Wtime();
		C.data = Gauss_Jordan(&A, &B, &world, worldSize, myRank);
		MPI_Barrier(world);
		t2 = MPI_Wtime();
		
		MPI_Barrier(world);
		if (myRank == 0)
		{
			free(A.data);
			free(B.data);
			free(C.data);
		}

		if (myRank == 0)
		{
			printf("Matrix Size: A = %dx%d, B = %dx%d, time taken: %f seconds.\n\n", testcases[i], testcases[i], testcases[i], 1, t2-t1);
		}
		
		MPI_Barrier(world);
	}
	MPI_Finalize();
	return 0;


}
