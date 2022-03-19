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
		printf("\nBeginning Eigenvector test cases!\n");
		printf("----------------------------------------------\n\n");
	}

	int testcases[11] = {50, 100, 250, 400, 800, 1400, 2000, 3000, 4000, 8000, 10000};
	int testworks[4] = {2, 5, 8, 10};
	int tests = 11;
	int tests2 = 4;
	double t1, t2, dif;
	
	matrix A, B, C;
	
	if (myRank == 0)
	{
		printf("Testing some 'Make Sure It Works' cases:\n\n");
	}	
	
	/* Bevins Question
	if (myRank == 0)
	{
		zeroMatrix(&A, 2, 2);
		A.data[0] = 0;
		A.data[1] = 1;
		A.data[2] = -2;
		A.data[3] = -3;
		zeroMatrix(&B,2,1);
	}
	else
	{
		initNullMatrix(&A, 2, 2);
		initNullMatrix(&B, 2, 1);
	}
	
	writeToFile(&A, 2, 2, world, worldSize, myRank);
	B.data = EigenVector(2, 2, world, worldSize, myRank);
	
	if (myRank == 0)
	{
		printMatrix(&A);
		printMatrix(&B);
	}
	*/	

	for (int i = 0; i < tests2; i++)
	{
		if (myRank == 0)
		{
			randMatrix(&A, testworks[i], testworks[i], 1);
			randMatrix(&B, testworks[i], 1, 1);
			zeroMatrix(&C, testworks[i], 1);
		}
		else
		{
			initNullMatrix(&A, testworks[i], testworks[i]);
			initNullMatrix(&B, testworks[i], 1);
			initNullMatrix(&C, testworks[i], 1);
		}

		MPI_Barrier(world);
		t1 = MPI_Wtime();
		double normofB = norm(&B, world, worldSize, myRank);			
		MPI_Barrier(world);
		t2 = MPI_Wtime();
		MPI_Barrier(world);
		dif = t2 - t1;

		if(myRank == 0)
		{
			printf("Matrix A before Eigen Vector:\n\n");
			printMatrix(&A);
		}	
		writeToFile(&A, testworks[i], testworks[i], world, worldSize, myRank);
		
		MPI_Barrier(world);
		t1 = MPI_Wtime();
		C.data = EigenVector(testworks[i], testworks[i], world, worldSize, myRank);
		MPI_Barrier(world);
		t2 = MPI_Wtime();
		MPI_Barrier(world);

		if (myRank == 0)
		{
			printf("Eigen Vector of A:\n\n");
			printMatrix(&C);
			printf("Matrix Size: A = %dx%d, time taken to find Eigen Vector: %f seconds.\n\n", A.rows, A.cols, t2-t1);

			printf("Matrix B:\n\n");
			printMatrix(&B);
			printf("Matrix Size: B = %dx%d, time taken to find norm = %f: %f seconds.\n", B.rows, B.cols, normofB, dif);
		}
		//TODO ADD IN RETURN EIGENVALUE
		if (myRank == 0)
		{
			free(A.data);
			free(B.data);
			free(C.data);
		}
	}
	
	for (int i = 0; i < tests; i++)
	{
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
		double normofB = norm(&B, world, worldSize, myRank);			
		MPI_Barrier(world);
		t2 = MPI_Wtime();
		MPI_Barrier(world);
		dif = t2-t1;

		writeToFile(&A, testcases[i], testcases[i], world, worldSize, myRank);
		
		MPI_Barrier(world);
		t1 = MPI_Wtime();
		C.data = EigenVector(testcases[i], testcases[i], world, worldSize, myRank);
		MPI_Barrier(world);
		t2 = MPI_Wtime();
		MPI_Barrier(world);

		if (myRank == 0)
		{
			printf("Matrix Size: A = %dx%d, time taken to find Eigen Vector: %f seconds.\n\n", A.rows, A.cols, t2-t1);
			printf("Matrix Size: B = %dx%d, time taken to find norm of B: %f seconds.\n", B.rows, B.cols, dif);
		}

		if (myRank == 0)
		{
			free(A.data);
			free(B.data);
			free(C.data);
		}
	}
	MPI_Finalize();
	return 0;


}
