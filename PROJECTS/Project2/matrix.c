//Program: matrix.c
//Author: Mitchell Krystofiak
//Description: Implementation of the matrix class
//Date: October 31, 2021

#include "matrix.h"

/*
 * Remember: the structure of matrix typdef is:
 * 	int rows, int cols, double * data
 */

/*
 * Write Matrix Data to a file to read in later.
 */
void writeToFile(matrix *A, int rows, int cols, MPI_Comm world, int worldSize, int myRank)
{
	MPI_File fh;

	int *sendCounts = malloc(sizeof(int) * worldSize);
	int *dsplc = malloc(sizeof(int) * worldSize);

	int tmp = 0;
	for (int i = 0; i < worldSize; i++)
	{
		sendCounts[i] = (rows * cols) / worldSize;
		dsplc[i] = tmp;
		tmp += sendCounts[i];
	}

	if ((rows * cols) % worldSize > 0)
	{
		sendCounts[worldSize - 1] += (rows * cols) % worldSize;
	}

	double *Atemp = malloc(sizeof(double) * sendCounts[myRank]);

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
double *readFromFile(int rows, int cols, MPI_Comm world, int worldSize, int myRank)
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

	int *sendCounts = malloc(sizeof(int) * worldSize);
	int *dsplc = malloc(sizeof(int) * worldSize);

	int tmp = 0;
	for (int i = 0; i < worldSize; i++)
	{
		sendCounts[i] = (rows * cols) / worldSize;
		dsplc[i] = tmp;
		tmp += (rows * cols) / worldSize;
	}

	if ((rows * cols) % worldSize > 0)
	{
		sendCounts[worldSize - 1] += (rows * cols) % worldSize;
	}

	double *Atemp = malloc(sizeof(double) * sendCounts[myRank]);

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

/*
 * Initializes a matrix with random values.
 NOTE I CHANGED THIS TO BE ABLE TO TEST PAGERANK CHANGED RAND CHANGE BACK ONCE
 */
void randMatrix(matrix *A, int rows, int cols, int c)
{
	srand(time(0));
	A->rows = rows;
	A->cols = cols;
	A->data = malloc(A->rows * A->cols * sizeof(double));

	if (c == 1) //int
	{
		for (int i = 0; i < A->rows; i++)
		{
			for (int j = 0; j < A->cols; j++)
			{
				ACCESS(A, i, j) = rand() % 2;
			}
		}
	}
	else if (c == 0) //double
	{
		for (int i = 0; i < A->rows; i++)
		{
			for (int j = 0; j < A->cols; j++)
			{
				ACCESS(A, i, j) = (((double)(rand() % 100)) / ((double)(rand() % 50) + 1));
			}
		}
	}
}

/*
 * Initializaes a matrix with zero values.
 */
void zeroMatrix(matrix *A, int rows, int cols)
{
	A->rows = rows;
	A->cols = cols;
	A->data = malloc(A->rows * A->cols * sizeof(double));
	for (int i = 0; i < A->rows; i++)
	{
		for (int j = 0; j < A->cols; j++)
		{
			ACCESS(A, i, j) = 0.0;
		}
	}
}

/*
 * Initializes a Identity Matrix.
 */
void idenMatrix(matrix *A, int rows, int cols)
{
	A->rows = rows;
	A->cols = cols;
	A->data = malloc(A->rows * A->cols * sizeof(double));
	for (int i = 0; i < A->rows; i++)
	{
		for (int j = 0; j < A->cols; j++)
		{
			if (i == j)
			{
				ACCESS(A, i, j) = 1.0;
			}
			else
			{
				ACCESS(A, i, j) = 0.0;
			}
		}
	}
}

/*
 * Prints the matrix out.
 */
void printMatrix(matrix *A)
{
	for (int i = 0; i < A->rows; i++)
	{
		for (int j = 0; j < A->cols; j++)
		{
			printf("%f ", ACCESS(A, i, j));
		}
		printf("\n");
	}
	printf("\n");
}

/*
 * Copies matrix B into A.
 */
void copy(matrix *A, matrix *B)
{
	if (A->data && A->rows * A->cols != B->rows * B->cols)
	{
		free(A->data);
	}

	A->rows = B->rows;
	A->cols = B->cols;

	if (B->data)
	{
		A->data = malloc(sizeof(double) * B->rows * B->cols);
		for (int i = 0; i < B->rows * B->cols; i++)
		{
			A->data[i] = B->data[i];
		}
	}
}

/*
 * Initializes a Null Matrix meant to hold the memory for non-root nodes.
 */
void initNullMatrix(matrix *A, int rows, int cols)
{
	A->rows = rows;
	A->cols = cols;
	A->data = NULL;
}

/*
 * Compares the values of two matrices and returns 0 or 1.
 */
int isEqual(matrix *A, matrix *B)
{
	if (A->rows != B->rows || A->cols != B->cols)
	{
		return 0;
	}
	for (int i = 0; i < A->rows; i++)
	{
		for (int j = 0; j < A->cols; j++)
		{
			if (ACCESS(A, i, j) != ACCESS(B, i, j))
			{
				return 0;
			}
		}
	}
	return 1;
}

/*
 * Uses MPI_Scatter to compute addition in parallel.
 */
double *add(matrix *A, matrix *B, MPI_Comm world, int worldSize, int myRank)
{
	if ((A->rows != B->rows) || (A->cols != B->cols))
	{
		printf("Matrices must have same dimensions!\n");
		return 0;
	}

	int length = A->rows * A->cols; // total data in A and B, since dimensions are same
	int block = length / worldSize;
	int *sendcounts = malloc(worldSize * sizeof(int));
	int *displs = malloc(worldSize * sizeof(int));

	for (int i = 0; i < worldSize; i++)
	{
		sendcounts[i] = block;
	}
	for (int i = 0; i < length % worldSize; i++)
	{
		sendcounts[i]++;
	}
	displs[0] = 0;
	for (int i = 1; i < worldSize; i++)
	{
		displs[i] = displs[i - 1] + sendcounts[i - 1];
	}

	double *C = malloc(length * sizeof(double));	 //new matrix
	double *Atemp = malloc(length * sizeof(double)); //temporary variable for A's values
	double *Btemp = malloc(length * sizeof(double)); //temporary variable for B's values
	double *sum = malloc(length * sizeof(double));

	MPI_Scatterv(
		A->data, sendcounts, displs, MPI_DOUBLE,
		Atemp, sendcounts[myRank], MPI_DOUBLE,
		0, world);

	MPI_Scatterv(
		B->data, sendcounts, displs, MPI_DOUBLE,
		Btemp, sendcounts[myRank], MPI_DOUBLE,
		0, world);

	int blocksize = sendcounts[myRank];
	int blockbytes = blocksize * sizeof(double);
	int cache = blockbytes / 32000;
	if (blockbytes < 32000)
	{
		for (int i = 0; i < blocksize; i++)
		{
			sum[i] = Atemp[i] + Btemp[i];
		}
	}
	else
	{
		for (int i = 0; i < blockbytes * cache; i += cache)
		{
			for (int j = i; j < i + cache && j < blocksize; j++)
			{
				sum[i] = Atemp[i] + Btemp[i];
			}
		}
	}

	MPI_Gatherv(
		sum, sendcounts[myRank], MPI_DOUBLE,
		C, sendcounts, displs, MPI_DOUBLE,
		0, world);

	free(sendcounts);
	free(displs);
	free(Atemp);
	free(Btemp);
	free(sum);
	return C;
}

/*
 * Uses MPI_Scatter to compute subtraction in parallel.
 */

double *subtract(matrix *A, matrix *B, MPI_Comm world, int worldSize, int myRank)
{
	if ((A->rows != B->rows) || (A->cols != B->cols))
	{
		printf("Matrices must have same dimensions!\n");
		return 0;
	}

	int length = A->cols * A->rows;
	int block = length / worldSize;

	int *sendcounts = malloc(worldSize * sizeof(int));
	int *displs = malloc(worldSize * sizeof(int));

	for (int i = 0; i < worldSize; i++)
	{
		sendcounts[i] = block;
	}
	for (int i = 0; i < length % worldSize; i++)
	{
		sendcounts[i]++;
	}
	displs[0] = 0;
	for (int i = 1; i < worldSize; i++)
	{
		displs[i] = displs[i - 1] + sendcounts[i - 1];
	}

	double *C = malloc(length * sizeof(double));	 //new matrix
	double *Atemp = malloc(length * sizeof(double)); //temporary variable for A's values
	double *Btemp = malloc(length * sizeof(double)); //temporary variable for B's values
	double *dif = malloc(length * sizeof(double));

	MPI_Scatterv(
		//send data: A's contents, # of elements each node gets, indicies, data type, dsplc
		A->data, sendcounts, displs, MPI_DOUBLE,
		//receive data: receive buffer, # of elements in buffer, data type
		Atemp, sendcounts[myRank], MPI_DOUBLE,
		//settings: root node, communicator
		0, world);

	MPI_Scatterv(
		//send data: B's contents, # of elements each node gets, indicies, data type, dsplc
		B->data, sendcounts, displs, MPI_DOUBLE,
		//receive data: receive buffer, # of elements in buffer, data type
		Btemp, sendcounts[myRank], MPI_DOUBLE,
		//settings: root node, communicator
		0, world);

	int blocksize = sendcounts[myRank];
	int blockbytes = blocksize * sizeof(double);
	int cache = blockbytes / 32000;
	if (blockbytes < 32000)
	{
		for (int i = 0; i < blocksize; i++)
		{
			dif[i] = Atemp[i] - Btemp[i];
		}
	}
	else
	{
		for (int i = 0; i < blockbytes * cache; i += cache)
		{
			for (int j = i; j < i + cache && j < blocksize; j++)
			{
				dif[i] = Atemp[i] - Btemp[i];
			}
		}
	}
	MPI_Gatherv(
		dif, sendcounts[myRank], MPI_DOUBLE,
		C, sendcounts, displs, MPI_DOUBLE,
		0, world);

	free(sendcounts);
	free(displs);
	free(Atemp);
	free(Btemp);
	free(dif);
	return C;
}

/*
 * Uses MPI_Scatter to compute multiplication in parallel.
 */

double *multiply(matrix *A, matrix *B, MPI_Comm world, int worldSize, int myRank)
{
	if ((A->cols != B->rows))
	{
		printf("First matrix columns must match second matrix rows.\n");
		return 0;
	}

	matrix rv;
	matrix cv;

	zeroMatrix(&rv, 1, A->cols);
	zeroMatrix(&cv, B->rows, 1);

	double *total = NULL;

	if (myRank == 0)
	{
		total = malloc(A->rows * B->cols * sizeof(double));
	}

	for (int i = 0; i < A->rows; i++)
	{
		for (int j = 0; j < B->cols; j++)
		{
			if (myRank == 0)
			{
				for (int k = 0; k < A->cols; k++)
				{
					rv.data[k] = ACCESS(A, i, k);
					//printf("rv.data[%d] = %0.2f\n", k, rv.data[k]);
				}
			}
			if (myRank == 0)
			{
				for (int k = 0; k < B->rows; k++)
				{
					cv.data[k] = ACCESS(B, k, j);
					//printf("cv.data[%d] = %0.2f\n", k, cv.data[k]);
				}
			}
			double product = innerProduct(&rv, &cv, world, worldSize, myRank);
			if (myRank == 0)
			{
				total[INDEX(B, i, j)] = product;
				//printf("Product of B[%d,%d] = %0.2f\n", i, j, product);
			}
		}
	}

	free(rv.data);
	free(cv.data);
	return total;
}

/*
 * Uses MPI_Scatter to compute inner products in parallel.
 */
double innerProduct(matrix *A, matrix *B, MPI_Comm world, int worldSize, int myRank)
{
	if (A->rows != 1 && A->cols != 1)
	{
		printf("First entry is not a vector!\n");
		return 1;
	}
	if (B->rows != 1 && B->cols != 1)
	{
		printf("Second entry is not a vector!\n");
		return 1;
	}
	int length = A->rows * A->cols;
	int block = length / worldSize;
	int *sendcounts = malloc(worldSize * sizeof(int));
	int *displs = malloc(worldSize * sizeof(int));

	for (int i = 0; i < worldSize; i++)
	{
		sendcounts[i] = block;
	}
	for (int i = 0; i < length % worldSize; i++)
	{
		sendcounts[i]++;
	}
	displs[0] = 0;
	for (int i = 1; i < worldSize; i++)
	{
		displs[i] = displs[i - 1] + sendcounts[i - 1];
	}

	double *Atemp = malloc(sendcounts[myRank] * sizeof(double));
	double *Btemp = malloc(sendcounts[myRank] * sizeof(double));

	MPI_Scatterv(A->data, sendcounts, displs, MPI_DOUBLE, Atemp, sendcounts[myRank], MPI_DOUBLE, 0, world);
	MPI_Scatterv(B->data, sendcounts, displs, MPI_DOUBLE, Btemp, sendcounts[myRank], MPI_DOUBLE, 0, world);

	double sum = 0;
	double totalSum = 0;

	int blocksize = sendcounts[myRank];
	int blockbytes = blocksize * sizeof(double);
	int cache = blockbytes / 32000;
	if (blockbytes < 32000)
	{
		for (int i = 0; i < blocksize; i++)
		{
			sum += Atemp[i] * Btemp[i];
		}
	}
	else
	{
		for (int i = 0; i < blockbytes * cache; i += cache)
		{
			for (int j = i; j < i + cache && j < blocksize; j++)
			{
				sum += Atemp[i] * Btemp[i];
			}
		}
	}

	MPI_Reduce(&sum, &totalSum, 1, MPI_DOUBLE, MPI_SUM, 0, world);

	free(sendcounts);
	free(displs);
	free(Atemp);
	free(Btemp);
	if (myRank == 0)
	{
		return totalSum;
	}
	return 1;
}

/*
 * Solves the system of equations Ax = b using Gauss-Jordan Elimination.
 */
double *Gauss_Jordan(matrix *A, matrix *B, MPI_Comm *world, int worldSize, int myRank)
{
	if (A->rows != B->rows)
	{
		printf("Need to have matching row sizes!\n");
		return NULL;
	}

	matrix a, b;
	if (myRank == 0)
	{
		zeroMatrix(&a, A->rows, B->rows);
		zeroMatrix(&b, B->rows, B->cols);
	}
	else
	{
		initNullMatrix(&a, A->rows, A->rows);
		initNullMatrix(&b, B->rows, B->cols);
	}
	matrix *Acopy = &a;
	matrix *Bcopy = &b;
	if (myRank == 0)
	{
		copy(Acopy, A);
		copy(Bcopy, B);
	}
	int sendCountsA[worldSize];
	int sendCountsB[worldSize];
	int dsplcA[worldSize];
	int dsplcB[worldSize];

	for (int i = 0; i < worldSize; i++)
	{
		sendCountsA[i] = (Acopy->rows / worldSize) * Acopy->cols;
		sendCountsB[i] = (Bcopy->rows / worldSize) * Bcopy->cols;
	}
	for (int i = 0; i < (Acopy->rows % worldSize); i++)
	{
		sendCountsA[i] += Acopy->cols;
		sendCountsB[i] += Bcopy->cols;
	}

	dsplcA[0] = dsplcB[0] = 0;
	for (int i = 1; i < worldSize; i++)
	{
		dsplcA[i] = dsplcA[i - 1] + sendCountsA[i - 1];
		dsplcB[i] = dsplcB[i - 1] + sendCountsB[i - 1];
	}

	double *Atemp = (double *)malloc(sendCountsA[myRank] * sizeof(double));
	double *Btemp = (double *)malloc(sendCountsB[myRank] * sizeof(double));

	double l[Acopy->cols];
	double Arow[Acopy->cols];
	double Brow[Bcopy->cols];
	MPI_Scatterv(Acopy->data, sendCountsA, dsplcA, MPI_DOUBLE, Atemp, sendCountsA[myRank], MPI_DOUBLE, 0, *world);
	MPI_Scatterv(Bcopy->data, sendCountsB, dsplcB, MPI_DOUBLE, Btemp, sendCountsB[myRank], MPI_DOUBLE, 0, *world);
	for (int k = 0; k < Acopy->rows; k++)
	{
		if (myRank == 0)
		{
			for (int i = 0; i < Acopy->rows; i++)
			{
				l[i] = ((double)ACCESS(Acopy, i, k)) / ((double)ACCESS(Acopy, k, k));
			}
			for (int i = 0; i < Acopy->cols; i++)
			{
				Arow[i] = ACCESS(Acopy, k, i);
			}
			for (int i = 0; i < Bcopy->cols; i++)
			{
				Brow[i] = ACCESS(Bcopy, k, i);
			}
		}

		MPI_Bcast(&l, Acopy->cols, MPI_DOUBLE, 0, *world);
		MPI_Bcast(&Arow, Acopy->cols, MPI_DOUBLE, 0, *world);
		MPI_Bcast(&Brow, Bcopy->cols, MPI_DOUBLE, 0, *world);

		int offset = dsplcA[myRank] / Acopy->cols;

		for (int r = 0; r < (sendCountsA[myRank] / Acopy->cols); r++)
		{
			if (k == r + offset)
			{
				continue;
			}
			for (int c = 0; c < Acopy->cols; c++)
			{
				Atemp[INDEX(Acopy, r, c)] = Atemp[INDEX(Acopy, r, c)] - (l[r + offset] * Arow[c]);
			}
			for (int c = 0; c < Bcopy->cols; c++)
			{
				Btemp[INDEX(Bcopy, r, c)] = Btemp[INDEX(Bcopy, r, c)] - (l[r + offset] * Brow[c]);
			}
		}

		if (myRank == 0)
		{
			free(Acopy->data);
			free(Bcopy->data);
			Acopy->data = (double *)malloc(Acopy->rows * Acopy->cols * sizeof(double));
			Bcopy->data = (double *)malloc(Bcopy->rows * Bcopy->cols * sizeof(double));
		}
		MPI_Gatherv(Atemp, sendCountsA[myRank], MPI_DOUBLE, Acopy->data, sendCountsA, dsplcA, MPI_DOUBLE, 0, *world);
		MPI_Gatherv(Btemp, sendCountsB[myRank], MPI_DOUBLE, Bcopy->data, sendCountsB, dsplcB, MPI_DOUBLE, 0, *world);
	}
	if (myRank == 0)
	{
		double L[Acopy->cols];
		for (int i = 0; i < Acopy->rows; i++)
		{
			L[i] = ACCESS(Acopy, i, i);
		}
		for (int i = 0; i < Acopy->rows; i++)
		{
			for (int j = 0; j < Acopy->cols; j++)
			{
				ACCESS(Acopy, i, j) = ACCESS(Acopy, i, j) / L[i];
			}
		}
		for (int i = 0; i < Bcopy->rows; i++)
		{
			for (int j = 0; j < Bcopy->cols; j++)
			{
				ACCESS(Bcopy, i, j) = ACCESS(Bcopy, i, j) / L[i];
			}
		}
	}

	free(Acopy->data);
	free(Atemp);
	free(Btemp);

	if (myRank == 0)
	{
		return Bcopy->data;
	}
	return NULL;
}

/*
 * Normalizes a vector x ==> x = x/||x|| = x/sqrt(x[0]^2 + x[1]^2 + ... + x[n-1]^2)
 */
double norm(matrix *A, MPI_Comm world, int worldSize, int myRank)
{
	if (A->rows > 1 && A->cols != 1)
	{
		printf("Input must be a vector!\n");
		return -1;
	}
	if (A->cols > 1 && A->rows != 1)
	{
		printf("Input must be a vector!\n");
		return -1;
	}

	if (A->rows == 1)
	{
		double *temp = transpose(A);
		free(A->data);
		A->data = temp;
	}

	double TotalSum = 0;
	double Sum = 0;

	int sendCounts[worldSize];
	int dsplc[worldSize];

	int s = 0;
	for (int i = 0; i < worldSize; i++)
	{
		sendCounts[i] = (A->cols * A->rows) / worldSize;
	}
	for (int i = 0; i < (A->cols * A->rows) % worldSize; i++)
	{
		sendCounts[i] += 1;
	}
	for (int i = 0; i < worldSize; i++)
	{
		dsplc[i] = s;
		s += sendCounts[i];
	}

	double values[sendCounts[myRank]];

	MPI_Scatterv(A->data, sendCounts, dsplc, MPI_DOUBLE, values, sendCounts[myRank], MPI_DOUBLE, 0, world);

	for (int i = 0; i < sendCounts[myRank]; i++)
	{
		Sum += pow(values[i], 2);
	}
	//Does it matter to have the norm on all processors?
	MPI_Reduce(&Sum, &TotalSum, 1, MPI_DOUBLE, MPI_SUM, 0, world);

	return sqrt(TotalSum);
}

/*
 * Repeatedly multiplies A*x = x/||x|| to find the Eigenvector of A within a tolerance of 10^-16.
 */
double *EigenVector(int rows, int cols, MPI_Comm world, int worldSize, int myRank)
{
	MPI_File fh;
	matrix *A;
	matrix test;
	if (myRank == 0)
	{
		zeroMatrix(&test, rows, cols);
	}
	else
	{
		initNullMatrix(&test, rows, cols);
	}
	A = &test;
	if (myRank == 0)
	{
		zeroMatrix(A, rows, cols);
	}
	else
	{
		initNullMatrix(A, rows, cols);
	}

	MPI_File_open(world, "Matrix.data", MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
	if (myRank == 0)
	{
		MPI_File_read(fh, A->data, rows * cols, MPI_DOUBLE, MPI_STATUS_IGNORE);
	}
	MPI_File_close(&fh);
	MPI_Barrier(world);

	matrix temp1, temp2;
	if (myRank == 0)
	{
		zeroMatrix(&temp1, A->rows, 1);
		zeroMatrix(&temp2, A->rows, 1);
		for (int i = 0; i < temp1.rows * temp1.cols; i++)
		{
			temp1.data[i] = 1;
			temp2.data[i] = 1;
		}
	}
	else
	{
		initNullMatrix(&temp1, A->rows, 1);
		initNullMatrix(&temp2, A->cols, 1);
	}
	int maxcount = 0;
	double tol = .0000000000000001;
	int done = 0;
	double *dif;
	double length;
	while (maxcount < 10000 && done == 0)
	{
		free(temp2.data);
		temp2.data = temp1.data;

		temp1.data = multiply(A, &temp1, world, worldSize, myRank);

		length = norm(&temp1, world, worldSize, myRank);

		MPI_Bcast(&length, 1, MPI_DOUBLE, 0, world);

		int sendCounts[worldSize];
		int dsplc[worldSize];

		for (int i = 0; i < worldSize; i++)
		{
			sendCounts[i] = (temp1.rows * temp1.cols) / worldSize;
		}
		for (int i = 0; i < (temp1.rows * temp1.cols) % worldSize; i++)
		{
			sendCounts[i] += 1;
		}

		int s = 0;
		for (int i = 0; i < worldSize; i++)
		{
			dsplc[i] = s;
			s += sendCounts[i];
		}

		double values[sendCounts[myRank]];
		MPI_Scatterv(temp1.data, sendCounts, dsplc, MPI_DOUBLE, values, sendCounts[myRank], MPI_DOUBLE, 0, world);

		for (int i = 0; i < sendCounts[myRank]; i++)
		{
			values[i] /= length;
		}

		MPI_Gatherv(values, sendCounts[myRank], MPI_DOUBLE, temp1.data, sendCounts, dsplc, MPI_DOUBLE, 0, world);

		dif = subtract(&temp1, &temp2, world, worldSize, myRank);

		if (myRank == 0)
		{
			done = 1;
			for (int i = 0; i < temp1.rows * temp2.cols; i++)
			{
				if ((dif[i] > 0 ? dif[i] : dif[i] * -1) > tol)
				{
					done = 0;
					break;
				}
			}
		}
		if (myRank == 0)
		{
			free(dif);
		}

		MPI_Bcast(&done, 1, MPI_DOUBLE, 0, world);

		maxcount++;
	}

	free(A->data);
	free(temp2.data);
	if (myRank == 0)
	{
		return temp1.data;
	}
	return NULL;
}

/*
 * Calculates the tranpose of a matrix.
 */
double *transpose(matrix *A)
{
	double *t = malloc(A->rows * A->cols * sizeof(double));

	for (int i = 0; i < A->cols; i++)
	{
		for (int j = 0; j < A->rows; j++)
		{
			t[i * A->rows + j] = ACCESS(A, j, i);
		}
	}
	return t;
}



void split_matrix(int *sendcts, int *displs, int total, int worldsize)
{
	int rem = total % worldsize;
	displs[0] = 0;
	for (int i = 0; i < worldsize; i++)
	{
		sendcts[i] = total / worldsize;
		if (rem > 0)
		{
			sendcts[i]++;
			rem--;
		}
		if (i > 0)
		{
			displs[i] = displs[i - 1] + sendcts[i - 1];
		}
	}
}

void create_ones_vector(matrix *vector, int rank, int worldsize, MPI_Comm world)
{
	if (rank == 0)
	{
		srand(time(NULL));
		for (int i = 0; i < vector->rows; i++)
		{
			for (int j = 0; j < vector->cols; j++)
			{

				vector->data[LENSINDEX(vector->rows, vector->cols, i, j)] = 1;

				//     printf("inserting %.2f, at index ( %d %d)\n",ins,i,j);

				//   printf("inserted %.2f, at index ( %d %d)\n",mat->data[INDEX(mat->rows, mat->cols, i, j)],i,j);
			}
			//  printf("\n");
		}
	}
}
void scalar_multiply(matrix *lhs, double scalar, MPI_Comm world, int worldSize, int myRank)
{
	/* 
	performs scalar multiplication in parallel with scatter v and gather v does not return changes matrix passed
	*/
	int *sendcts = malloc(sizeof(int) * worldSize);
	int *displs = malloc(sizeof(int) * worldSize);
	split_matrix(sendcts, displs, lhs->rows * lhs->cols, worldSize);
	double *receive = malloc(sendcts[myRank] * sizeof(double));

	MPI_Scatterv(lhs->data, sendcts, displs, MPI_DOUBLE, receive, sendcts[myRank], MPI_DOUBLE, 0, world);
	for (int i = 0; i < sendcts[myRank]; i++)
	{
		receive[i] *= scalar;
	}

	MPI_Gatherv(receive, sendcts[myRank], MPI_DOUBLE, lhs->data, sendcts, displs, MPI_DOUBLE, 0, world);
}

void matrix_vect_mult(matrix *lhs, matrix *vector, matrix *final, MPI_Comm world, int worldsize, int rank)
{
	int *sendcts = malloc(sizeof(int) * worldsize);
	int *displs = malloc(sizeof(int) * worldsize);
	split_matrix(sendcts, displs, lhs->rows, worldsize);
	for (int i = 0; i < worldsize; i++)
	{
		sendcts[i] *= lhs->cols;
		displs[i] *= lhs->cols;
	}
	//now each process has however many rows
	MPI_Barrier(world);
	double *lrow = malloc(sizeof(double) * sendcts[rank]);

	/*need to change to scatter and gather alg*/
	MPI_Scatterv(
		lhs->data,
		sendcts,
		displs,
		MPI_DOUBLE,
		lrow,
		sendcts[rank],
		MPI_DOUBLE,
		0,
		world);
	for (int i = 0; i < worldsize; i++)
	{
		sendcts[i] /= lhs->cols;
		displs[i] /= lhs->cols;
	}
	double *adder = malloc(sendcts[rank] * sizeof(double));

	MPI_Bcast(vector->data, vector->rows * vector->cols, MPI_DOUBLE, 0, world);

	for (int i = 0; i < sendcts[rank]; i++)
	{ //for each recieved row
		adder[i] = 0;
		for (int j = 0; j < lhs->cols; j++)
		{
			adder[i] += (lrow[LENSINDEX(sendcts[rank] / lhs->cols, lhs->cols, i, j)] * vector->data[LENSINDEX(vector->rows, 1, i, j)]);
		}
		//   printf("%lf\n", adder[i]);
	}
	MPI_Gatherv(adder, sendcts[rank], MPI_DOUBLE, final->data, sendcts, displs, MPI_DOUBLE, 0, world);

	free(sendcts);
	free(displs);
	free(lrow);
}
double *page_rank(matrix *M, MPI_Comm world, int worldSize, int myRank)
{
	/*todo: testing this
	for normalizing the adj matrix we should transpose than take the sum of the row using MPI then divide each index by the sum
	calcuolates pagerank by inverse(identity-sM)se s being the scalar random jump variable and e being the all ones vector
	allow s to be dynamic
	:endtodo
	*/
	matrix Identity, e, final;
	final.rows = M->rows;
	final.cols = 1;
	final.data = malloc(M->rows * sizeof(double));
	e.rows = M->rows;
	e.cols = 1;
	e.data = malloc(M->rows * sizeof(double));
	create_ones_vector(&e, myRank, worldSize, world);
	MPI_Bcast(M, M->rows * M->cols, MPI_DOUBLE, 0, world);
	idenMatrix(&Identity, M->rows, M->cols);
	int *norms = malloc(sizeof(int) * M->cols);
	// normalize it so that columns have a sum of 1, making it a stochastic matrix need to make this parallel
	for (int i = 0; i < M->cols; i++)
	{
		norms[i] = 0;
	}
	if (myRank == 0)
	{
		for (int i = 0; i < M->cols; i++)
		{

			for (int j = 0; j < M->rows; j++)
			{

				if (ACCESS(M, j, i) == 1)
				{
					norms[i]++;
				}
			}
			//	printf("norm for col %d: %d\n", i, norms[i]);
		}
		for (int i = 0; i < M->cols; i++)
		{
			for (int j = 0; j < M->rows; j++)
			{
				ACCESS(M, j, i) =ACCESS(M, j, i) / norms[i];
			}
		}
		// printf("M= after norm\n");
		// for (int i = 0; i < M->rows; i++)
		// {
		// 	printf("[");
			
		// 	for (int j = 0; j < M->cols; j++)
		// 	{
		// 		printf(" %lf ", ACCESS(M, i, j));
		// 	}
		// 				printf("]\n");

		// }

	}
	MPI_Bcast(M, M->rows * M->cols, MPI_DOUBLE, 0, world);
	double s = .85; //al

	scalar_multiply(M, s, world, worldSize, myRank); //now M=sM
	MPI_Bcast(M, M->rows * M->cols, MPI_DOUBLE, 0, world);
	M->data = subtract(&Identity, M, world, worldSize, myRank); //I-sm
	MPI_Barrier(world);
	M->data = Gauss_Jordan(M, &Identity, &world, worldSize, myRank);
	scalar_multiply(&e, s, world, worldSize, myRank);
	// if (myRank == 0)
	// {
	// 	printf("e=\n");
	// 	for (int i = 0; i < M->rows; i++)
	// 	{
	// 		printf("[");
	// 		printf(" %lf ", e.data[i]);
	// 		printf("]\n");
	// 	}
	// 	printf("M=\n");
	// 	for (int i = 0; i < M->rows; i++)
	// 	{
	// 		printf("[");
			
	// 		for (int j = 0; j < M->cols; j++)
	// 		{
	// 			printf(" %lf ", ACCESS(M, i, j));
	// 		}
	// 					printf("]\n");

	// 	}
	// }
	matrix_vect_mult(M, &e, &final, world, worldSize, myRank);

	if (myRank == 0)
	{
		// printf("M page ranks=\n");
		// for (int i = 0; i < M->rows; i++)
		// {
		// 	printf("[");
		// 	printf(" %lf ", final.data[i]);
		// 	printf("]\n");
		// }
		return final.data;
	}
	return NULL;
}
void create_random_matrix(matrix *mat, int limit, int rank, MPI_Comm world)
{
	/*todo make parallel */
	if (rank == 0)
	{
		srand(time(NULL));
		for (int i = 0; i < mat->rows; i++)
		{
			for (int j = 0; j < mat->cols; j++)
			{
				double ins = rand() % limit;
				//     printf("inserting %.2f, at index ( %d %d)\n",ins,i,j);

				mat->data[LENSINDEX(mat->rows, mat->cols, i, j)] = ins;
				//   printf("inserted %.2f, at index ( %d %d)\n",mat->data[INDEX(mat->rows, mat->cols, i, j)],i,j);
			}
			//  printf("\n");
		}
	}
}

/**
 * Reads in a matrix from the given fileName with dimensions rows X cols.
 * Uses MPI to divide the matrix into an even number of rows and then extract the non-zero elements along with their
 * row index and column index, and then writes these values to the file Dense.data in the following format:
 * 
 * Dense.data: | totalNonZero * sizeof(double) | totalNonZero * sizeof(int) | totalNonZero * sizeof(int) |
 *             |        Non-Zero Values        | Corresponding row Indices  | Corresponding col Indices  |
 *
 * Returns the number of non zero values from the original matrix so the output file can be parsed later
 * @author Olivia Rorke
*/
int toDenseMatrix(char* fileName, int rows, int cols, MPI_Comm world, int worldSize, int myRank) {
	
	MPI_File fh, fhDense;
	
	//Divide the matrix up as evenly as possibly by number of rows
	int* sendcts = malloc(worldSize * sizeof(int));
	int* offsets = malloc(worldSize * sizeof(int));
	
    //Determine how many rows to send to each node
    int minSend = rows / worldSize;
    for(int i = 0; i < worldSize; i++) {
        sendcts[i] = minSend;
		if(i < rows % worldSize) {
			sendcts[i]++;
		}
		sendcts[i] *= cols;
    }
	
	//Determine offset into file for reading
	offsets[0] = 0;
	for(int i = 1; i < worldSize; i++) {
		offsets[i] = offsets[i - 1] + sendcts[i - 1];
	}
	
	//Scatter the data
	int numRead, offset;
	MPI_Scatter(sendcts, 1, MPI_INT, &numRead, 1, MPI_INT, 0, world);
	MPI_Scatter(offsets, 1, MPI_INT, &offset, 1, MPI_INT, 0, world);
	
	//Read in the section of the matrix this node is responsible for
	double *temp = malloc(numRead * sizeof(double));
	MPI_File_open(world, fileName, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
	MPI_File_read_at(fh, offset * sizeof(double), temp, numRead, MPI_DOUBLE, MPI_STATUS_IGNORE);
	MPI_File_close(&fh);
	
	//Each node condenses their section of the matrix
	int rowsRecv = numRead / cols;
	int globalStartRow = offset / cols;
	int numNonZero = 0;
	
	for(int r = 0; r < rowsRecv; r++) {
		for(int c = 0; c < cols; c++) {
			double curr = temp[(r * cols) + c];
			if(curr != 0) {
				numNonZero++;
			}
		}
	}
	
	double *values = malloc(numNonZero * sizeof(double));
	int *rowIdx = malloc(numNonZero * sizeof(int));
	int *colIdx = malloc(numNonZero * sizeof(int));
	int counter = 0;
	
	for(int r = 0; r < rowsRecv; r++) {
		for(int c = 0; c < cols; c++) {
			double curr = temp[(r * cols) + c];
			if(curr != 0) {
				values[counter] = curr;
				rowIdx[counter] = r + globalStartRow;
				colIdx[counter] = c;
				counter++;
			}
		}
	}

	free(temp);
	
	//Report to the root how many non zero values were found
	int *nonZeroCounts = NULL;
	int totalNonZero;
	if(myRank == 0) {
		nonZeroCounts = malloc(worldSize * sizeof(int));
	}
	MPI_Gather(&numNonZero, 1, MPI_INT, nonZeroCounts, 1, MPI_INT, 0, world);
	
	//Let all nodes know how many non zero values were found
	MPI_Reduce(&numNonZero, &totalNonZero, 1, MPI_INT, MPI_SUM, 0, world);
	MPI_Bcast(&totalNonZero, 1, MPI_INT, 0, world);
	
	//Create offsets for the new dense file
	int* denseOffsets = malloc(worldSize * sizeof(int));
	denseOffsets[0] = 0;
	for(int i = 1; i < worldSize; i++) {
		denseOffsets[i] = denseOffsets[i - 1] + nonZeroCounts[i - 1];
	}
	int denseOffset;
	MPI_Scatter(denseOffsets, 1, MPI_INT, &denseOffset, 1, MPI_INT, 0, world);
	
	//Write the dense matrix contents to the new file
	MPI_File_open(world, "Dense.data", MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fhDense);
	MPI_File_write_at(fhDense, denseOffset, values, numNonZero, MPI_DOUBLE, MPI_STATUS_IGNORE);
	MPI_File_write_at(fhDense, denseOffset + (totalNonZero * sizeof(double)), rowIdx, numNonZero, MPI_INT, MPI_STATUS_IGNORE);
	MPI_File_write_at(fhDense, denseOffset + (totalNonZero * (sizeof(double) + sizeof(int))), colIdx, numNonZero, MPI_INT, MPI_STATUS_IGNORE);
	MPI_File_close(&fhDense);
	
	//Free data
	free(sendcts);
	free(offsets);
	free(denseOffsets);
	free(nonZeroCounts);
	free(values);
	free(rowIdx);
	free(colIdx);
	
	//Return the number of non-zero elements found in the matrix
	return totalNonZero;
}

/** 
 * Reads in the data from a condensed matrix
 * @param A - a cMatrix object with the numEntries and fileName fields populated
 * After completion the cMatrix passed will contain the values, rows, and columns 
 * of the sparse matrix data contained in the specified file.
 *
 * @author Oliva Rorke
 */
void readCMatrix(cMatrix *A, MPI_Comm world, int worldSize, int myRank) {

	MPI_File fh;
	
	//Divide the matrix up as evenly as possibly by number of rows
	int* sendcts = malloc(worldSize * sizeof(int));
	int* offsets = malloc(worldSize * sizeof(int));
	
    //Determine how many rows to send to each node
    int minSend = A->numEntries / worldSize;
    for(int i = 0; i < worldSize; i++) {
        sendcts[i] = minSend;
		if(i < A->numEntries % worldSize) {
			sendcts[i]++;
		}
    }
	
	//Determine offset into file for reading
	offsets[0] = 0;
	for(int i = 1; i < worldSize; i++) {
		offsets[i] = offsets[i - 1] + sendcts[i - 1];
	}
	
	//Scatter the data
	int numRead, offset;
	MPI_Scatter(sendcts, 1, MPI_INT, &numRead, 1, MPI_INT, 0, world);
	MPI_Scatter(offsets, 1, MPI_INT, &offset, 1, MPI_INT, 0, world);
	
	//Read in the section of the matrix this node is responsible for
	double *valuesLocal = malloc(numRead * sizeof(double));
	int *rowsLocal = malloc(numRead * sizeof(int));
	int *colsLocal = malloc(numRead * sizeof(int));
	
	MPI_File_open(world, A->fileName, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
	MPI_File_read_at(fh, offset * sizeof(double), valuesLocal, numRead, MPI_DOUBLE, MPI_STATUS_IGNORE);
	MPI_File_read_at(fh, (offset * sizeof(int)) + (A->numEntries * sizeof(double)), rowsLocal, numRead, MPI_INT, MPI_STATUS_IGNORE);
	MPI_File_read_at(fh, (offset * sizeof(int)) + (A->numEntries * (sizeof(int) + sizeof(double))), colsLocal, numRead, MPI_INT, MPI_STATUS_IGNORE);
	MPI_File_close(&fh);
	
	//Gather the data into root node
	double* values = NULL;
	int* rows = NULL;
	int* cols = NULL;
	if(myRank == 0) {
		values = malloc(A->numEntries * sizeof(double));
		rows = malloc(A->numEntries * sizeof(int));
		cols = malloc(A->numEntries * sizeof(int));
	}
	
	MPI_Gatherv(valuesLocal, numRead, MPI_DOUBLE, values, sendcts, offsets, MPI_DOUBLE, 0, world);
	MPI_Gatherv(rowsLocal, numRead, MPI_INT, rows, sendcts, offsets, MPI_INT, 0, world);
	MPI_Gatherv(colsLocal, numRead, MPI_INT, cols, sendcts, offsets, MPI_INT, 0, world);
	
	//Update the matrix being returned
	A->values = values;
	A->rows = rows;
	A->cols = cols;
	
	//Free data
	free(sendcts);
	free(offsets);
	free(valuesLocal);
	free(rowsLocal);
	free(colsLocal);
}

/*
 * Reads in a file and computes the adjacency matrix of it.
 * Note: This is using the "arxiv-citations.txt" file as a 
 *       reference.
 */
double * Adjacency(matrix *A, char * file, int papers)
{
	FILE * fh;
       	if ((fh = fopen(file, "r")) < 0)
	{
		printf("Error: File failed to open!\n");
		return NULL;
	}
	else if (count != A->rows)
	{
		printf("Error: Number of papers does not match the number of rows of A!\n");
		return NULL;
	}
	else if (A->rows || A->rows)
	{
		printf("Error: A must be a square matrix!\n");
		return NULL;
	}

	char *curr = malloc(sizeof(char) * 17); //max id length + 1
	char *prev = malloc(sizeof(char) * 17); //to check the previous lines
	char *idmap[papers];

	int lineCount = 0;

	//three types of lines to check for:
	//1. +++++ -> separates different authors
	//2. ----- -> author above cites author's below
	//3. authorid's
	
	int count = 0;

	while (fscanf(fh, "%[^\n]\n", curr))
	{
		if (strcmp(curr, "+++++") == 0) 
		{
			//we have a line of +'s, need to move to next iteration
			count = 0;
			continue;
		}	
		if (count == 0)
		{
			//makes an array such that: id -> {id, id, id, id, id, ..}
			idmap[count] = malloc(sizeof(char)*length); 
			strcpy(idmap[count], curr); //to preserve the 'mapping' id
			count++;
		}
		if (strcmp(curr, "-----") == 0)
		{
			//this means we are starting to find the citations
				
		}
		
		//What I'm trying to do:
		// 1. If we see a plus line, record it and skip.
		// 2. If we see an id line, and the previous line is a plus line, 
		//    we know that we have the vertex we're counting for.
		//    We then want to store it and find out what it's citations are.
		//    
		//    We will also know the next line is a dash line. If we have
		//    non-plus lines after the dash line, we want to store that id
		//    in a list of lists of ids.
		//    For example vertex[id-0] = [id-1, id-2, id-3]. 
		//    vertex[id-1] = []
		//    vertex[id-2] = []
		//    vertex[id-3] = [id-2]
		//    Which gives us A:
		//    0 1 1 1
		//    0 0 0 0
		//    0 0 0 0
		//    0 0 1 0


		strcpy(prev, curr);
		


		
	}
	
}



