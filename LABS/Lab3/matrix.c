//Program: matrix.c
//Author: Mitchell Krystofiak
//Description: Implementation of the matrix class
//Date: October 31, 2021

#include"matrix.h"

/*
 * Remember: the structure of matrix typdef is:
 * 	int rows, int cols, double * data
 */ 

/*
 * Initializes a matrix with random values.
 */

void randMatrix(matrix *A, int rows, int cols, int c)
{
	srand(time(0));
	A->rows = rows;
	A->cols = cols;
	A->data = malloc(A->rows*A->cols*sizeof(double)); 
	
	if (c == 1) //int
	{
		for (int i = 0; i < A->rows; i++)
		{
			for (int j = 0; j < A->cols; j++)
			{
				ACCESS(A,i,j) = rand() % 100; 
			}
		}
	}
	else if (c == 0) //double
	{
		for (int i = 0; i < A->rows; i++)
		{
			for (int j = 0; j < A->cols; j++)
			{
				ACCESS(A,i,j) = (((double)(rand() % 100))/((double)(rand() % 50) + 1));
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
	A->data = malloc(A->rows*A->cols*sizeof(double));
	for (int i = 0; i < A->rows; i++)
	{
		for (int j = 0; j < A->cols; j++)
		{
			ACCESS(A,i,j) = 0.0;
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
	A->data = malloc(A->rows*A->cols*sizeof(double));
	for (int i = 0; i < A->rows; i++)
	{
		for (int j = 0; j < A->cols; j++)
		{
			if (i == j)
			{
				ACCESS(A,i,j) = 1.0;
			}
			else
			{
				ACCESS(A,i,j) = 0.0;
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
			printf("%f ", ACCESS(A,i,j));
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
	if (A->data && A->rows*A->cols != B->rows*B->cols)
	{
		free(A->data);
	}

	A->rows = B->rows;
	A->cols = B->cols;

	if (B->data)
	{
		A->data = malloc(sizeof(double)*B->rows*B->cols);
		for (int i = 0; i < B->rows*B->cols; i++)
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
			if (ACCESS(A,i,j) != ACCESS(B,i,j))
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
double * add(matrix *A, matrix *B, MPI_Comm world, int worldSize, int myRank)
{
	if ((A->rows != B->rows) || (A->cols != B->cols))
	{
		printf("Matrices must have same dimensions!\n");
		return 0;
	}
	
	int length = A->rows * A->cols; // total data in A and B, since dimensions are same
	int block = length / worldSize;
	int *sendcounts = malloc(worldSize*sizeof(int));
	int *displs = malloc(worldSize*sizeof(int));	
	
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
		displs[i] = displs[i-1] + sendcounts[i-1];
	}
	
	double* C = malloc(length*sizeof(double));     //new matrix
	double* Atemp = malloc(length*sizeof(double)); //temporary variable for A's values
	double* Btemp = malloc(length*sizeof(double)); //temporary variable for B's values
	double* sum   = malloc(length*sizeof(double));
	
		
	MPI_Scatterv(
		A->data, sendcounts, displs, MPI_DOUBLE, 
		Atemp, sendcounts[myRank], MPI_DOUBLE,	 
		0, world);
	
	MPI_Scatterv(
		B->data, sendcounts, displs, MPI_DOUBLE, 
		Btemp, sendcounts[myRank], MPI_DOUBLE,       
		0, world);
	
	int blocksize = sendcounts[myRank];
	int blockbytes = blocksize*sizeof(double);
	int cache = blockbytes/32000;
	if ( blockbytes < 32000)
	{
		for (int i = 0; i < blocksize; i++)
		{
			sum[i] = Atemp[i] + Btemp[i];
		}
	}
	else
	{
		for (int i = 0; i <blockbytes*cache; i+=cache)
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

double * subtract(matrix *A, matrix *B, MPI_Comm world, int worldSize, int myRank)
{
	if ((A->rows != B->rows) || (A->cols != B->cols))
	{
		printf("Matrices must have same dimensions!\n");
		return 0;
	}

	int length = A->cols*A->rows;
	int block = length / worldSize;
	
	int *sendcounts = malloc(worldSize*sizeof(int));
	int *displs = malloc(worldSize*sizeof(int));
	
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
		displs[i] = displs[i-1] + sendcounts[i-1];
	}
	
	double* C = malloc(length*sizeof(double));     //new matrix
	double* Atemp = malloc(length*sizeof(double)); //temporary variable for A's values
	double* Btemp = malloc(length*sizeof(double)); //temporary variable for B's values
	double* dif   = malloc(length*sizeof(double));

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
	int blockbytes = blocksize*sizeof(double);
	int cache = blockbytes/32000;
	if ( blockbytes < 32000)
	{
		for (int i = 0; i < blocksize; i++)
		{
			dif[i] = Atemp[i] - Btemp[i];
		}
	}
	else
	{
		for (int i = 0; i <blockbytes*cache; i+=cache)
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

double * multiply(matrix *A, matrix *B, MPI_Comm world, int worldSize, int myRank)
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
	
	double* total = NULL;

	if (myRank == 0)
	{
		total = malloc(A->rows*B->cols*sizeof(double));
	}

	for (int i = 0; i < A->rows; i++)
	{
		for (int j = 0; j < B->cols; j++)
		{
			if (myRank == 0)
			{
				for (int k = 0; k < A->cols; k++)
				{
					rv.data[k] = ACCESS(A,i,k);
					//printf("rv.data[%d] = %0.2f\n", k, rv.data[k]);
				}
			}
			if (myRank == 0)
			{
				for (int k = 0; k < B->rows; k++)
				{
					cv.data[k] = ACCESS(B,k,j);
					//printf("cv.data[%d] = %0.2f\n", k, cv.data[k]);
				}
			}
			double product = innerProduct(&rv, &cv, world, worldSize, myRank);
			if (myRank == 0)
			{
				total[INDEX(B,i,j)] = product;
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
	int length = A->rows *A->cols;
	int block = length / worldSize;
	int *sendcounts = malloc(worldSize*sizeof(int));
	int *displs = malloc(worldSize*sizeof(int));
	
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
		displs[i] = displs[i-1] + sendcounts[i-1];
	}
	

	double *Atemp = malloc(sendcounts[myRank]*sizeof(double));
	double *Btemp = malloc(sendcounts[myRank]*sizeof(double));

	MPI_Scatterv(A->data, sendcounts, displs, MPI_DOUBLE, Atemp, sendcounts[myRank], MPI_DOUBLE, 0, world);
	MPI_Scatterv(B->data, sendcounts, displs, MPI_DOUBLE, Btemp, sendcounts[myRank], MPI_DOUBLE, 0, world);
	
	double sum = 0;
	double totalSum = 0;

	int blocksize = sendcounts[myRank];
	int blockbytes = blocksize*sizeof(double);
	int cache = blockbytes/32000;
	if ( blockbytes < 32000)
	{
		for (int i = 0; i < blocksize; i++)
		{
			sum += Atemp[i] * Btemp[i];
		}
	}
	else
	{
		for (int i = 0; i <blockbytes*cache; i+=cache)
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
double * Gauss_Jordan(matrix *A, matrix *B, MPI_Comm *world, int worldSize,  int myRank)
{
	if (A->rows != B->rows)
	{
		printf("Need to have matching row sizes!\n");
		return NULL;
	}

	matrix a,b;
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
	matrix* Acopy = &a;
	matrix* Bcopy = &b;
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
		sendCountsA[i] = (Acopy->rows/worldSize) * Acopy->cols;
		sendCountsB[i] = (Bcopy->rows/worldSize) * Bcopy->cols;
	}
	for (int i = 0; i < (Acopy->rows % worldSize); i++)
	{
		sendCountsA[i] += Acopy->cols;
		sendCountsB[i] += Bcopy->cols;
	}

	dsplcA[0] = dsplcB[0] = 0;
	for (int i = 1; i < worldSize; i++)
	{
		dsplcA[i] = dsplcA[i-1] + sendCountsA[i-1];
		dsplcB[i] = dsplcB[i-1] + sendCountsB[i-1];
	}
	
	double* Atemp = (double*)malloc(sendCountsA[myRank]*sizeof(double));
	double* Btemp = (double*)malloc(sendCountsB[myRank]*sizeof(double));

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
				l[i] = ((double)ACCESS(Acopy,i,k))/((double)ACCESS(Acopy,k,k));
			}
			for (int i = 0; i < Acopy->cols; i++)
			{
				Arow[i] = ACCESS(Acopy,k,i);
			}
			for (int i = 0; i < Bcopy->cols; i++)
			{
				Brow[i] = ACCESS(Bcopy,k,i);
			}
		}
		
		MPI_Bcast(&l, Acopy->cols, MPI_DOUBLE, 0, *world);
		MPI_Bcast(&Arow, Acopy->cols, MPI_DOUBLE, 0, *world);
		MPI_Bcast(&Brow, Bcopy->cols, MPI_DOUBLE, 0, *world);

		int offset = dsplcA[myRank]/Acopy->cols;

		for (int r = 0; r < (sendCountsA[myRank] / Acopy->cols); r++)
		{
			if (k == r + offset)
			{
				continue;
			}
			for (int c = 0; c < Acopy->cols; c++)
			{
				Atemp[INDEX(Acopy,r,c)] = Atemp[INDEX(Acopy,r,c)] - (l[r+offset] * Arow[c]);
			}
			for (int c = 0; c < Bcopy->cols; c++)
			{
				Btemp[INDEX(Bcopy,r,c)] = Btemp[INDEX(Bcopy,r,c)] - (l[r+offset] * Brow[c]);
			}
		}

		if (myRank == 0)
		{
			free(Acopy->data);
			free(Bcopy->data);
			Acopy->data = (double*)malloc(Acopy->rows*Acopy->cols*sizeof(double));
			Bcopy->data = (double*)malloc(Bcopy->rows*Bcopy->cols*sizeof(double));
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
				ACCESS(Acopy,i,j) = ACCESS(Acopy,i,j) / L[i];
			}
		}
		for (int i = 0; i < Bcopy->rows; i++)
		{
			for (int j = 0; j < Bcopy->cols; j++)
			{
				ACCESS(Bcopy,i,j) = ACCESS(Bcopy,i,j) / L[i];
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
double * normalize(matrix *A, MPI_Comm world, int worldSize, int myRank)
{

}

/*
 * Repeatedly multiplies A*x = x/||x|| to find the Eigenvector of A within a tolerance of 10^-16.
 */
double * power_method(matrix *A, MPI_Comm * world, int worldSize, int myRank)
{
	
}

/*
 * Calculates the tranpose of a matrix.
 */
matrix transpose(matrix *A)
{
	matrix T;
	T.rows = A->cols;
	T.cols = A->rows;
	T.data = malloc(A->rows*A->cols*sizeof(double));
	
	for (int i = 0; i < A->cols; i++)
	{
		for (int j = 0; j < A->rows; j++)
		{
			T.data[i*A->rows + j] = ACCESS(A, j, i);
		}
	}
	return T;
}
