//Program: matrix.c
//Author: Mitchell Krystofiak
//Description: Implementation of the matrix class
//Date: September 24, 2021

#include"matrix.h"

/*
 * Remember: the structure of matrix typdef is:
 * 	int rows, int cols, double * data
 *
 */ 

/*
 * Initializes a matrix with random values.
 */

void randMatrix(matrix *A, int rows, int cols)
{
	A->rows = rows;
	A->cols = cols;
	A->data = malloc(A->rows*A->cols*sizeof(double)); 
	for (int i = 0; i < A->rows; i++)
	{
		for (int j = 0; j < A->cols; j++)
		{
			ACCESS(A,i,j) = rand() % 100 + 1;
		}
	}	
}

/*
 * Initializaed a matrix with zero values.
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
			ACCESS(A,i,j) = 0;
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
			printf("%0.2f ", ACCESS(A,i,j));
		}
		printf("\n");
	}
	printf("\n");
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
	
	//printf("Rank %d SHOULD receive %d things\n", myRank, sendcounts[myRank]);
		
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

	for (int i = 0; i < sendcounts[myRank]; i++)
	{
		//printf("Rank %d adding %f + %f\n", myRank, Atemp[i], Btemp[i]);
		sum[i] = Atemp[i] + Btemp[i];
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

	for (int i = 0; i < sendcounts[myRank]; i++)
	{
		dif[i] = Atemp[i] - Btemp[i];
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
	for (int i = 0; i < sendcounts[myRank]; i++)
	{
		sum += Atemp[i] * Btemp[i];
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
