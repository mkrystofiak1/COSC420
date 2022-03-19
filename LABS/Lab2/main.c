//Author: Mitchell Krystofiak
//Class: COSC 420 - Lab2
//Date: September 24, 2021
//Description: A program to use scatter and gather to perform the dot product between two vectors.

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<mpi.h>
#include<time.h>
#include"matrix.h"

/* 
 * Function just for initializing the NULL matrices for non-root nodes
*/
void initNullMatrix(matrix *A, int rows, int cols)
{
	A->cols = cols;
	A->rows = rows;
	A->data = NULL;
}

int main(int argc, char ** argv) 
{
	srand(time(NULL));

	//world variable declarations

	MPI_Init(&argc, &argv);
	MPI_Comm world = MPI_COMM_WORLD;
	int worldSize;
	int myRank;

	MPI_Comm_size(world, &worldSize);
	MPI_Comm_rank(world, &myRank);

	//large, general inner products
	
	matrix v1, v2, v3, v4, v5, v6, v7, v8, v9, v10;

	//general operations
	
	matrix A, B, C, D, E, F, G, H, I, J, K;	

	//large operations
	
	matrix L1, L2, L3, L4, L5, L6, L7;

	//timing variables
	
	double t1, t2;
	if (myRank == 0)
	{	

		//Large vector, inner product declarations

		printf("Initializing a 20x1 vector, v1.\n\n");
		randMatrix(&v1, 20, 1);
		printMatrix(&v1);
		
		printf("Initializing a 20x1 vector, v2.\n\n");
		randMatrix(&v2, 20, 1);
		printMatrix(&v2);
		
		printf("Initializing a 500x1 vector, v3.\n\n");
		randMatrix(&v3, 500, 1);
			
		printf("Initializing a 500x1 vector, v4.\n\n");
		randMatrix(&v4, 500, 1);
		
		printf("Initializing a 3563x1 vector, v5.\n\n");
		randMatrix(&v5, 3563, 1);
		
		printf("Initializing a 3563x1 vector, v6.\n\n");
		randMatrix(&v6, 3563, 1);
		
		printf("Initializing a 7693x1 vector, v7.\n\n");
		randMatrix(&v7, 7693, 1);
		
		printf("Initializing a 7693x1 vector, v8.\n\n");
		randMatrix(&v8, 7693, 1);
		
		printf("Initializing a 15000x1 vector, v9.\n\n");
		randMatrix(&v9, 15000, 1);
		
		printf("Initializing a 15000x1 vector, v10.\n\n");
		randMatrix(&v10, 15000, 1);
		
		//General testing matrix declarations
		
		
		printf("Initializing a 3x3 matrix, A.\n\n");	
		zeroMatrix(&A, 3, 3);
		printMatrix(&A);
		
		printf("Initializing a 3x3 matrix, B.\n\n");
		randMatrix(&B, 3, 3);
		printMatrix(&B);		
		
		printf("Initializing a 3x3 matrix, C.\n\n");
		zeroMatrix(&C, 3, 3);
		printMatrix(&C);
		
		printf("Initializing a 3x3 matrix, D.\n\n");
		zeroMatrix(&D, 3, 3);
		printMatrix(&D);
		
		printf("Initializing a 2x4 matrix, E.\n\n");
		randMatrix(&E, 2, 4);
		printMatrix(&E);
		
		printf("Initializing a 4x3 matrix, F.\n\n");
		randMatrix(&F, 4, 3);
		printMatrix(&F);
		
		printf("Initializing a 10x10 matrix, G.\n\n");
		randMatrix(&G, 10, 10);
		printMatrix(&G);
		
		printf("Initializing a 10x1 matrix, H.\n\n");
		randMatrix(&H, 10, 1);
		printMatrix(&H);
		
		printf("Initializing a 2x3 matrix, I.\n\n");
		zeroMatrix(&I, 2, 3);
		printMatrix(&I);
		
		printf("Initializing a 10x1 matrix, J.\n\n");
		zeroMatrix(&J, 10, 1);
		printMatrix(&J);
		
		printf("Initializing a 10x1 matrix, K.\n\n");
		randMatrix(&K, 10, 1);
		printMatrix(&K);
		
		//Large matrix declarations - root node

		printf("Initialzing a 3000x3000 matrix, L1.\n\n");
		randMatrix(&L1, 3000, 3000);
		
		printf("Initialzing a 3000x3000 matrix, L2.\n\n");
		randMatrix(&L2, 3000, 3000);
		
		printf("Initialzing a 4500x6783 matrix, L3.\n\n");
		randMatrix(&L3, 4500, 6783);
		
		printf("Initialzing a 150452x15987 matrix, L4.\n\n");
		randMatrix(&L4, 15452, 15987);
		
		printf("Initialzing a 15987x15765 matrix, L5.\n\n");
		randMatrix(&L5, 15987, 15765);
		
		printf("Initialzing a 3000x3000 matrix, L6.\n\n");
		randMatrix(&L6, 3000, 3000);

		printf("Initialzing a 15452x15765 matrix, L7.\n\n");
		randMatrix(&L7, 15452, 15765);

	}
	else
	{
		
		//Large vector, inner product declarations - non root nodes
		
		initNullMatrix(&v1, 20, 1);
		initNullMatrix(&v2, 20, 1);
		initNullMatrix(&v3, 500, 1);
		initNullMatrix(&v4, 500, 1);
		initNullMatrix(&v5, 3563, 1);
		initNullMatrix(&v6, 3563, 1);
		initNullMatrix(&v7, 7693, 1);
		initNullMatrix(&v8, 7693, 1);
		initNullMatrix(&v9, 15000, 1);
		initNullMatrix(&v10, 15000, 1);
		
		//General testing matrix declarations - non root nodes

		initNullMatrix(&A, 3, 3);
		initNullMatrix(&B, 3, 3);
		initNullMatrix(&C, 3, 3);
		initNullMatrix(&D, 3, 3);
		initNullMatrix(&E, 2, 4);
		initNullMatrix(&F, 4, 3);
		initNullMatrix(&G, 10, 10);
		initNullMatrix(&H, 10, 1);
		initNullMatrix(&I, 2, 3);
		initNullMatrix(&J, 10, 1);
		initNullMatrix(&K, 10, 1);

		//Large testing

		initNullMatrix(&L1, 3000, 3000);
		initNullMatrix(&L2, 3000, 3000);
		initNullMatrix(&L3, 4500, 6783);
		initNullMatrix(&L4, 15452, 15987);
		initNullMatrix(&L5, 15987, 15765);
		initNullMatrix(&L6, 3000, 3000);
		initNullMatrix(&L7, 15452, 15765);
	}
	MPI_Barrier(world);
	t1 = MPI_Wtime();
	C.data = add(&A,&B,world, worldSize, myRank);
	MPI_Barrier(world);
	t2 = MPI_Wtime();
	if (myRank == 0)
	{
		printf("Time taken A+B: %f seconds.\n\n", t2-t1);
	}
	
	MPI_Barrier(world);
	t1 = MPI_Wtime();
	D.data = subtract(&A, &B, world, worldSize, myRank);
	MPI_Barrier(world);
	t2 = MPI_Wtime();
	if (myRank == 0)
	{
		printf("Time taken A-B: %f seconds.\n\n", t2-t1);
	}
	
	MPI_Barrier(world);
	t1 = MPI_Wtime();
	I.data = multiply(&E, &F, world, worldSize, myRank);	
	MPI_Barrier(world);
	t2 = MPI_Wtime();
	if (myRank == 0)
	{
		printf("Time taken E*F: %f seconds.\n\n", t2-t1);
	}

	MPI_Barrier(world);
	t1 = MPI_Wtime();
	J.data = multiply(&G, &H, world, worldSize, myRank);
	MPI_Barrier(world);
	t2 = MPI_Wtime();
	if (myRank == 0)
	{
		printf("Time taken G*H: %f seconds.\n\n", t2-t1);
	}

	// Vector inner products
	
	MPI_Barrier(world);
	t1 = MPI_Wtime();
	double ip1 = innerProduct(&v1, &v2, world, worldSize, myRank);
	MPI_Barrier(world);
	t2 = MPI_Wtime();
	if (myRank == 0)
	{
		printf("Time taken v1*v2: %f seconds.\n\n", t2-t1);
	}

	MPI_Barrier(world);
	t1 = MPI_Wtime();
	double ip2 = innerProduct(&v3, &v4, world, worldSize, myRank);
	MPI_Barrier(world);
	t2 = MPI_Wtime();
	if (myRank == 0)
	{
		printf("Time taken v3*v4: %f seconds.\n\n", t2-t1);
	}

	MPI_Barrier(world);
	t1 = MPI_Wtime();
	double ip3 = innerProduct(&v5, &v6, world, worldSize, myRank);
	MPI_Barrier(world);
	t2 = MPI_Wtime();
	if (myRank == 0)
	{
		printf("Time taken v5*v6: %f seconds.\n\n", t2-t1);
	}

	MPI_Barrier(world);
	t1 = MPI_Wtime();
	double ip4 = innerProduct(&v7, &v8, world, worldSize, myRank);
	MPI_Barrier(world);
	t2 = MPI_Wtime();
	if (myRank == 0)
	{
		printf("Time taken v7*v8: %f seconds.\n\n", t2-t1);
	}

	MPI_Barrier(world);
	t1 = MPI_Wtime();
	double ip5 = innerProduct(&v9, &v10, world, worldSize, myRank);
	MPI_Barrier(world);
	t2 = MPI_Wtime();
	if (myRank == 0)
	{
		printf("Time taken v9*v10: %f seconds.\n\n", t2-t1);
	}

	if (myRank == 0)
	{
		printf("Testing the functionality of the inner product function.\n\n");
		printf("v1 * v2 = %0.2f\n", ip1);
		printf("v3 * v4 = %0.2f\n", ip2);
		printf("v5 * v6 = %0.2f\n", ip3);
		printf("v7 * v8 = %0.2f\n", ip4);
		printf("v9 * v10 = %0.2f\n\n\n", ip5);
	}		
	if (myRank == 0)
	{
		t1 = MPI_Wtime();
		E = transpose(&E);
		t2 = MPI_Wtime();
		printf("Time taken E^T: %f\n\n", t2-t1);

		printf("Testing operations.\n\n");

		printf("Matrix C = A + B:\n\n");
		printMatrix(&C);
		printf("\nMatrix D = A - B:\n\n");
		printMatrix(&D);
		printf("\nMatrix E = E Transpose:\n\n");
		printMatrix(&E);
		printf("\nMatrix I = E * F:\n\n");
		printMatrix(&I);
		printf("\nMatrix J = G * H:\n\n"); 
		printMatrix(&J);
	}
	if (myRank == 0)
	{
		printf("Testing large operations! :)\n");
		printf("These won't be printed out, but 'Done' will print when complete.\n\n");
	}

	MPI_Barrier(world);
	t1 = MPI_Wtime();
	L6.data = add(&L1, &L2, world, worldSize, myRank);
	MPI_Barrier(world);
	t2 = MPI_Wtime();
	if (myRank == 0)
	{
		printf("L1 + L2, Done.\n\n");
		printf("Time taken L1+L2: %f seconds.\n\n", t2-t1);
	}
	
	MPI_Barrier(world);
	t1 = MPI_Wtime();
	L6.data = subtract(&L1, &L2, world, worldSize, myRank);
	MPI_Barrier(world);
	t2 = MPI_Wtime();
	if (myRank == 0)
	{
		printf("L1 - L2, Done.\n\n");
		printf("Time taken L1-L2: %f seconds.\n\n", t2-t1);
	}
	
	if (myRank == 0)
	{
		t1 = MPI_Wtime();
		L3 = transpose(&L3);
		t2 = MPI_Wtime();
		printf("L3 transpose, Done.\n\n");
		printf("Time taken L3^T: %f seconds.\n\n", t2-t1);
	}

	MPI_Barrier(world);
	t1 = MPI_Wtime();
	//L7.data = multiply(&L4, &L5, world, worldSize, myRank);
	MPI_Barrier(world);
	t2 = MPI_Wtime();
	if (myRank == 0)
	{
		printf("L4 * L5, Done.\n\n");
		printf("Time taken L4*L5: %f seconds.\n\n", t2-t1);
	}
	


		
	free(A.data);
	free(B.data);
	free(C.data);
	free(D.data);
	free(E.data);
	free(F.data);
	free(G.data);
	free(H.data);
	free(I.data);
	free(J.data);
	free(K.data);
	free(v1.data);
	free(v2.data);
	free(v3.data);
	free(v4.data);
	free(v5.data);
	free(v6.data);
	free(v7.data);
	free(v8.data);
	free(v9.data);
	free(v10.data);
	free(L1.data);
	free(L2.data);
	free(L3.data);
	free(L4.data);
	free(L5.data);
	free(L6.data);
	free(L7.data);

	MPI_Finalize();
	return 0;
}
