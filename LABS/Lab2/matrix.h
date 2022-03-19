//Program: matrix.h
//Author: Mitchell Krystofiak
//Description: Definition of the matrix class.
//Date: September 24, 2021

#ifndef MATRIX_H
#define MATRIX_H

#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<mpi.h>

#define INDEX(A,i,j) A->cols*i + j
#define ACCESS(A,i,j) A->data[INDEX(A,i,j)]

typedef struct {
	int rows;
	int cols;
	double* data;
} matrix;

//matrix operations are going to take in the MPI_Comm world variable
//and the int myRank variable to communicate with other nodes

void randMatrix(matrix *A, int rows, int cols);
void zeroMatrix(matrix *A, int rows, int cols);
void printMatrix(matrix *A);

double * add(matrix *A, matrix *B, MPI_Comm world, int worldSize, int myRank);
double * subtract(matrix *A, matrix *B, MPI_Comm world, int worldSize, int myRank);
double * multiply(matrix *A, matrix *B, MPI_Comm world, int worldsize, int myRank);
double innerProduct(matrix *A, matrix *B, MPI_Comm world, int worldSize, int myRank);
matrix transpose(matrix *A);


#endif
