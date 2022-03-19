//Program: matrix.h
//Author: Mitchell Krystofiak
//Description: Definition of the matrix class.
//Date: October 31, 2021

#ifndef MATRIX_H
#define MATRIX_H

#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>
#include<mpi.h>
#include<unistd.h>
#include<string.h>

#define INDEX(A,i,j) A->cols*i + j
#define ACCESS(A,i,j) A->data[INDEX(A,i,j)]

#define MIN(x,y) (((x) < (y)) ? (x) : (y))
#define MAX(x,y) (((x) > (u)) ? (x) : (y))

typedef struct {
	int rows;
	int cols;
	double* data;
} matrix;

//matrix operations are going to take in the MPI_Comm world variable
//and the int myRank variable to communicate with other nodes

//Lab2

void randMatrix(matrix *A, int rows, int cols, int c);
void zeroMatrix(matrix *A, int rows, int cols);
void idenMatrix(matrix *A, int rows, int cols);
void printMatrix(matrix *A);
void copy(matrix *A, matrix *B);
void initNullMatrix(matrix *A, int rows, int cols);
int isEqual(matrix *A, matrix *B);

double * add(matrix *A, matrix *B, MPI_Comm world, int worldSize, int myRank);
double * subtract(matrix *A, matrix *B, MPI_Comm world, int worldSize, int myRank);
double * multiply(matrix *A, matrix *B, MPI_Comm world, int worldsize, int myRank);
double innerProduct(matrix *A, matrix *B, MPI_Comm world, int worldSize, int myRank);
double * transpose(matrix *A);

//Lab3

double * Gauss_Jordan(matrix *A, matrix *B, MPI_Comm * world, int worldSize, int myRank); 

//Lab4

double norm(matrix *A, MPI_Comm world, int worldSize, int myRank);
double * EigenVector(int rows, int cols, MPI_Comm world, int worldSize, int myRank);

#endif
