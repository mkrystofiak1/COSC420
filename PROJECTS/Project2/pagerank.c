#include "matrix.h"
int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    int worldSize, myRank;
    MPI_Comm world = MPI_COMM_WORLD;
    MPI_Comm_size(world, &worldSize);
    MPI_Comm_rank(world, &myRank);
    matrix M;
    M.rows=3;
    M.cols=3;
    M.data=malloc(sizeof(double)*M.rows*M.cols);
    matrix *AM=&M;
    create_random_matrix(&M,2,myRank,world);
    MPI_Barrier(world);
    if (myRank == 0)
    {
        printf("\nM=\n");
        for (int i = 0; i < M.rows; i++)
        {

            printf("[");
            for (int j = 0; j < M.cols; j++)
            {
                printf(" %lf ", ACCESS(AM,i,j));
            }

            printf("]\n");
        }
    }
    double * fin=NULL;
    fin=page_rank(&M, world, worldSize, myRank);
    MPI_Bcast( fin , M.rows , MPI_DOUBLE, 0, world);
   // MPI_Barrier(world);
    if (myRank == 0)
    {

        printf("\nFinal M page ranks=\n");
        for (int i = 0; i < M.rows; i++)
        {
            printf("\n[");
            printf(" %lf ", fin[i]);
            printf("]\n");
        }
    }
    free(M.data);
    free(AM);
    MPI_Finalize();
    return 0;
}