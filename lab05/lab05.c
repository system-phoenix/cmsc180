#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#define n 10000
#define p 1

int m[n][n /p];

void print_matrix(int matrix[n][n]) {
        int i, j;
        for(i = 0; i < n; i++) {
                for(j = 0; j < n; j++) {
                        printf("%d\t", matrix[i][j]);
                }
                printf("\n");
        }
}

void broadcast(int matrix[n][n]) {
        int i, j;
        int np = n / p;
}

void column_sum(int matrix[n][n]) {
        int i, j;
}

int main(int argc, char** argv) {
        int myrank, nprocs;
        int matrix[n][n];
        MPI_Init(&argc, &argv);
        MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
        MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

        printf("Hello from processor %d of %d\n", myrank, nprocs);

        srand(0);
        if(myrank == 0) {
                int i, j;
                for(i = 0; i < n; i++) {
                        for(j = 0; j < n; j++) {
                                matrix[i][j] = (rand() % 10) + 1;
                        }
                }
        }

        column_sum(matrix);

        MPI_Finalize();
        return 0;
}
