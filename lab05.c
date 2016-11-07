#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#define n 10
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

void print_m(int matrix[n][n/p]) {
        int i, j;
        for(i = 0; i < n; i++) {
                for(j = 0; j < n/p; j++) {
                        printf("%d\t", matrix[i][j]);
                }
                printf("\n");
        }
}

void broadcast(int matrix[n][n], int rank) {
        int i, j;
        int np = n / p;
        int sendCounts[p], displs[p];
        int recvcount = n * (n / p);
        if(rank == 0) {
                for(i = 0, j = 0; i < p; i++) {
                        sendCounts[i] = np;
                        displs[i] = j;
                        j += np;
                }
        }
        // int MPI_Scatterv(void *sendbuf, int *sendcounts, int *displs, MPI_Datatype sendtype, void *recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm);
        MPI_Scatterv(matrix, sendCounts, displs, MPI_INT, m, recvcount, MPI_INT, 0, MPI_COMM_WORLD);
        if(rank == 0) {
                MPI_Scatterv(matrix, sendCounts, displs, MPI_INT, m, recvcount, MPI_INT, 0, MPI_COMM_WORLD);
        }
}

void column_sum(int matrix[n][n]) {
        int i, j;
}

int main(int argc, char** argv) {
        int rank, nprocs;
        int matrix[n][n];
        MPI_Init(&argc, &argv);
        MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);

        printf("Hello from processor %d of %d\n", rank, nprocs);

        srand(0);
        if(rank == 0) {
                int i, j;
                for(i = 0; i < n; i++) {
                        for(j = 0; j < n; j++) {
                                matrix[i][j] = (rand() % 10) + 1;
                        }
                }
                print_matrix(matrix);
                broadcast(matrix, rank);
        } else {
                broadcast(matrix, rank);
        }
        print_m(m);

        // column_sum(matrix);

        MPI_Finalize();
        return 0;
}
