#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define n 10000

void printMatrix(int *matrix, int row, int col, int last) {
    printf("==================================\n");
    int i, j;
    for(i = 0; i < row; i++) {
        for(j = 0; j < col; j++) {
            if(!last) {
                printf("%d\t", matrix[i * col + j]);
            } else {
                if(i + 1 == row) {
                    printf("%d\t", matrix[i * col + j]);
                }
            }
        }
        if(!last) {
            printf("\n");
        } else {
            if(i + 1 == row) {
                printf("\n");
            }
        }
    }
    printf("==================================\n");
}

void printCounts(int * array, int proc) {
    int i;
    for(i = 0; i < proc; i++) {
        printf("%d\t", array[i]);
    }
    printf("\n");
}

void transpose(int **matrix, int row, int col) {
    int * tMatrix = (int *) malloc(sizeof(int) * row * col);
    int i, j;
    for(i = 0; i < row; i++) {
        for(j = 0; j < col; j++) {
            tMatrix[j * row + i] = (*matrix)[i * col + j];
        }
    }
    free((*matrix));
    (*matrix) = tMatrix;
}

int main(int argc, char** argv) {
    int rank, p;
    int i = 0;
    int displ = 0;
    int *absMatrix = NULL;
    int *relMatrix = NULL;
    int *sendCounts = NULL;
    int *displs = NULL;

    clock_t begin, end;
    double time_spent;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    srand(time(NULL));

    // initialize the local matrix;
    relMatrix = (int *) malloc(sizeof(int) * (n * n / p));

    if(rank == 0) {
        // if rank is 0, or the root,
        absMatrix = (int *) malloc(sizeof(int) * (n * n));  // initialize the main matrix to be distributed
        sendCounts = (int *) malloc(sizeof(int) * p);       // initialize the container for different sending sizes
        displs = (int *) malloc(sizeof(int) * p);           // initialize the displacement values

        for(i = 0; i < n * n; i++) {
            absMatrix[i] = rand() % 10 + 1;                 //randomize value for the matrix
        }
        // transpose(&absMatrix, n, n);

        for(i = 0; i < p; i++) {
            displs[i] = displ;
            sendCounts[i] = round((float) (n / p) + i * (n / p)) - round((float) i * (n / p));
            displ += sendCounts[i];
        }
    }

    begin = clock();

    MPI_Scatterv(absMatrix, sendCounts, displs, MPI_INT, relMatrix, n * n / p, MPI_INT, 0, MPI_COMM_WORLD);

    for(i = n / p; i < n * n / p; i++) {
        relMatrix[i] = relMatrix[i] + relMatrix[i - n / p];
    }

    MPI_Gatherv(relMatrix, n / p, MPI_INT, absMatrix, sendCounts, displs, MPI_INT, 0, MPI_COMM_WORLD);

    end = clock();
    time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

    printf("n: %d with p: %d\nTime elapsed: %lf\n", n, p, time_spent);

    MPI_Finalize();
    return 0;
}
