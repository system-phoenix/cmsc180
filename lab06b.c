#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define n 4008

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

void transform(int **matrix, int row, int col) {
    int * tMatrix = (int *) malloc(sizeof(int) * row * col);
    int i, j;
    for(i = 0; i < row; i++) {
        for(j = 0; j < col; j++) {
            if(i == j) {
                tMatrix[i * col + j] = 0;
            } else {
                tMatrix[i * col + j] = tMatrix[j * row + i] = (*matrix)[i * col + j];
            }
        }
    }
    free((*matrix));
    (*matrix) = tMatrix;
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

void matrix_multiply(int **result, int **A, int ay, int ax, int **B, int by, int bx) {
    if(ax == by) {
        (*result) = (int *) malloc(sizeof(int) * ay * bx);
        int i, j, k;
        for(i = 0; i < ay; i++) {
            for(j = 0; j < bx; j++) {
                int sum = 0;
                for(k = 0; k < ax; k++) {
                    sum += (*A)[i * bx + k] * (*B)[j * ay + k];
                }
                (*result)[i * bx + j] = sum;
            }
        }
    }
}

void matrix_divide(int ** A, int ay, int ax, int B) {
    int i, j;
    for(i = 0; i < ay; i++) {
        for(j = 0; j < ax; j++) {
            (*A)[i * ax + j] /= B;
        }
    }
}

long void_columns(int * A, int * B) {
    int i, j, squares = 0;
    for(i = 0; i < n; i++) {
        if((A)[i * n + i]) {
            for(j = 0; j < n; j++) {
                if(i == j) continue;
                else {
                    if((B)[i * n + j]) {
                        (A)[j * n + j] = 0;
                    }
                }
            }

            squares += (A)[i * n + i];
        }
    }
    return squares;
}

int main(int argc, char** argv) {
    int rank, p;
    int i = 0;
    int *absMatrix = NULL, *resMatrix = NULL, *revMatrix = NULL;
    int *A = NULL, *B = NULL;
    int *temp = NULL;
    long squares;

    clock_t begin, end;
    double time_spent;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    srand(time(NULL));

    // initialize the local matrices;
    A = (int *) malloc(sizeof(int) * (n * n / p));
    resMatrix = (int *) malloc(sizeof(int) * (n * n));
    revMatrix = (int *) malloc(sizeof(int) * (n * n));
    printf("Creating matrix...\n");
    absMatrix = (int *) malloc(sizeof(int) * (n * n));  // initialize the main matrix to be distributed
    
    for(i = 0; i < n * n; i++) {
        resMatrix[i] = 0;
    }

    if(rank == 0) {
        // if rank is 0, or the root,

        printf("Randomizing the contents of the matrix...\n");
        for(i = 0; i < n * n; i++) {
            absMatrix[i] = rand() % 2;                 //randomize value for the matrix
        }
        transform(&absMatrix, n, n);

        // absMatrix[0] = 0;   absMatrix[1] = 1;   absMatrix[2] = 0;   absMatrix[3] = 0;   absMatrix[4] = 1;

        // absMatrix[5] = 1;   absMatrix[6] = 0;   absMatrix[7] = 1;   absMatrix[8] = 1;   absMatrix[9] = 1;

        // absMatrix[10] = 0;  absMatrix[11] = 1;  absMatrix[12] = 0;  absMatrix[13] = 1;  absMatrix[14] = 0;

        // absMatrix[15] = 0;  absMatrix[16] = 1;  absMatrix[17] = 1;  absMatrix[18] = 0;  absMatrix[19] = 1;

        // absMatrix[20] = 1;  absMatrix[21] = 1;  absMatrix[22] = 0;  absMatrix[23] = 1;  absMatrix[24] = 0;


        // absMatrix[0] = 0;   absMatrix[1] = 1;   absMatrix[2] = 0;   absMatrix[3] = 1;

        // absMatrix[4] = 1;   absMatrix[5] = 0;   absMatrix[6] = 1;   absMatrix[7] = 1;

        // absMatrix[8] = 0;   absMatrix[9] = 1;   absMatrix[10] = 0;  absMatrix[11] = 1;

        // absMatrix[12] = 1;  absMatrix[13] = 1;  absMatrix[14] = 1;  absMatrix[15] = 0;
        // printf("absolute matrix\n");
        // printMatrix(absMatrix, n, n, 0);

        printf("Clock begun!\n");
        begin = clock();
    }

    MPI_Bcast(absMatrix, n * n, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatter(absMatrix, n * n / p, MPI_INT, A, n * n / p, MPI_INT, 0, MPI_COMM_WORLD);
    matrix_multiply(&resMatrix, &absMatrix, n, n, &A, n, n / p);
    MPI_Gather(resMatrix, n * n / p, MPI_INT, revMatrix, n * n / p, MPI_INT, 0, MPI_COMM_WORLD);

    MPI_Bcast(revMatrix, n * n, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatter(revMatrix, n * n / p, MPI_INT, A, n * n / p, MPI_INT, 0, MPI_COMM_WORLD);
    matrix_multiply(&resMatrix, &revMatrix, n, n, &A, n, n / p);
    MPI_Gather(resMatrix, n * n / p, MPI_INT, revMatrix, n * n / p, MPI_INT, 0, MPI_COMM_WORLD);
    
    if(rank == 0) {
        // printf("final matrix\n");
        // printMatrix(revMatrix, n, n, 0);

        matrix_divide(&revMatrix, n, n, 8);
        squares = void_columns(revMatrix, absMatrix);

        end = clock();
        printf("Process done!\n");
        time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

        printf("Squares: %ld\n", squares);

        printf("n: %d with p: %d\n", n, p);
        printf("Time elapsed: %lf\n", time_spent);
    }

    MPI_Finalize();
    return 0;
}
