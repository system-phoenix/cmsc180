#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <time.h>
#define n 102

typedef struct node{
    int root, target, data;
    struct node * next;
} Path;

int rank, p;
int sMin, sMax; float pMin, pMax, gAvg;
int * absMatrix = NULL, * displc = NULL;
Path * paths = NULL;

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

void printCounts() {
    Path *ptr;
    printf("All of the data paths found:\n");
    for(ptr = paths; ptr; ptr = ptr->next) {
        printf("(%d, %d): %d\n", ptr->root, ptr->target, ptr->data);
    }
}

void transform(int **matrix, int row, int col) {
    int * tMatrix = (int *) malloc(sizeof(int) * row * col);
    int i, j;
    for(i = 0; i < row; i++) {
        for(j = 0; j < col; j++) {
            if(j < i) {
                continue;
            } else if(i == j) {
                tMatrix[i * col + j] = 0;
            } else if(i == 0) {
                tMatrix[i * col + j] = tMatrix[j * row + i] = 1;
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

void pathAdd(int root, int target, int data) {
    Path *path, *ptr;
    path = malloc(sizeof(Path));
    path->root = root;
    path->target = target;
    path->data = data;
    path->next = NULL;
    if(paths == NULL) {
        paths = path;
    } else {
        for(ptr = paths; ptr->next; ptr = ptr->next) {
            if(ptr->root == path->root && ptr->target == path->target && ptr->data == path->data) {
                free(path);
                path = NULL;
                break;
            }
        }
        if(path) {
            ptr->next = path;
            // printf("Successfully added path!\n");
        } else {
            // printf("Duplicate!\n");
        }
    }
}

void removeRedundancies() {
    Path *ptr, *stk;
    for(stk = paths, ptr = paths->next; ptr;) {
        if(ptr->root == stk->root && ptr->target == stk->target && ptr->data == stk->data) {
            stk->next = ptr->next;
            free(ptr);
            ptr = stk->next;
            // printf("Removed redundant data!\n");
            continue;
        }
        stk = stk->next; ptr = ptr->next;
    }
}

void DFS(int root, int target, int i, int path, int * visited) {
    int j;
    // printf("visited nodes so far... ");
    // for(j = 0; j < n; j++) {
    //     if(visited[j]) {
    //         printf("%d\t", j);
    //     }
    // }
    // printf("\n");
    if(i == target) {
        // printf("Path found! Length: %d\n", path);
        pathAdd(root, target, path);
        return;
    }
    for(j = 0; j < n; j++) {
        if(j == root || j == i || visited[j]) continue;
        else {
            if(absMatrix[i * n + j]) {
                visited[j] = 1;
                DFS(root, target, j, path + 1, visited);
                visited[j] = 0;
            }
        }
    } 
}

void preDFS() {
    int i, j, k;
    int *visited = NULL;
    // printf("Setting up...\n");
    for(i = displc[rank * 2]; i < displc[rank * 2 + 1]; i++) {
        // printf("Current root: %d\n", i);
        for(j = 0; j < n; j++) {
            if(j <= i) continue;
            else {
                // printf("\n\nFinding path from %d to %d...\n", i, j);
                visited = malloc(sizeof(int) * n);
                for(k = 0; k < n; k++) {
                    visited[k] = 0;
                }
                DFS(i, j, i, 0, visited);
                free(visited);
            }
        }
    }
}

void postDFS() {
    Path *ptr;
    for(; paths != NULL;) {
        ptr = paths;
        paths = paths->next;
        free(ptr);
    }
}

void computeSeries() {
    Path *ptr;
    int root = 0, target = 1, sum = 0;
    int min = -1, max = -1, cnt = 0;
    float avg = 0;
    for(ptr = paths; ptr;) {
        if(ptr->root != root || ptr->target != target) {
            root = ptr->root;
            target = ptr->target;
            if(min < 0 || sum < min) min = sum;
            if(max < 0 || sum > max) max = sum;
            sum = 0;
        } else {
            avg += ptr->data;
            sum += ptr->data;
            ptr = ptr->next;
            cnt++;
        }
    }
    sMin = min;
    sMax = max;
    avg /= (float)cnt;
    gAvg = avg;
    // printf("Series -- min: %d, max: %d, avg: %4f\n", min, max, avg);
}

void computeParallel() {
    Path *ptr;
    int root = 0, target = 1, cnt = 0;
    float min = -1, max = -1, avg = 0, sum = 0;
    for(ptr = paths; ptr;) {
        if(ptr->root != root || ptr->target != target) {
            sum = 1 / sum;
            root = ptr->root;
            target = ptr->target;
            if(min < 0 || sum < min) min = sum;
            if(max < 0 || sum > max) max = sum;
            sum = 0;
        } else {
            avg += ptr->data;
            sum += 1 / (float) ptr->data;
            ptr = ptr->next;
            cnt++;
        }
    }
    pMin = min;
    pMax = max;
    avg /= (float)cnt;
    if(gAvg != avg) {
        gAvg = (gAvg + avg) / 2;
    }
    // printf("Parallel -- min: %f, max: %f, avg: %4f\n", min, max, avg);
}

int main(int argc, char** argv) {
    int i = 0;
    int * seriesMin = NULL, *seriesMax = NULL;
    float * parallelMin = NULL, *parallelMax = NULL, *averages = NULL;

    clock_t begin, end;
    double time_spent;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    srand(time(NULL));

    printf("n: %d, p: %d\n", n, p);

    //initialize the main matrix
    absMatrix = (int *) malloc(sizeof(int) * (n * n));
    
    //configure the displacements
    displc = (int *) malloc(sizeof(int) * p * 2);
    for(i = 0; i < p; i++) {
        displc[i * 2]        = i * n / p;
        displc[i * 2 + 1]    = i * n / p + n / p - 1;
    }


    if(rank == 0) {
        // if rank is 0, or the root,

        printf("Randomizing the contents of the matrix...\n");
        for(i = 0; i < n * n; i++) {
            absMatrix[i] = rand() % 2;                 //randomize value for the matrix
        }
        transform(&absMatrix, n, n);

        seriesMin = (int *) malloc(sizeof(int) * p);
        seriesMax = (int *) malloc(sizeof(int) * p);
        parallelMin = (float *) malloc(sizeof(float) * p);
        parallelMax = (float *) malloc(sizeof(float) * p);
        averages = (float *) malloc(sizeof(float) * p);

        printf("Process begun!\n");
        begin = clock();
    }

    MPI_Bcast(absMatrix, n * n, MPI_INT, 0, MPI_COMM_WORLD);

    preDFS();
    removeRedundancies();
    computeSeries();
    computeParallel();

    MPI_Gather(&sMin, 1, MPI_INT, seriesMin, p, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather(&sMax, 1, MPI_INT, seriesMax, p, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather(&pMin, 1, MPI_INT, parallelMin, p, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather(&pMax, 1, MPI_INT, parallelMax, p, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gather(&gAvg, 1, MPI_INT, averages, p, MPI_INT, 0, MPI_COMM_WORLD);
    
    if(rank == 0) {
        sMin = sMax = -1; pMin = pMax = -1; gAvg = 0;
        for(i = 0; i < p; i++) {
            gAvg += averages[i];
            if(sMin < 0 || seriesMin[i] < sMin) sMin = seriesMin[i];
            if(sMax < 0 || seriesMax[i] > sMax) sMax = seriesMax[i];
            if(pMin < 0 || parallelMin[i] < pMin) pMin = parallelMin[i];
            if(pMax < 0 || parallelMax[i] > pMax) pMax = parallelMax[i];
        }
        end = clock();
        printf("Process done!\n");
        time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
        // printCounts();

        printf("n: %d with p: %d\n", n, p);
        printf("Time elapsed: %lf\n", time_spent);
    }

    printf("Cleaning up...\n");
    postDFS();
    free(absMatrix);
    free(displc);
    free(seriesMax);
    free(seriesMin);
    free(parallelMax);
    free(parallelMin);
    free(averages);
    MPI_Finalize();
    return 0;
}
