#include <pthread.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>

pthread_barrier_t barrier;

//returns a square matrix of size d, populated by random numbers between 0 and 1
double** createMatrix(int d){
    //allocate memory for matrix
    double* values = malloc((long unsigned) d * (long unsigned) d * sizeof(double));
    double** matrix = malloc((long unsigned) d * sizeof(double*));
    int i;
    for (i = 0; i < d; i++){
        matrix[i] = values + i*d;
    }

    //populate with random numbers, dependent on the dimension
    srand((unsigned) d);
    int j;
    for (i = 0; i < d; i++){
        for (j = 0; j < d; j++){
            matrix[i][j] = ((rand() / (double) RAND_MAX));
        }
    }
    return matrix;
}


//returns a square matrix of size d, used for storage of values after averaging
double** createCopyMatrix(double** matrix, int d){
    //allocate memory for matrix
    double* copyValues = malloc((long unsigned) d * (long unsigned) d * sizeof(double));
    double** copyMatrix = malloc((long unsigned) d * sizeof(double*));
    int i, j;
    for (i = 0; i < d; i++){
        copyMatrix[i] = copyValues + i*d;
    }
    for (i = 0; i < d; i++){
        for (j = 0; j < d; j++){
            copyMatrix[i][j] = matrix[i][j];
        }
    }

    return copyMatrix;
}

void freeMatrices(double** matrix, double** copyMatrix){
    //free necessary memory
    free(*matrix);
    free(matrix);

    free(*copyMatrix);
    free(copyMatrix);
}

//sequential implementation
void sequential (double** matrix, double** copyMatrix, int d, double precision){
    double** tmp;
    bool done = false;

    //while the difference between the old and new values for all elements of the matrix are more than the precision
    while (done == false){
        int count = 0;

        //calculate average, and check whether the difference is less than the average. If it is, record it as an increment on count.
        int i, j;
        for (i = 1; i < d-1; i++){
            for (j = 1; j < d-1; j++){
                copyMatrix[i][j] = (matrix[i][j-1] + matrix[i][j+1] + matrix[i - 1][j] + matrix[i + 1][j]) / 4;
                if(fabs(copyMatrix[i][j] - matrix[i][j]) < precision){
                    count++;
                }
            }
        }

        //swap pointers to arrays
        tmp = matrix;
        matrix = copyMatrix;
        copyMatrix = tmp;

        //if the number of elements that have a difference less than the precision is equal to the total number of elements being processed, we're done
        if (count == ((d * d) - (4 * d) + 4)){
            done = true;
        }
    }
}

//contains range of rows to process
typedef struct Range {
    int start;
    int end;
}Range;

//arguments to be passed in to function run by each thread
typedef struct Arguments {
    int id;
    double** matrix;
    int dimension;
    double** copyMatrix;
    double precision;
    int nThreads;
    bool* allDone;
	//pthread_barrier_t *barrier;
    Range* ranges;
}Arguments;

//checks if all threads are done
bool ifAllDone(bool* a, int n){
    int i;
    for (i = 0; i < n; i++){
        if (a[i] == false){
            return false;
        }
    }
    return true;
}

void* avg(void* arguments){
    Arguments *args = (Arguments *) arguments;
    double** matrix = (*args).matrix;
    double precision = (*args).precision;
    int d = (*args).dimension;
    double** copyMatrix = (*args).copyMatrix;
    int threadID = (*args).id;
    bool* allDone = (*args).allDone;
    double** tmp;

    //while all threads have not reported they are done
    while (ifAllDone(allDone, (*args).nThreads) == false){

        int count = 0;
        int i, j;
        //calculate average, and check whether the difference is less than the average. If it is, record it as an increment on count.
        for (i = ((*args).ranges[threadID].start); i <= ((*args).ranges[threadID].end); i++){
            for (j = 1; j < d-1; j++){
                copyMatrix[i][j] = (matrix[i][j-1] + matrix[i][j+1] + matrix[i - 1][j] + matrix[i + 1][j]) / 4;
                if(fabs(copyMatrix[i][j] - matrix[i][j]) < precision){
                    count++;
                }
            }
        }

        //wait for all threads to finish prcossing their given rows
        pthread_barrier_wait(&barrier); 

        //swap pointers to arrays
        tmp = matrix;
        matrix = copyMatrix;
        copyMatrix = tmp;

        //if the number of elements that have a difference less than the precision is equal to the total number of elements being processed, this thread is done
        if (count == ((d - 2) * (1 + ((*args).ranges[threadID].end - (*args).ranges[threadID].start)))){
            allDone[threadID] = true;
        }

        //wait for all threads to report if they're done
        pthread_barrier_wait(&barrier);
    }
    return NULL;
}

void parallel(double** matrix, double** copyMatrix, int d, double precision, int n){

    //distributing rows
    Range range[n];
    int nRowsEach = (d - 2) / n;
    int leftover = d - 2 - (nRowsEach * n);
    int start = 1;
    int end = nRowsEach;
    int i;
    for (i = 0; i < n; i++){
        range[i].start = start;
        //give leftovers along the way, if there are any
        if (leftover > 0){
            end++;
            leftover--;
        }
        range[i].end = end;

        start = end + 1;
        end = start + nRowsEach - 1;
    }

    //contains all threads
    pthread_t threads[n];
    //initialise barrier
    pthread_barrier_init(&barrier, NULL, (unsigned) n);
    //contains whether all threads are done
    bool allDone[n];
    for(i = 0; i < n; i++){
        allDone[i] = 0;
    }
    struct Arguments threadData[n];

    for(i = 0; i < n; i++){
        threadData[i].id = i;
        threadData[i].matrix = matrix;
        threadData[i].copyMatrix = copyMatrix;
        threadData[i].dimension = d;
        threadData[i].precision = precision;
        threadData[i].ranges = range;
        threadData[i].nThreads = n;
        threadData[i].allDone = allDone;
    }

    //creation of threads
    for (i = 0; i < n; i++){
        pthread_create(&threads[i], NULL, &avg, (void *)&threadData[i]);
    }

    //join threads
    for(i = 0; i < n; i++){
        pthread_join(threads[i], NULL);
    }

    pthread_barrier_destroy(&barrier);
}


int main(void){
    clock_t start, end;
    start = clock();

    int dimension = 300;
    double precision = 0.001;
    int numberOfThreads = 16;
    double** matrix = createMatrix(dimension);
    double** copyMatrix = createCopyMatrix(matrix, dimension);
    //sequential(matrix, copyMatrix, dimension, precision);
    parallel(matrix, copyMatrix, dimension, precision, numberOfThreads);
    freeMatrices(matrix, copyMatrix);

    end = clock();
    printf("%f\n", (double) (end - start) / (CLOCKS_PER_SEC * numberOfThreads));

    return 0;
}