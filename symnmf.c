#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define CHECK_MEMORY_ALLOCATION(ptr) \
    do { \
        if (!(ptr)) { \
            fprintf(stderr, "An error has occurred: Memory allocation failed.\n"); \
            exit(EXIT_FAILURE); \
        } \
    } while (0)


typedef struct {
    double *coordinates;
} point;

double** sym(point* points, int n, int d);
double euclideanDistance(double* p, double* q, int d);
void printMatrix(double** matrix, int n);




double euclideanDistance(double* p, double* q, int d) {
    double diff, sum = 0.0;
    int i;
    for (i = 0; i < d; i++) {
        diff = p[i] - q[i];
        sum += diff * diff;
    }
    return sqrt(sum);
}
double** sym(point* points, int n, int d) {
    // Allocate memory for the similarity matrix A
    double** A = (double**)malloc(n * sizeof(double*));
    CHECK_MEMORY_ALLOCATION(A);
    for (int i = 0; i < n; i++) {
        A[i] = (double*)malloc(n * sizeof(double));
        CHECK_MEMORY_ALLOCATION(A[i]);
    }

    // Calculate the similarity matrix A
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j) {
                A[i][j] = 0.0;  // Diagonal elements are 0
            } else {
                double distance = euclideanDistance(points[i].coordinates, points[j].coordinates, d);
                A[i][j] = exp(-distance * distance / 2.0);
            }
        }
    }

    return A;
}


void printMatrix(double** matrix, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%f ", matrix[i][j]);
        }
        printf("\n");
    }
}

