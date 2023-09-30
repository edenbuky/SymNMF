#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "symnmf.h"

#define CHECK_MEMORY_ALLOCATION(ptr) \
    do { \
        if (!(ptr)) { \
            fprintf(stderr, "An error has occurred: Memory allocation failed.\n"); \
            exit(EXIT_FAILURE); \
        } \
    } while (0)

int main(int argc, char* argv[]) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <goal> <filename>\n", argv[0]);
        return EXIT_FAILURE;
    }

    const char* goal = argv[1];
    const char* filename = argv[2];

    int numPoints, dimensions;
    point* points = readPointsFromFile(filename, &numPoints, &dimensions);
    if (!points) {
        fprintf(stderr, "An Error Has Occurred\n");
        return EXIT_FAILURE;
    }

    double** matrix = NULL;

    if (strcmp(goal, "sym") == 0) {
        matrix = sym(points, numPoints, dimensions);
    } else if (strcmp(goal, "ddg") == 0) {
        double** A = sym(points, numPoints, dimensions);
        matrix = ddg(A, numPoints);
        // Free the memory for A after using it
        for (int i = 0; i < numPoints; i++) {
            free(A[i]);
        }
        free(A);
    } else if (strcmp(goal, "norm") == 0) {
        double** A = sym(points, numPoints, dimensions);
        double** D = ddg(A, numPoints);
        matrix = norm(A, D, numPoints);
        // Free the memory for A and D after using them
        for (int i = 0; i < numPoints; i++) {
            free(A[i]);
            free(D[i]);
        }
        free(A);
        free(D);
    } else {
        fprintf(stderr, "Invalid goal specified\n");
        freePoints(points, numPoints);
        return EXIT_FAILURE;
    }

    // Print the resulting matrix
    for (int i = 0; i < numPoints; i++) {
        for (int j = 0; j < numPoints; j++) {
            printf("%.4f", matrix[i][j]);
            if (j < numPoints - 1) {
                printf(",");
            }
        }
        printf("\n");
    }

    // Free the memory for the matrix and points
    freeMatrix(matrix,numPoints);
    freePoints(points, numPoints);

    return EXIT_SUCCESS;
}
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

double** ddg(double** A, int n) {
    // Allocate memory for the diagonal degree matrix D
    double** D = (double**)malloc(n * sizeof(double*));
    CHECK_MEMORY_ALLOCATION(D);
    for (int i = 0; i < n; i++) {
        D[i] = (double*)calloc(n, sizeof(double)); // Initialize with zeros
        CHECK_MEMORY_ALLOCATION(D[i]);
    }

    // Calculate the diagonal degree matrix D
    for (int i = 0; i < n; i++) {
        double sum = 0.0;
        for (int j = 0; j < n; j++) {
            sum += A[i][j];
        }
        D[i][i] = sum;  // Set the diagonal element
    }

    return D;
}

double** norm(double** A, double** D, int n) {
    // Allocate memory for the normalized similarity matrix W
    double** W = (double**)malloc(n * sizeof(double*));
    CHECK_MEMORY_ALLOCATION(W);
    for (int i = 0; i < n; i++) {
        W[i] = (double*)malloc(n * sizeof(double));
        CHECK_MEMORY_ALLOCATION(W[i]);
    }

    // Calculate the normalized similarity matrix W
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            W[i][j] = A[i][j] / (sqrt(D[i][i]) * sqrt(D[j][j]));
        }
    }

    return W;
}


void printMatrix(double** matrix, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n-1; j++) {
            printf("%f,", matrix[i][j]);
        }
        printf("%f\n", matrix[i][n]);
    }
}

point* readPointsFromFile(const char* filename, int* numPoints, int* dimensions) {
    char line[256];
    FILE* file;
    int count, pointIndex, coordinateIndex;
    char* token;
    point* points;
    *numPoints = 0;
    *dimensions = 0;
    file = fopen(filename, "r");
    if (file == NULL) {
        printInvalidInputError("An Error Has Occurred\n");
        return NULL;
    }
    while (fgets(line, sizeof(line), file)) {
        (*numPoints)++;
        count = 0;
        token = strtok(line, ",");
        while (token != NULL) {
            count++;
            token = strtok(NULL, ",");
        }
        if (*dimensions == 0) {
            *dimensions = count;
        } else if (*dimensions != count) {
            fclose(file);
            printInvalidInputError("An Error Has Occurred\n");
            return NULL;
        }
    }
    
    points = (point*)malloc(*numPoints * sizeof(point));
    CHECK_MEMORY_ALLOCATION(points);

    rewind(file);

    pointIndex = 0;
    while (fgets(line, sizeof(line), file)) {
        points[pointIndex].coordinates = (double*)malloc(*dimensions * sizeof(double));
        if (points[pointIndex].coordinates == NULL) {
            freePoints(points, pointIndex);
            printInvalidInputError("An Error Has Occurred\n");
            return NULL;
        }

        coordinateIndex = 0;
        token = strtok(line, ",");
        while (token != NULL) {
            points[pointIndex].coordinates[coordinateIndex] = atof(token);
            coordinateIndex++;
            token = strtok(NULL, ",");
        }
        ++pointIndex;
    }
    fclose(file);
    return points;
}
void freePoints(point* points, int numPoints) {
    int i;
    if (points == NULL) {
        return;
    }
    for (i = 0; i < numPoints; i++) {
        if (points[i].coordinates != NULL){
            free(points[i].coordinates);
        }
    }

    free(points);
}
void freeMatrix(ouble** matrix, int numPoints) {
    int i;
    if (matrix == NULL) {
        return;
    }
    for (i = 0; i < numPoints; i++) {
        free(matrix[i]);
    }
    free(matrix);
}
void printInvalidInputError(const char* message) {
    printf("%s\n", message);
    exit(1);
}
double squaredFrobeniusNorm(double** matrix, int m, int n) {
    double sum = 0.0;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            sum += matrix[i][j] * matrix[i][j];
        }
    }
    return sum;
}

