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

    // Calculate the normalized similarity matrix W
    matrix * mA = createMatrix(A, n, n);
    matrix * mD = createMatrix(D, n, n);
    diagPow(mD);
    matrix * tmp = matMul(mD, mA);
    matrix * W = matMul(tmp, mD);

    return W -> data;
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
void freeMatrix(double** matrix, int numPoints) {
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
double squaredFrobeniusNorm(matrix* m1, matrix* m2) {
    double sum = 0.0;
    for (int i = 0; i < m1->r; i++) {
        for (int j = 0; j < m1->c; j++) {
            double tmp = ((m1->data)[i][j] - (m2->data)[i][j]);
            sum += (tmp * tmp);
        }
    }
    return sum;
}

////////////////////////////////////////////////////////////////////////////////////////////////////

struct matrix{
    int r;
    int c;
    double ** data;
};
typedef struct matrix matrix;

double get(matrix* mat, int i, int j){
    return (mat->data)[i][j];
}

void set(matrix* mat, int i, int j, double val){
    (mat->data)[i][j] = val;
}

matrix * create_matrix(double ** arr, int r, int c){
    matrix * new = (matrix*)malloc(sizeof(matrix));
    new -> r = r;
    new -> c = c;
    new -> data = arr;
    return new;
}

double mulSum(double* arr1, double* arr2, int l){
   int i;
   int cnt = 0;
   for(i = 0; i < l; i++){
       cnt += arr1[i] * arr2[i];
   }
   return cnt;
}

double * getCol(matrix* mat, int index){
    int i;
    int len = mat->r;
    double * col = (double*)malloc(len * sizeof(double));
    for(i = 0; i < len; i++){
        col[i] = get(mat, i, index);
    }
    return col;
}

void printMat(matrix * m){
    int i,j;
    for(i = 0; i < m->r; i++){
        for(j = 0; j< m->c; j++){
            printf("%.4f \t", (m->data)[i][j]);
        }
        printf("\n");
    }
}

matrix* matMul(matrix* m1, matrix* m2){
    int i, j;
    double * col;
    int r1 = m1->r;
    int r2 = m2->r;
    int c1 = m1->c;
    int c2 = m2->c;
    
    if(c1 != r2){
        printf("sizes dont match");
    }
    double** ans = (double**)malloc(r1 * sizeof(double*));
    for (i = 0; i < r1; i++){
        ans[i] = (double*)malloc(c2 * sizeof(double));
    }
    
    for (i = 0; i < r1; i++){
        for (j = 0; j < c2; j++){
            col = getCol(m2, j);
            ans[i][j] = mulSum((m1->data)[i], col, r2);
        }
    }
    return create_matrix(ans, r1, c2);
}

void diagPow(matrix * m){
    //power of -0.5 to diagonal matrix
    int i, tmp;
    for(i = 0; i < m->c; i++){
        tmp = 1/sqrt((m->data)[i][i]);
        (m->data)[i][i] = tmp;
    }

}

matrix* transpose(matrix* m){
    int r = m->c;
    int c = m->r;
    double** arr = (double**)malloc(r * sizeof(double*));
    for (i = 0; i < r; i++){
        arr[i] = (double*)malloc(c * sizeof(double));
    }
     for (i = 0; i < r; i++)
        for (j = 0; j < c; j++)
            arr[i][j] = (m->data)[j][i];
    
    return create_matrix(arr, r, c);
    
}

matrix* oneIter(matrix * H, matrix * W){
    double b = 0.5;
    double tmp;
    int r = H->r;
    int c = H->c;

    matrix* mone = matMul(W,H);
    matrix* H_t = transpose(H);
    matrix* mehane = matMul(H, H_t);
    mehane = matMul(mehane, H);

    double ** top = mone -> data;
    double ** bottom = mehane -> data;
    
    double** arr = (double**)malloc(r * sizeof(double*));
    for (i = 0; i < r; i++){
        arr[i] = (double*)malloc(c * sizeof(double));
    }
     for (i = 0; i < r; i++){
        for (j = 0; j < c; j++){
            tmp = b * top[i][j] / bottom[i][j];
            arr[i][j] = 1 - b + tmp;
        }
     }
     return create_matrix(arr, r, c);

}

matrix * updateH(matrix *H, matrix *W. int iter, double epsilon){
    matrix* H_old = &H;
     matrix * H_new = NULL;
    do{
        H_new = oneIter(H_old, W);
        double fNorm =  squaredFrobeniusNorm(H_new, H_old);
        iter --;
        H_old = H_new;
    }
    while((sqrt(fNorm) >= epsilon*epsilon) || (iter > 0));

    return H_new;

}


