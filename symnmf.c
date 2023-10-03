#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "symnmf.h"



int main(int argc, char* argv[]) {
    const char* goal;
    const char* filename;
    int numPoints, dimensions;
    point* points;
    matrix* mat;
    matrix* A;

    if (argc != 3) {
        fprintf(stderr, "Usage: %s <goal> <filename>\n", argv[0]);
        return EXIT_FAILURE;
    }

    goal = argv[1];
    filename = argv[2];

    points = readPointsFromFile(filename, &numPoints, &dimensions);
    if (!points) {
        fprintf(stderr, "An Error Has Occurred\n");
        return EXIT_FAILURE;
    }

    mat = NULL;

    if (strcmp(goal, "sym") == 0) {
        mat = sym(points, numPoints, dimensions);
    } else if (strcmp(goal, "ddg") == 0) {
        A = sym(points, numPoints, dimensions);
        mat = ddg(A);
        /* Free the memory for A after using it */
        freeMatrix(A);
    } else if (strcmp(goal, "norm") == 0) {
        matrix * A = sym(points, numPoints, dimensions);
        matrix* D = ddg(A);
        mat = norm(A, D);
        /* Free the memory for A and D after using them */
        freeMatrix(A);
        freeMatrix(D);
    } else {
        fprintf(stderr, "Invalid goal specified\n");
        freePoints(points, numPoints);
        return EXIT_FAILURE;
    }

    /* Print the resulting matrix */
    printMat(mat);

    /* Free the memory for the matrix and points */
    freeMatrix(mat);
    freePoints(points, numPoints);

    return EXIT_SUCCESS;
}

int checkMemoryAllocation(void* ptr) {
    if (!ptr) {
        fprintf(stderr, "An error has occurred: Memory allocation failed.\n");
        return 0;
    } return 1;
}

double ** build2Darray(int r, int c){
    int i,j;
    double** A = (double**)malloc(r * sizeof(double*));
    if(checkMemoryAllocation(A) == 0){
        return NULL;
    }
    for (i = 0; i < r; i++) {
        A[i] = (double*)malloc(c * sizeof(double));
        if(checkMemoryAllocation(A[i]) == 0){
            for(j = 0; j < i; j++){
                free(A[j]);
            }
            free(A);
            return NULL;
        }
    }
    return A;
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
matrix * sym(point* points, int n, int d) {
    /* Allocate memory for the similarity matrix A*/
    int i,j;
    double** A = build2Darray(n,n);

    /* Calculate the similarity matrix A */
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (i == j) {
                A[i][j] = 0.0;  /*Diagonal elements are 0*/
            } else {
                double distance = euclideanDistance(points[i].coordinates, points[j].coordinates, d);
                A[i][j] = exp(-distance * distance / 2.0);
            }
        }
    }

    return create_matrix(A, n, n);
}

matrix * ddg(matrix* A) {
    int i,j,n;
    double ** D;
    double sum;
    n = A->c;
    D = build2Darray(n,n);
    if(D == NULL){
        freeMatrix(A);
        exit(1);
    }
    for (i = 0; i < n; i++) {
         sum = 0.0;
        for (j = 0; j < n; j++) {
            sum += (A->data)[i][j];
        }
        D[i][i] = sum;  
        for(j=0; j<n; j++){
            if (j != i){
                D[i][j] = 0.0;
            }
        }
    }

    return create_matrix(D, n, n);
}

matrix* norm(matrix* A, matrix* D) {
    matrix * tmp;
    matrix * W;
    diagPow(D);
    tmp = matMul(D, A);
    W = matMul(tmp, D);

    return W;
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
        printf("An Error Has Occurred\n");
        exit(1);
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
            printf("An Error Has Occurred\n");
            exit(1);
        }
    }
    
    points = (point*)malloc(*numPoints * sizeof(point));
    if(checkMemoryAllocation(points) == 0){
        freePoints(points, *numPoints);
        exit(1);

    }

    rewind(file);

    pointIndex = 0;
    while (fgets(line, sizeof(line), file)) {
        points[pointIndex].coordinates = (double*)malloc(*dimensions * sizeof(double));
        if (points[pointIndex].coordinates == NULL) {
            freePoints(points, pointIndex);
            printf("An Error Has Occurred\n");
            exit(1);
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
void freeArray(double** matrix, int numPoints) {
    int i;
    if (matrix == NULL) {
        return;
    }
    for (i = 0; i < numPoints; i++) {
        free(matrix[i]);
    }
    free(matrix);
}

void freeMatrix(matrix * m){
    freeArray(m->data, m->r);
    free(m);
}


double squaredFrobeniusNorm(matrix* m1, matrix* m2) {
    double sum = 0.0;
    int i,j;
    for (i = 0; i < m1->r; i++) {
        for (j = 0; j < m1->c; j++) {
            double tmp = ((m1->data)[i][j] - (m2->data)[i][j]);
            sum += (tmp * tmp);
        }
    }
    return sum;
}

double get(matrix* mat, int i, int j){
    return (mat->data)[i][j];
}

void set(matrix* mat, int i, int j, double val){
    (mat->data)[i][j] = val;
}

matrix * create_matrix(double ** arr, int r, int c){
    matrix * new = (matrix*)malloc(sizeof(matrix));
    if(checkMemoryAllocation(new) == 0){
        freeArray(arr, r);
        exit(1);
    }
    new -> r = r;
    new -> c = c;
    new -> data = arr;
    return new;
}

double mulSum(double* arr1, double* arr2, int l){
   int i;
   double cnt = 0.0;
   for(i = 0; i < l; i++){
       cnt += arr1[i] * arr2[i];
   }
   return cnt;
}

double * getCol(matrix* mat, int index){
    int i;
    int len = mat->r;
    double * col = (double*)malloc(len * sizeof(double));
    if(checkMemoryAllocation(col)==0){
        freeMatrix(mat);
        exit(1);
    }
    for(i = 0; i < len; i++){
        col[i] = get(mat, i, index);
    }
    return col;
}

void printMat(matrix * m){
    int i,j,n;
    n =m->c;
    for(i = 0; i < m->r; i++){
        for(j = 0; j< n-1; j++){
            printf("%.4f,", (m->data)[i][j]);
        }
        printf("%.4f\n", (m->data)[i][n-1]);
    }
}

/*matrix* matMul2(matrix* m1, matrix* m2){
    int i, j;
    double * col;
    double** ans;
    int r1 = m1->r;
    int r2 = m2->r;
    int c1 = m1->c;
    int c2 = m2->c;
    
    if(c1 != r2){
        printf("sizes dont match");
    }
    ans = build2Darray(r1,c2);
    if(!ans){
        freeMatrix(m1);
        freeMatrix(m2);
        return NULL;
    }
    
    for (i = 0; i < r1; i++){
        for (j = 0; j < c2; j++){
            col = getCol(m2, j);
            ans[i][j] = mulSum((m1->data)[i], col, r2);
        }
    }
    return create_matrix(ans, r1, c2);
}
*/

matrix* matMul(matrix* m1, matrix* m2){
    int i, j;
    double * col;
    double** ans;
    int r1 = m1->r;
    int r2 = m2->r;
    int c1 = m1->c;
    int c2 = m2->c;
  
    ans = build2Darray(r1,c2);
    if(!ans){
        freeMatrix(m1);
        freeMatrix(m2);
        return NULL;
    }
  
    for (int i = 0; i < r1; i++) {
        for (int j = 0; j < c2; j++) {
            ans[i][j] = 0;
  
            for (int k = 0; k < r2; k++) {
                ans[i][j] += (m1->data)[i][k] * (m2->data)[k][j];
            }
  
        }
    }
}


void diagPow(matrix * m){
    int i;
    double x, tmp;
    for(i = 0; i < m->c; i++){
        x = get(m, i, i);
        tmp = 1/sqrt(x);
        set(m, i, i, tmp);
    }

}

matrix* transpose(matrix* m){
    int i,j;
    int r = m->c;
    int c = m->r;
    double** arr = build2Darray(r,c);
    if(!arr){
        freeMatrix(m);
        return NULL;
    }
     for (i = 0; i < r; i++)
        for (j = 0; j < c; j++)
            arr[i][j] = (m->data)[j][i];
    
    return create_matrix(arr, r, c);
    
}

matrix* oneIter(matrix * H, matrix * W){
    int i,j;
    double b = 0.5;
    double tmp;
    matrix* mone;
    matrix* H_t;
    matrix* mehane;
    double ** bottom;
    double ** top;
    double ** arr;
    int r = H->r;
    int c = H->c;
    
    mone = matMul(W,H);
    H_t = transpose(H);
    if(!H_t){
        freeMatrix(W);
        exit(1);
    }
    mehane = matMul(H, H_t);
    if(!mehane){
        freeMatrix(W);
        exit(1);
    }
    mehane = matMul(mehane, H);
    if(!mehane){
        freeMatrix(W);
        exit(1);
    }

    top = mone -> data;
    bottom = mehane -> data;
    
    arr = build2Darray(r,c);
    if(!arr){
        freeMatrix(H);
        freeMatrix(W);
        freeMatrix(H_t);
        exit(1);
    }
     for (i = 0; i < r; i++){
        for (j = 0; j < c; j++){
            tmp = b * top[i][j] / bottom[i][j];
            arr[i][j] = (H->data)[i][j] * (1 - b + tmp);
        }
     }
     return create_matrix(arr, r, c);

}

matrix * updateH(matrix *H, matrix *W){
    double fNorm;
    int iter = ITER;
    matrix* H_old = H;
     matrix * H_new = NULL;
    do{
        H_new = oneIter(H_old, W);
        fNorm =  squaredFrobeniusNorm(H_new, H_old);
        iter --;
        H_old = H_new;
    }
    while((sqrt(fNorm) >= EPSILON) || (iter > 0));

    return H_new;

}


