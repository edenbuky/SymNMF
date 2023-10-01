#ifndef symnmf_h
#define symnmf_h
#define ITER 300
#define EPSILON 0.0001

typedef struct {
    double *coordinates;
}point;
struct matrix{
    int r;
    int c;
    double ** data;
};
typedef struct matrix matrix;
point* readPointsFromFile(const char* filename, int* numPoints, int* dimensions);
double euclideanDistance(double* p, double* q, int d);
matrix * sym(point* points, int n, int d);
matrix * ddg(matrix* A);
matrix* norm(matrix* A, matrix* D);
void printMatrix(double** matrix, int n);
void printInvalidInputError(const char* message);
void freePoints(point* points, int numPoints);
void freeMatrix(matrix * m);
double squaredFrobeniusNorm(matrix* m1, matrix* m2);
matrix * updateH(matrix *H, matrix *W);
matrix* oneIter(matrix * H, matrix * W);
matrix* transpose(matrix* m);
void diagPow(matrix * m);
matrix* matMul(matrix* m1, matrix* m2);
void printMat(matrix * m);
double * getCol(matrix* mat, int index);
double mulSum(double* arr1, double* arr2, int l);
matrix * create_matrix(double ** arr, int r, int c);
void set(matrix* mat, int i, int j, double val);
double get(matrix* mat, int i, int j);
int checkMemoryAllocation(void* ptr);
double ** build2Darray(int r, int c);
void freeArray(double** matrix, int numPoints);

#endif
