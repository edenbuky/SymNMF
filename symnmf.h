#ifndef symnmf_h
#define symnmf_h

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
double** sym(point* points, int n, int d);
double** ddg(double** A, int n);
double** norm(double** A, double** D, int n);
void printMatrix(double** matrix, int n);
void printInvalidInputError(const char* message);
void freePoints(point* points, int numPoints);
void freeMatrix(double** matrix, int numPoints);
double squaredFrobeniusNorm(double** matrix, int m, int n);
matrix * updateH(matrix *H, matrix *W. int iter, double epsilon);
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


#endif
