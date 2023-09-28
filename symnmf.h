#ifndef symnmf_h
#define symnmf_h

typedef struct {
    double *coordinates;
}point;

point* readPointsFromFile(const char* filename, int* numPoints, int* dimensions);
double euclideanDistance(double* p, double* q, int d);
double** sym(point* points, int n, int d);
double** ddg(double** A, int n);
double** norm(double** A, double** D, int n);
void printMatrix(double** matrix, int n);
void printInvalidInputError(const char* message);
void freePoints(point* points, int numPoints);
//static PyObject *fit(PyObject *self, PyObject *args);

#endif
