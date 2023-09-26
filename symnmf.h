#ifndef symnmf.h
#define symnmf.h
double euclideanDistance(double* p, double* q, int d);
double** sym(point* points, int n, int d);
double** ddg(double** A, int n);
double** norm(double** A, double** D, int n);
static PyObject *fit(PyObject *self, PyObject *args);
typedef struct {
    double *coordinates;
} point;
#endif //symnmf.h
#define CHECK_DOUBLE_POINTER_MEMORY_ALLOCATION(ptr, size) \
    do { \
        if (!(ptr)) { \
            fprintf(stderr, "An error has occurred: Memory allocation for double pointer failed.\n"); \
            exit(EXIT_FAILURE); \
        } \
        for (int i = 0; i < (size); i++) { \
            if (!(ptr)[i]) { \
                fprintf(stderr, "An error has occurred: Memory allocation for inner pointer failed.\n"); \
                exit(EXIT_FAILURE); \
            } \
        } \
    } while (0)
