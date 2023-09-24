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
