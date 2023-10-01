#define PY_SSIZE_T_CLEAN 
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include "symnmf.h"


point* convertPyListToPoints(PyObject *pyList, int n, int d) {
    point* points = (point*)malloc(n * sizeof(point));

    /*Iterate over the outer list*/
    for (int i = 0; i < n; i++) {
        PyObject *innerList = PyList_GetItem(pyList, i);

        /*Allocate memory for the coordinates of the point*/
        points[i].coordinates = (double*)malloc(d * sizeof(double));

        /*Iterate over the inner list*/
        for (int j = 0; j < d; j++) {
            PyObject *coordObj = PyList_GetItem(innerList, j);
            points[i].coordinates[j] = PyFloat_AsDouble(coordObj);
        }
    }

    return points;
}
point* convertPyListToArray(PyObject *pyList, int numRows, int numCols) {
    *numRows = PyList_Size(pyList);
    if (*numRows <= 0) {
        return NULL;
    }

    *numCols = PyList_Size(PyList_GetItem(pyList, 0));
    if (*numCols <= 0) {
        return NULL;
    }
    double** array = (double**)malloc(*numRows * sizeof(double*));
    CHECK_MEMORY_ALLOCATION(array);
    /*Iterate over the outer list*/
    for (int i = 0; i < numRows; i++) {
        PyObject *innerList = PyList_GetItem(pyList, i);
        array[i] = (double*)malloc(*numCols * sizeof(double));
        CHECK_MEMORY_ALLOCATION(array[i]);

        /*Iterate over the inner list*/
        for (int j = 0; j < numCols; j++) {
            PyObject *coordObj = PyList_GetItem(innerList, j);
            array[i][j] = PyFloat_AsDouble(coordObj);
        }
    }
    return array;
}

static PyObject* convertCMatrixToPyList(double **matrix, int numPoints)
{
    PyObject* py_matrix = PyList_New(numPoints);
    for (int i = 0; i < numPoints; i++) {
        PyObject* row = PyList_New(numPoints);
        for (int j = 0; j < numPoints; j++) {
            PyList_SetItem(row, j, PyFloat_FromDouble(matrix[i][j]));
        }
        PyList_SetItem(py_matrix, i, row);
    }
    return py_matrix;
}

static PyObject* py_sym(PyObject *self, PyObject *args)
{
    PyObject *dataPoints;
    int numPoints, dimensions;
    PyObject* py_matrix;
    point* points;
    double** A;
    if(!PyArg_ParseTuple(args, "O", &dataPoint)) {
        return NULL; /* In the CPython API, a NULL value is never valid for a
                        PyObject* so it is used to signal that an error has occurred. */
    }
    numPoints = PyList_size(dataPoints);
    dimensions = PyList_size(dataPoints[i]);
    points = convertPyListToPoints(dataPoints,numPoints, dimensions);
    A = sym(points, numPoints, dimensions);
    py_matrix = convertCMatrixToPyList(A, numPoints);
    freeArray(A, numPoints);
    freePoints(points, numPoints);
    return py_matrix;
}

static PyObject* py_ddg(PyObject *self, PyObject *args)
{
    PyObject *dataPoints;
    int numPoints, dimensions;
    double** A, D;
    PyObject* py_matrix;
    point* points;
    if(!PyArg_ParseTuple(args, "O", &dataPoint)) {
        return NULL; /* In the CPython API, a NULL value is never valid for a
                        PyObject* so it is used to signal that an error has occurred. */
    }
    numPoints = PyList_size(dataPoints);
    dimensions = PyList_size(dataPoints[i]);
    points = convertPyListToPoints(dataPoints, numPoints, dimensions);
    A = sym(points, numPoints, dimensions);
    D = ddg(A, numPoints);
    py_matrix = convertCMatrixToPyList(D, numPoints);
    freeArray(A, numPoints);
    freeArray(D, numPoints);
    freePoints(points, numPoints);
    return py_matrix;
}
static PyObject* py_norm(PyObject *self, PyObject *args)
{
    PyObject *dataPoints;
    int numPoints, dimensions;
    double** A, D, N;
    PyObject* py_matrix;
    point* points;
    if(!PyArg_ParseTuple(args, "O", &dataPoint)) {
        return NULL; /* In the CPython API, a NULL value is never valid for a
                        PyObject* so it is used to signal that an error has occurred. */
    }
    numPoints = PyList_size(dataPoints);
    dimensions = PyList_size(dataPoints[i]);
    points = convertPyListToPoints(dataPoints, numPoints, dimensions);
    A = sym(points, numPoints, dimensions);
    D = ddg(A, numPoints);
    N = norm(A, D, numPoints);
    py_matrix = convertCMatrixToPyList(N, numPoints);
    freeArray(A, numPoints);
    freeArray(D, numPoints);
    freeArray(N, numPoints);
    freePoints(points, numPoints);
    return py_matrix;
}

static PyObject* py_symnmf(PyObject *self, PyObject *args)
{
    PyObject *W, *H;
    double** arrayH, arrayW;
    int k, numRows, numColums;
    double eps;
    matrix* matH, matW newH;
    /* This parses the Python arguments into a double (d)  variable named z and int (i) variable named n*/
    if(!PyArg_ParseTuple(args, "iiOiO", &numRows, &numColums , &H, &W)) {
        return NULL; /* In the CPython API, a NULL value is never valid for a
                        PyObject* so it is used to signal that an error has occurred. */
    }
    arrayH = convertPyListToArray(*H,numRows,numColums);
    arrayW = convertPyListToArray(*W,numRows,numRows);
    
    matH = create_matrix(*arrayH, numRows, numColums);
    matW = create_matrix(*arrayW, numRows, numRows);
    convertCMatrixToPyList(arrayH);
    convertCMatrixToPyList(arrayW);
    newH = updateH(*matH, *matW, ITER, EPSILON);

    // convert newH to pyList 
    
    //PyObject* py_matrix = convertCMatrixToPyList(newH ,numPoints);
    // Free the C matrix and points
    freeArray(updatedH, numRows);
    return py_matrix;   
}

static PyMethodDef SymNMFMethods[] = {
    {"sym",(PyCFunction) py_sym, METH_VARARGS, pyDoc_STR("Compute the similarity matrix")},
    {"ddg",(PyCFunction) py_ddg, METH_VARARGS, pyDoc_STR("Compute the Diagonal Degree Matrix")},
    {"norm",(PyCFunction) py_norm, METH_VARARGS, pyDoc_STR("Compute the normalized similarity matrix")},
    {"symnmf", py_symnmf, METH_VARARGS, pyDoc_STR("Perform the full symNMF and return H")},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef symnmfmodule = {
    PyModuleDef_HEAD_INIT,
    "mysymnmf",
    NULL,
    -1,
    SymNMFMethods
};

PyMODINIT_FUNC PyInit_mysymnmf(void) {
    PyObject *m;
    m = PyModule_Create(&symnmfmodule);
    if (!m){
        return NULL;
    }
    return m;
}