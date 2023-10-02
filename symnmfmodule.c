#define PY_SSIZE_T_CLEAN 
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include "symnmf.h"


point* convertPyListToPoints(PyObject *pyList, int n, int d) {
    int i,j;
    point* points = (point*)malloc(n * sizeof(point));
    if(checkMemoryAllocation(points) == 0){
        exit(1);
    }
    /*Iterate over the outer list*/
    for (i = 0; i < n; i++) {
        PyObject *innerList = PyList_GetItem(pyList, i);

        /*Allocate memory for the coordinates of the point*/
        points[i].coordinates = (double*)malloc(d * sizeof(double));
        if(checkMemoryAllocation(points[i].coordinates) == 0){
            for(j=0; j < i; j++){
                free(points[j].coordinates);
            }
            free(points);
        }

        /*Iterate over the inner list*/
        for (j = 0; j < d; j++) {
            PyObject *coordObj = PyList_GetItem(innerList, j);
            points[i].coordinates[j] = PyFloat_AsDouble(coordObj);
        }
    }
    return points;
}


double** convertPyListToArray(PyObject *pyList, int numRows, int numCols) {
    int i,j;
    numRows = PyList_Size(pyList);
    if (numRows <= 0) {
        return NULL;
    }

    numCols = PyList_Size(PyList_GetItem(pyList, 0));
    if (numCols <= 0) {
        return NULL;
    }
    double** array = (double**)malloc(numRows * sizeof(double*));
    if(checkMemoryAllocation(array)==0){
        return NULL;
    }
    /*Iterate over the outer list*/
    for (i = 0; i < numRows; i++) {
        PyObject *innerList = PyList_GetItem(pyList, i);
        array[i] = (double*)malloc(numCols * sizeof(double));
        if(checkMemoryAllocation(array[i]) == 0){
            for(j = 0; j < i; j++){
                free(array[j]);
            }
            free(array);
            return NULL;
        }

        /*Iterate over the inner list*/
        for (j = 0; j < numCols; j++) {
            PyObject *coordObj = PyList_GetItem(innerList, j);
            array[i][j] = PyFloat_AsDouble(coordObj);
        }
    }
    return array;
}

static PyObject* convertCMatrixToPyList(matrix * m)
{
    int i,j;
    int r = m -> r;
    int c = m->c;

    PyObject* py_matrix = PyList_New(r);
    for (i = 0; i < r; i++) {
        PyObject* row = PyList_New(c);
        for (j = 0; j < c; j++) {
            PyList_SetItem(row, j, PyFloat_FromDouble((m->data)[i][j]));
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
    matrix* A;
    if(!PyArg_ParseTuple(args, "O", &dataPoints)) {
        return NULL; /* In the CPython API, a NULL value is never valid for a
                        PyObject* so it is used to signal that an error has occurred. */
    }
    numPoints = PyList_Size(dataPoints);
    dimensions = PyList_Size(PyList_GetItem(dataPoints, 0));
    points = convertPyListToPoints(dataPoints,numPoints, dimensions);
    A = sym(points, numPoints, dimensions);
    py_matrix = convertCMatrixToPyList(A);
    freeMatrix(A);
    freePoints(points, numPoints);
    return py_matrix;
}

static PyObject* py_ddg(PyObject *self, PyObject *args)
{
    PyObject *dataPoints;
    int numPoints, dimensions;
    matrix * A;
    matrix * D;
    PyObject* py_matrix;
    point* points;
    if(!PyArg_ParseTuple(args, "O", &dataPoints)) {
        return NULL; /* In the CPython API, a NULL value is never valid for a
                        PyObject* so it is used to signal that an error has occurred. */
    }
    numPoints = PyList_Size(dataPoints);
    dimensions = PyList_Size(PyList_GetItem(dataPoints, 0));
    points = convertPyListToPoints(dataPoints, numPoints, dimensions);
    A = sym(points, numPoints, dimensions);
    D = ddg(A);
    py_matrix = convertCMatrixToPyList(D);
    freeMatrix(A);
    freeMatrix(D);
    freePoints(points, numPoints);
    return py_matrix;
}
static PyObject* py_norm(PyObject *self, PyObject *args)
{
    PyObject *dataPoints;
    int numPoints, dimensions;
    matrix* A;
    matrix* D;
    matrix* N;
    PyObject* py_matrix;
    point* points;
    if(!PyArg_ParseTuple(args, "O", &dataPoints)) {
        return NULL; /* In the CPython API, a NULL value is never valid for a
                        PyObject* so it is used to signal that an error has occurred. */
    }
    numPoints = PyList_Size(dataPoints);
    dimensions = PyList_Size(PyList_GetItem(dataPoints, 0));
    points = convertPyListToPoints(dataPoints, numPoints, dimensions);
    A = sym(points, numPoints, dimensions);
    D = ddg(A);
    N = norm(A, D);
    py_matrix = convertCMatrixToPyList(N);
    freeMatrix(A);
    freeMatrix(D);
    freeMatrix(N);
    freePoints(points, numPoints);
    return py_matrix;
}

static PyObject* py_symnmf(PyObject *self, PyObject *args)
{
    PyObject *W, *H;
    double** arrayH;
    double** arrayW;
    int numRows, numColums;
    matrix* matH;
    matrix* matW;
    matrix* newH;
    /* This parses the Python arguments into a double (d)  variable named z and int (i) variable named n*/
    if(!PyArg_ParseTuple(args, "iiOO", &numRows, &numColums , &H, &W)) {
        return NULL; /* In the CPython API, a NULL value is never valid for a
                        PyObject* so it is used to signal that an error has occurred. */
    }
    arrayH = convertPyListToArray(H,numRows,numColums);
    arrayW = convertPyListToArray(W,numRows,numRows);
    
    matH = create_matrix(arrayH, numRows, numColums);
    matW = create_matrix(arrayW, numRows, numRows);
    newH = updateH(matH, matW);

    // convert newH to pyList 
    
    PyObject* py_matrix = convertCMatrixToPyList(newH);
    // Free the C matrix and points
    freeMatrix(newH);
    freeMatrix(matH);
    freeMatrix(matW);
    return py_matrix;   
}

static PyMethodDef SymNMFMethods[] = {
    {"sym",(PyCFunction) py_sym, METH_VARARGS, "Compute the similarity matrix"},
    {"ddg",(PyCFunction) py_ddg, METH_VARARGS, "Compute the Diagonal Degree Matrix"},
    {"norm",(PyCFunction) py_norm, METH_VARARGS, "Compute the normalized similarity matrix"},
    {"symnmf", py_symnmf, METH_VARARGS, "Perform the full symNMF and return H"},
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
