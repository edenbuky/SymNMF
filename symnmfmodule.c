#define PY_SSIZE_T_CLEAN 
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include "symnmf.h"



static PyObject* fit(PyObject *self, PyObject *args)
{
    PyObject *dataPoints, *W, *H;
    int k;
    const char *goal;
    double** matrix;
    int numPoints, dimensions;
    /* This parses the Python arguments into a double (d)  variable named z and int (i) variable named n*/
    if(!PyArg_ParseTuple(args, "sO|iOO", &goal, &dataPoint,&k, &H, &W)) {
        return NULL; /* In the CPython API, a NULL value is never valid for a
                        PyObject* so it is used to signal that an error has occurred. */
    }
    if(strcmp(goal, "sym") == 0){
        matrix = sym(points, numPoints, dimensions);
    }else if (strcmp(goal, "ddg") == 0) {
        double** A = sym(points, numPoints, dimensions);
        matrix = ddg(A, numPoints);
        for (int i = 0; i < numPoints; i++) {
            free(A[i]);
        }
        free(A);
    }else if (strcmp(goal, "norm") == 0) {
        double** A = sym(points, numPoints, dimensions);
        double** D = ddg(A, numPoints);
        matrix = norm(A, D, numPoints);
        // Free A and D after using
        for (int i = 0; i < numPoints; i++) {
            free(A[i]);
            free(D[i]);
        }
        free(A);
        free(D);
    }
    PyObject* py_matrix = PyList_New(numPoints);
    for (int i = 0; i < numPoints; i++) {
        PyObject* row = PyList_New(numPoints);
        for (int j = 0; j < numPoints; j++) {
            PyList_SetItem(row, j, PyFloat_FromDouble(matrix[i][j]));
        }
        PyList_SetItem(py_matrix, i, row);
    }

    // Free the C matrix and points
    for (int i = 0; i < numPoints; i++) {
        free(matrix[i]);
    }
    free(matrix);
    freePoints(points, numPoints);

    return py_matrix;   
}

static PyMethodDef SymNMFMethods[] = {
    {"fit",
     (PyCFunction) fit, 
     METH_VARARGS, pyDoc_STR("Execute the SymNMF function based on the goal")},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef symnmfmodule = {
    PyModuleDef_HEAD_INIT,
    "symnmf",
    NULL,
    -1,
    SymNMFMethods
};

PyMODINIT_FUNC PyInit_symnmf(void) {
    PyObject *m;
    m = PyModule_Create(&symnmfmodule);
    if (!m){
        return NULL;
    }
    return m;
}
