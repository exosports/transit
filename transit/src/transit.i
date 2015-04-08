%module transit_module
%{
#define SWIG_FILE_WITH_INIT
extern int mytest(int argc, char **argv);
%}
%typemap(in) char ** {
  /* Check if is a list */
  if (PyList_Check($input)) {
    int size = PyList_Size($input);
    int i = 0;
    $1 = (char **) malloc((size+1)*sizeof(char *));
    for (i = 0; i < size; i++) {
      //PyObject *o = PyList_GetItem($input,i);
       $1[i] = PyString_AsString(PyList_GetItem($input,i));
    }
    $1[i] = 0;
  } else {
    PyErr_SetString(PyExc_TypeError,"not a list");
    return NULL;
  }
}

// This cleans up the char ** array we mallocd before the function call
%typemap(freearg) char ** {
  free((char *) $1);
}
#include "transit.h"
extern int mytest(int argc, char **argv);
