%module transit_module
%{
#define SWIG_FILE_WITH_INIT
//extern struct transit transit;
//extern long itr;
//extern struct timeval tv;
//extern double t0;
//extern int argc;
//extern char ** argv;
//extern int init_run;

extern void transit_init(int argc, char **argv);
extern int  get_no_samples(void);
extern void get_waveno_arr(double * waveno_arr, int waveno);
extern void run_transit(double * re_input, int transint, double *\
transit_out,int transit_out_size);
extern void free_memory(void);
%}

%include "numpy.i"

%init %{
        import_array();
%}

%apply (double* INPLACE_ARRAY1,int DIM1) {(double* waveno_arr, int waveno)}
%apply (double* INPLACE_ARRAY1,int DIM1) {(double* transit_out, int transit_out_size)}
%apply (double* IN_ARRAY1, int DIM1) {(double* re_input, int transint)}
/*%exception
{
     errno = 0;
     $action
 
     if (errno != 0)
     {
         switch(errno)
         {
             case ENOMEM:
                 PyErr_Format(PyExc_MemoryError, "Failed malloc()");
                 break;
             default:
                 PyErr_Format(PyExc_Exception, "Unknown exception");
         }
         SWIG_fail;
     }
}
*/

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
//extern struct transit transit;
//extern long itr;
//extern struct timeval tv;
//extern double t0;
//extern int argc;
//extern char ** argv;
//extern int init_run;

extern void transit_init(int argc, char **argv);
extern int  get_no_samples(void);
extern void get_waveno_arr(double * waveno_arr, int waveno);
extern void run_transit(double * re_input, int transint, double *\
transit_out,int transit_out_size);
extern void free_memory(void);

