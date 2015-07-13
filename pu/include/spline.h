#ifndef _SPLINE_H
#define _SPLINE_H

#include <stdio.h>
#include <stdlib.h>

#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* spline.c                                                                 */
extern inline double * tri P_((double *a, double *d, double *c, double *b,
                               double *e, long n));
extern inline double * spline3 P_((double *xi, double *yi, double *x,
                                   double *z,  double *h,  double *y,
                                   long nx, long N));
extern double * splinterp P_((long N, double *xi, double *yi, long nx,
                                     double *xout, double *yout));
extern double splinterp_pt P_((double *z, long N, double *x, double *y,
                               double xout, double yout));
extern double * spline_init P_((double *z, double *x, double *y, long N));

#undef P_

#endif  /* _SPLINE_H                                                        */
