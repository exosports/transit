// Copyright (C) 2015-2016 University of Central Florida. All rights reserved.
// Transit is under an open-source, reproducible-research license (see LICENSE).

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
extern inline void tri P_((double *h, double *y, double *z, long n));
extern inline void spline3 P_((double *xi, double *yi, double *x,
                                   double *z,  double *h,  double *y,
                                   long nx, long N));
extern void splinterp P_((long N, double *xi, double *yi, long nx,
                                     double *xout, double *yout));
extern double splinterp_pt P_((double *z, long N, double *x, double *y,
                               double xout));
extern void spline_init P_((double *z, double *x, double *y, long N));

#undef P_

#endif  /* _SPLINE_H                                                        */
