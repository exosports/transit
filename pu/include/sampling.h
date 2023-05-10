// Copyright (C) 2015-2016 University of Central Florida. All rights reserved.
// Transit is under an open-source, reproducible-research license (see LICENSE).

#ifndef _SAMPLING_H
#define _SAMPLING_H

#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>

#define SAMP_BITS     0x0000000f
#define SAMP_LINEAR   0x00000001
#define SAMP_SPLINE   0x00000002

#define INTERP_LINEAR 0x0000001

//function definition (sampling.h)
#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* sampling.c */
extern int lineinterpol P_((int ndat, double *x, double *y, int n, long *indx,
                            float *t, double *yout, double *dbgout));
extern inline void natcubsplinecoef P_((long n, double *x, double *y,
                                        double *h, double *D));
extern inline double interp P_((double refx, double *x, double *y, long n,
                                int intkind));
extern double lineinterp P_((double refx, double *x, double *y, long n));

#undef P_

#endif /* _SAMPLING_H */
