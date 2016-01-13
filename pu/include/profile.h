// Copyright (C) 2015-2016 University of Central Florida. All rights reserved.
// Transit is under an open-source, reproducible-research license (see LICENSE).

#ifndef _PROFILE_H
#define _PROFILE_H

#define PREC_VOIGT float
#define PREC_VOIGTP double

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define VOIGT_QUICK 0x00001   //Quick integration.
extern int _voigt_maxelements;

//function definition (voigt.h)
#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* voigt.c */
extern inline int voigtf P_((int nwn, float *wn, float wn0, double alphaL, double alphaD, float *vpro, double eps));
extern int voigtn P_((int nwn, double dwn, double alphaL, double alphaD, float **vpro, double eps, int flags));
extern inline int voigtn2 P_((int m, int nwn, double dwn, double alphaL, double alphaD, float **vpro, double eps, int flags));


#undef P_

#endif /* _PROFILE_H */

