// Copyright (C) 2015-2016 University of Central Florida. All rights reserved.
// Transit is under an open-source, reproducible-research license (see LICENSE).

#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* src/crosssec.c */
extern int readcs P_((struct transit *tr));
extern int interpcs P_((struct transit *tr));
extern int bicubicinterpolate P_((double **res, double **src, double *x1, long nx1, double *x2, long nx2, double *t1, long nt1, double *t2, long nt2));
extern void cserr P_((int max, char *name, int line));
extern int freemem_cs P_((struct cross *cross, long *pi));

#undef P_
