// Copyright (C) 2015-2016 University of Central Florida. All rights reserved.
// Transit is under an open-source, reproducible-research license (see LICENSE).

#ifndef _READATM_H
#define _READATM_H

#include <transit.h>

extern char *atmfilename;

#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* src/atmosphere/readatm.c */
extern int freemem_atmosphere P_((struct atm_data *at, long *pi));
extern double checkaddmm P_((double *mm, PREC_NREC r, prop_mol *molec, struct molecules *mol, int n, short mass));
extern void telldefaults P_((struct isotopes *iso, struct atm_data *at));
extern int getatm P_((struct transit *tr));
extern int getmnfromfile P_((FILE *fp, struct atm_data *at, struct transit *tr));
extern int readatmfile P_((FILE *fp, struct transit *tr, struct atm_data *at, prop_samp *rads, int nrad));
extern void storename P_((struct atm_data *at, char *line));
extern void getmoldata P_((struct atm_data *at, struct molecules *mol, char * filename));
extern int reloadatm P_((struct transit *tr, double *input));
extern int radpress P_((double g, double p0, double r0, double *temp,
                        double *mu, double *pressure, double *radius,
                        int nlayer, double rfct));
#undef P_

#endif
