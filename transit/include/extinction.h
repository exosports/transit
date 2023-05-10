// Copyright (C) 2015-2016 University of Central Florida. All rights reserved.
// Transit is under an open-source, reproducible-research license (see LICENSE).

#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* src/extinction.c */
extern int getprofile P_((float **pr,         double dwn, float dop,
                                 float lor, float ta, int nwave));
extern void savefile_extinct P_((char *filename, double **e, short *c,
                                 long nrad, long nwav));
extern void restfile_extinct P_((char *filename, double **e, short *c,
                                 long nrad, long nwav));
extern void outputinfo P_((char *outfile, long w, long dw, long ln, long dln,
                           double **kiso, double timesalpha, int fbinvoigt,
                           double temp, double rad));
extern int extwn P_((struct transit *tr));
extern void printone P_((struct transit *tr));
extern int freemem_extinction P_((struct extinction *ex, long *pi));
extern int restextinct P_((FILE *in, PREC_NREC nrad, short niso, PREC_NREC nwn,
                           struct extinction *ex));
extern int computemolext P_((struct transit *tr, PREC_RES **kiso,
                   PREC_ATM temp, PREC_ATM *density, double *Z, int permol));
extern int interpolmolext P_((struct transit *tr, PREC_NREC r, PREC_RES **kiso));
extern void computeextscat P_((double *e, long n, 
                               struct extscat *sc,
                               double *pressure,
                               double *temp,
                               struct molecules *mol,
                               double wn));

extern void computeextcloud P_((double *e, long n,
                                struct extcloud *cl, double *pressure,
                                double *temp, double tfct,
                                double *density,
                                double wn));
#undef P_
