// Copyright (C) 2015-2016 University of Central Florida. All rights reserved.
// Transit is under an open-source, reproducible-research license (see LICENSE).

#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* src/tau.c */
extern int init_optdepth P_((struct transit *tr));
extern int tau P_((struct transit *tr));
extern void outdebtauex P_((char *name, double **e, prop_samp *ip, double **t, long rn, long w));
extern void outdebex P_((char *name, double **e, double *r, long rn, long wi, long wf));
extern void outdebtau P_((char *name, prop_samp *ip, double **t, long wi, long wf));
extern void printtoomuch P_((char *file, struct optdepth *tau, prop_samp *wn, prop_samp *rad));
extern int freemem_tau P_((struct optdepth *tau, long *pi));
extern int detailout P_((prop_samp *wn, prop_samp *rad, struct detailfld *det, double **arr, short flag));
extern void print1dArrayDouble P_((FILE *outf, double *array, int noColumns, char *format));
extern void print2dArrayDouble P_((FILE *outf, double **array, int noRows, int noColumns, char *format, prop_samp *wn));
extern void savemolExtion P_((struct transit *tr, long ri));
extern void savetau P_((struct transit *tr));
extern void saveCIA P_((struct transit *tr));
void save1Darray P_((struct transit *tr, FILE *myFile, PREC_RES *array1d, int nrad, long wi));
FILE * openFile P_((char *filename, char *header));
void closeFile P_((FILE *myFile));

#undef P_
