#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* src/makesample.c */
extern int  makesample P_((prop_samp *samp, prop_samp *hint,
                           prop_samp *ref, const long fl));
extern int  makesample1 P_((prop_samp *samp, prop_samp *ref, const long fl));
extern int  makewnsample P_((struct transit *tr));
extern int  makeipsample P_((struct transit *tr));
extern int  makeradsample P_((struct transit *tr));
extern int  maketempsample P_((struct transit *tr));
extern void savesample P_((FILE *out, prop_samp *samp));
extern void savesample_arr P_((FILE *out, prop_samp *samp));
extern int restsample P_((FILE *in, prop_samp *samp));
extern int restsample_arr P_((FILE *in, prop_samp *samp));
extern int outsample P_((struct transit *tr));
extern void freemem_samp P_((prop_samp *samp));
/* Spline functions:                                                        */
extern inline double * tri P_((double *a, double *d, double *c, double *b,
                               double *e, long n));
extern inline double * spline3 P_((double *xi, double *yi, double *x, 
                                   double *z,  double *h,  double *y,
                                   long nx, long N));
extern inline double * splinterp P_((long N, double *xi, double *yi, long nx, 
                                     double *xout, double *yout));
extern double splinterp_pt P_((double *z, long N, double *x, double *y,
                               double xout, double yout));
extern double * spline_init P_((double *z, double *x, double *y, long N));

#undef P_
