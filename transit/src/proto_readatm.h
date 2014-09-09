#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* src/atmosphere/readatm.c */
extern int freemem_atmosphere P_((struct atm_data *at, long *pi));
extern double checkaddmm P_((double *mm, long r, prop_mol *molec, struct molecules *mol, int n, short mass));
extern void telldefaults P_((struct isotopes *iso, struct atm_data *at));
extern int getatm P_((struct transit *tr));
extern int getmnfromfile P_((FILE *fp, struct atm_data *at, struct transit *tr, PREC_ZREC *f_remainder));
extern int readatmfile P_((FILE *fp, struct transit *tr, struct atm_data *at, prop_samp *rads, int nrad, PREC_ZREC *f_remainder));
extern void storename P_((struct atm_data *at, char *line));
extern void getmoldata P_((struct atm_data *at, struct molecules *mol));
extern int reloadatm P_((struct transit *tr, double *input));
#undef P_




