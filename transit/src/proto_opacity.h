#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* src/opacity.c */
extern int opacity P_((struct transit *tr));
extern int readopacity P_((struct transit *tr, FILE *fp));
extern int calcopacity P_((struct transit *tr, FILE *fp));
extern int extinction P_((struct transit *tr, int r, int t));
extern int freemem_opacity P_((struct opacity *op, long *pi));

#undef P_
