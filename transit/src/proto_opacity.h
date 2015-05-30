#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* src/opacity.c */
extern int opacity P_((struct transit *tr));
extern int readopacity P_((struct transit *tr, FILE *fp));
extern int calcprofiles P_((struct transit *tr));
extern int calcopacity P_((struct transit *tr, FILE *fp));
extern int freemem_opacity P_((struct opacity *op, long *pi));

#undef P_
