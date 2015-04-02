/* Revision        March 19th,   2014 Jasmina Blecic
                   implemented switch eclipse/transit, added new functions  */
/* Revision        April 26th,   2014 Jasmina Blecic
                   implemented intensity grid and flux, added new functions */

#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* src/eclipse.c */
extern void printintens P_((struct transit *tr));
extern int emergent_intens P_((struct transit *tr));
extern int flux P_((struct transit *tr));
extern void printflux P_((struct transit *tr));
extern int freemem_intensityGrid P_((struct grid *intens, long *pi));
extern void freemem_localeclipse P_((void));

#undef P_
