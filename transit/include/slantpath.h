#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* src/slantpath.c */
extern static PREC_RES totaltau1 P_((PREC_RES b, PREC_RES *rad, PREC_RES refr,
                                     PREC_RES *ex, long nrad));
extern static PREC_RES totaltau2 P_((PREC_RES b, PREC_RES *rad, PREC_RES refr,
                                     PREC_RES *ex, long nrad));
extern static PREC_RES transittau P_((struct transit *tr, PREC_RES b, 
                                      PREC_RES *ex));
extern int modulation P_((struct transit *tr));
extern static PREC_RES modulation1 P_((PREC_RES *tau, long last, double toomuch,
                                       prop_samp *ip, struct geometry *sg));
inline static PREC_RES modulationm1 P_((PREC_RES *tau, long last,
                                        double toomuch, prop_samp *ip,
                                        struct geometry *sg));
extern static PREC_RES modulationperwn P_((struct transit *tr, PREC_RES *tau,
                                           PREC_RES w, long last, 
                                           double toomuch, prop_samp *ip));
extern void printmod P_((struct transit *tr));
extern int freemem_outputray P_((struct outputray *out, long *pi));

#undef P_

