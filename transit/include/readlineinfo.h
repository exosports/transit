// Copyright (C) 2015-2016 University of Central Florida. All rights reserved.
// Transit is under an open-source, reproducible-research license (see LICENSE).

#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* src/readlineinfo.c                                                       */
extern void datafileBS P_((FILE *fp, PREC_NREC offs, PREC_NREC nfields,
                           PREC_LNDATA target, PREC_NREC *resultp,
                           int reclength, int up));
extern int readtli_bin P_((FILE *fp, struct transit *tr, struct lineinfo *li));
extern int setimol P_((struct transit *tr));
extern int checkrange P_((struct transit *tr, struct lineinfo *li));
extern int readinfo_tli P_((struct transit *tr, struct lineinfo *li));
extern int readdatarng P_((struct transit *tr, struct lineinfo *li));
extern int readlineinfo P_((struct transit *tr));
extern int freemem_isotopes P_((struct isotopes *iso, long *pi));
extern int freemem_lineinfo P_((struct lineinfo *li, long *pi));
extern int freemem_linetransition P_((struct line_transition *lt, long *pi));
extern void saveline P_((FILE *fp, struct lineinfo *li));

#undef P_
