// Copyright (C) 2015-2016 University of Central Florida. All rights reserved.
// Transit is under an open-source, reproducible-research license (see LICENSE).

#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* src/opacity.c */
extern int opacity P_((struct transit *tr));
extern int calcprofiles P_((struct transit *tr));
extern int calcopacity P_((struct transit *tr, FILE *fp));
extern int readopacity P_((struct transit *tr, FILE *fp));
extern int shareopacity P_((struct transit *tr, FILE *fp));
extern int attachopacity P_((struct transit *tr));
extern int mountopacity P_((struct transit *tr));
extern int freemem_opacity P_((struct opacity *op, long *pi));

#undef P_
