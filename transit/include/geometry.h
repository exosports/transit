// Copyright (C) 2015-2016 University of Central Florida. All rights reserved.
// Transit is under an open-source, reproducible-research license (see LICENSE).

#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* src/geometry.c */
extern int setgeomhint P_((struct transit *tr));
extern int setgeom P_((struct geometry *sg, double time, long *flags));
extern inline double starvariation P_((double x, double y, double radius));

#undef P_
