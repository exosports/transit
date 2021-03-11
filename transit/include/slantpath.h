// Copyright (C) 2015-2016 University of Central Florida. All rights reserved.
// Transit is under an open-source, reproducible-research license (see LICENSE).

#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* src/slantpath.c */
extern int modulation P_((struct transit *tr));
extern void printmod P_((struct transit *tr));
extern int freemem_outputray P_((struct outputray *out, long *pi));

#undef P_

