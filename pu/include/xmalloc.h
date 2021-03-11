// Copyright (C) 2015-2016 University of Central Florida. All rights reserved.
// Transit is under an open-source, reproducible-research license (see LICENSE).

#ifndef _XMALLOC_H
#define _XMALLOC_H

#define malloc(n) xmalloc(n)
#define realloc(p,n) xrealloc(p,n)
#define calloc(n,s) xcalloc(n,s)
//#define strdup(p) xstrdup(p)

//functin definition (xmalloc.h)
#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* xmalloc.c */
extern void *xmalloc P_((size_t n));
extern void *xcalloc P_((size_t n, size_t s));
extern void *xrealloc P_((void *p, size_t n));
extern char *xstrdup P_((char *str));

#undef P_


#endif /* _XMALLOC_H */
