// Copyright (C) 2015-2016 University of Central Florida. All rights reserved.
// Transit is under an open-source, reproducible-research license (see LICENSE).

#ifndef _PROCOPT_H
#define _PROCOPT_H

#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <iomisc.h>
#include <xmalloc.h>

#define HELPTITLE    0xff
#define ADDPARAMFILE 0xfe

extern unsigned short procopt_debug;

struct optdocs {
  char *name;
  int val;
  int has_arg;
  char *def;
  char *prm;
  char *doc;
};


struct optcfg {
  const char *prg;
  const char *options;
  const char *nonopt;
  const char *contact;
  const char *intro;
  const char *intro2;
  const char *contintro;
  const char *noinfo;
  const char *enddocchar;
  const char *defword;
  const char *postdefword;
  const char *pretitle;
  const char *posttitle;
  const char *postoption;
  _Bool freecfg;
  _Bool usestderr;
  char argmode;
  char endpadchar;
  char *files;			/* Configuration file */
  int indentdoc;
  int columns;
  short helpmode;		/* currently 0 or 1. help display mode */
};


//Functions declaration
#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* procopt.c */
extern int procopt P_((int argc, char **argv, struct optdocs *opts, struct optcfg *cfg, int *longidxp));
extern int getopt_long_files P_((int argc, char **argv, char *shortopts, struct option *getopts, int *longidxp, char *paramfilelist));
extern void prochelp P_((int status));
extern void procopt_free P_((void));

#undef P_


#endif /* _PROCOPT_H */
