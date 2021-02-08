// Copyright (C) 2015-2016 University of Central Florida. All rights reserved.
// Transit is under an open-source, reproducible-research license (see LICENSE).

/* Required for compilation with shared memory libraries.                   */
#define _XOPEN_SOURCE 700

#ifndef _TRANSIT_H
#define _TRANSIT_H

#include <stdarg.h>
#include <math.h>
#include <errno.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/types.h>
#include <unistd.h>
#include <sampling.h>
#include <profile.h>
#include <iomisc.h>
#include <numerical.h>
#include <spline.h>
#include <xmalloc.h>
#include <strings.h>
#include <stdlib.h>
#include <stdio.h>
#include <alloca.h>

#define compattliversion 6

#include <flags_tr.h>
#include <constants_tr.h>
#include <types_tr.h>

extern struct transit transit;
extern long itr;
extern struct timeval tv;
extern double t0;
extern int argc;
extern char ** argv;
extern int init_run;

extern void transit_init(int argc, char **argv);
extern int  get_no_samples(void);
extern void get_waveno_arr(double * waveno_arr, int waveno);
extern void set_radius(double refradius);
extern void run_transit(double * re_input, int transint, double *\
		transit_out,int transit_out_size);


/*****   Macros   *****/
static __inline__ void
printextprogress(long wi, long wnn){
}

/* FINDME: Take this function out of here                                   */
static __inline__ double
stateeqnford(_Bool mass,  /* Abundance by mass (1) or by number (0)         */
             double q,    /* Abundance                                      */
             double ma,   /* Mean molecular mass (in AMU)                   */
             double mi,   /* Molecular mass of the species (in AMU)         */
             double p,    /* Pressure                                       */
             double t){   /* Temperature                                    */
  const double rho = AMU * q * p / KB / t;
  if(mass)
    return rho * ma;
  return rho * mi;
}

#define tr_output(level, ...) do { \
  if (((level) & TOUT_VERBMASK) <= verblevel) tr_output_fcn(level, __FILE__, \
    __LINE__, __VA_ARGS__); \
} while(0)

#define DEBUG_ERROR
#ifdef  DEBUG_ERROR
#define DBGERR | TERR_DBG
#else
#define DBGERR
#endif /* DEBUG_ERROR */


#define transiterror(flag, ...) \
    transiterror_fcn(flag DBGERR, __FILE__, __LINE__, __VA_ARGS__)
#define vtransiterror(flag, ...) \
    vtransiterror_fcn(flag DBGERR, __FILE__, __LINE__, __VA_ARGS__)


long fw_status;
#define fw(fcn, failurecondition, ...) do{               \
    if((fw_status=fcn(__VA_ARGS__)) failurecondition) {  \
      tr_output(TOUT_ERROR,                              \
                   #fcn "() returned error code %li\n",  \
                   fw_status);                           \
      exit(EXIT_FAILURE);                                \
    }                                                    \
                                       }while(0)

/* Add flag to transit if present in hint and change hint flag value:       */
#define transitacceptflag(transit, hint, flag) do{                          \
        transit |= hint&flag; hint &= ~(flag); }while(0)

#define transitallocerror(nmb) do { \
  tr_output(TOUT_ERROR | TOUT_LOCATE | TOUT_BANNER, \
    "Allocation failed for %i allocation units. Impossible to continue.\n", \
    nmb); \
  exit(EXIT_FAILURE); \
} while(0)

#define free_null(x) do{free(x);x=NULL;}while(0)

#ifdef NODEBUG_TRANSIT
#define transitASSERT(...) ((void)0)
#else
#define free(x)  do{free(x);x=NULL;}while(0)
#define transitASSERT(a,...) do { \
  if(a) { \
    tr_output(TOUT_ERROR | TOUT_LOCATE, __VA_ARGS__); \
    exit(EXIT_FAILURE); \
  } \
} while(0)
#endif

/* Maximum name length: */
#define maxeisoname 20

#ifndef HAVE__BOOL
#define _Bool short
#endif

extern int verblevel;      /* Verbose level > 4 is only for debuging       */
extern int maxline;        /* Initialized in transitstd                     */
extern int version;
extern int revision;

#include <structures_tr.h>

extern const ray_solution slantpath, eclipsepath;

/* Transit Prototypes                                                       */
#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* src/transit.c */
/* FINDME: Main commented out to make MPI to work: */
//extern int main P_((int argc, char **argv));
extern void freemem_transit P_((struct transit *tr));

#undef P_

/***** Prototypes *****/
#include <readlineinfo.h>
#include <readatm.h>
#include <transitstd.h>
#include <makesample.h>
#include <extinction.h>
#include <opacity.h>
#include <idxrefraction.h>
#include <tau.h>
#include <argum.h>
#include <geometry.h>
#include <crosssec.h>
#include <eclipse.h>
#include <slantpath.h>
#endif /* _TRANSIT_H */
