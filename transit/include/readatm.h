/*
 * readatm.h - Read atmospheric info common headers. Component of the
 *             Transit program.
 *
 * Copyright (C) 2004-2006 Patricio Rojo (pato@astro.cornell.edu)
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of version 2 of the GNU General 
 * Public License as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 * 02111-1307, USA.
 */

#ifndef _READATM_H
#define _READATM_H

#include <transit.h>

extern char *atmfilename;

#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* src/atmosphere/readatm.c */
extern int freemem_atmosphere P_((struct atm_data *at, long *pi));
extern double checkaddmm P_((double *mm, long r, prop_mol *molec, struct molecules *mol, int n, short mass));
extern void telldefaults P_((struct isotopes *iso, struct atm_data *at));
extern int getatm P_((struct transit *tr));
extern int getmnfromfile P_((FILE *fp, struct atm_data *at, struct transit *tr));
extern int readatmfile P_((FILE *fp, struct transit *tr, struct atm_data *at, prop_samp *rads, int nrad));
extern void storename P_((struct atm_data *at, char *line));
extern void getmoldata P_((struct atm_data *at, struct molecules *mol, char * filename));
extern int reloadatm P_((struct transit *tr, double *input));
extern int radpress P_((double g, double p0, double r0, double *temp,
                        double *mu, double *pressure, double *radius,
                        int nlayer, double rfct));
#undef P_

#endif
