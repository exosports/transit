/*
 * numerical.h - Miscellaneous numerical utilities header file
 *
 * Copyright (C) 2004 Patricio Rojo (pato@astro.cornell.edu)
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
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

#ifndef _NUMERICAL_H
#define _NUMERICAL_H

#include <stdio.h>
#include <stdlib.h>

#if __STDC__ || defined(__cplusplus)
#define P_(s) s
#else
#define P_(s) ()
#endif

/* numerical.c */
extern inline int binsearchie P_((double *arr, long i, long f, double val));
extern inline int binsearchei P_((double *arr, long i, long f, double val));
extern int binsearch P_((double *arr, long i, long f, double val));
extern double integ_trasim P_((double dx, double *y, long n));
extern double integ_trapz  P_((double *x, double *y, long n));
extern double interp_parab P_((double *x, double *y, double xr));
extern double interp_line P_((double *x, double *y, double xr));
extern int  downsample P_((double *input, double *out, int n, int scale));
extern double powi P_((double x, int n));
extern _Bool fixedcmp P_((double d1, double d2, int prec));
/* Simpson integration functions:                                          */
extern void geth P_((double *h, double *hsum, double *hratio, double *hfactor,
                     int n));
extern double simps P_((double *y, double *h, double *hsum, double *hratio,
                        double *hfactor, int n));
extern void makeh P_((double *x, double *h, int n));
extern inline double simpson P_((double *y, double *hsum, double *hratio,
                                 double *hfactor, int n));

#undef P_

#endif /* _NUMERICAL_H */
