/*
 * idxref.c
 * idxref.txc - Finds the index of refraction. Component of the Transit
 *              program.
 *
 * Copyright (C) 2003 Patricio Rojo (pato@astro.cornell.edu)
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

#include <transit.h>

/* \fcnfh
   Calculates the index of refraction. Right now it only gives 1 to all
   levels.

   @returns 0 on success
 */
int
idxrefrac(struct transit *tr)
{
  static struct idxref idx;
  long r;

  transitcheckcalled(tr->pi,"idxrefrac",2,
		     "getatm",TRPI_GETATM,
		     "makeradsample",TRPI_MAKERAD
		     );

  tr->ds.ir=&idx;

  idx.n=(PREC_RES *)calloc(tr->rads.n,sizeof(PREC_RES));
  for(r=0;r<tr->rads.n;r++)
    idx.n[r]=1;

  tr->pi|=TRPI_IDXREFRAC;
  return 0;
}
