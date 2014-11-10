/****************************** START LICENSE ******************************
Transit, a code to solve for the radiative-transifer equation for
planetary atmospheres.

This project was completed with the support of the NASA Planetary
Atmospheres Program, grant NNX12AI69G, held by Principal Investigator
Joseph Harrington. Principal developers included graduate students
Patricio E. Cubillos and Jasmina Blecic, programmer Madison Stemm, and
undergraduate Andrew S. D. Foster.  The included
'transit' radiative transfer code is based on an earlier program of
the same name written by Patricio Rojo (Univ. de Chile, Santiago) when
he was a graduate student at Cornell University under Joseph
Harrington.

Copyright (C) 2014 University of Central Florida.  All rights reserved.

This is a test version only, and may not be redistributed to any third
party.  Please refer such requests to us.  This program is distributed
in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.

Our intent is to release this software under an open-source,
reproducible-research license, once the code is mature and the first
research paper describing the code has been accepted for publication
in a peer-reviewed journal.  We are committed to development in the
open, and have posted this code on github.com so that others can test
it and give us feedback.  However, until its first publication and
first stable release, we do not permit others to redistribute the code
in either original or modified form, nor to publish work based in
whole or in part on the output of this code.  By downloading, running,
or modifying this code, you agree to these conditions.  We do
encourage sharing any modifications with us and discussing them
openly.

We welcome your feedback, but do not guarantee support.  Please send
feedback or inquiries to:

Joseph Harrington <jh@physics.ucf.edu>
Patricio Cubillos <pcubillos@fulbrightmail.org>
Jasmina Blecic <jasmina@physics.ucf.edu>

or alternatively,

Joseph Harrington, Patricio Cubillos, and Jasmina Blecic
UCF PSB 441
4111 Libra Drive
Orlando, FL 32816-2385
USA

Thank you for using transit!
******************************* END LICENSE ******************************/


/* List of functions defined:
int idxrefrac(struct transit *tr)
   Calculates the index of refraction.  Currently, it sets an index of
   refraction of 1.0 at all levels (no light bending).

int freemem_idexrefrac(struct idxref *ir, long *pi)
   Free index of refraction array

int restidxref(FILE *in, PREC_NREC nrad, struct idxref *ir)
   Restore hints structure, the structure needs to have been
   allocated before.

void saveidxref(FILE *out, PREC_NREC nrad, struct idxref *ir)
   Write index of refraction values to file pointed by out.
*/


#include <transit.h>

/* \fcnfh
   Calculates the index of refraction.  Currently, it sets an index of
   refraction of 1.0 at all levels (no light bending).

   Return: 0 on success                                              */
int
idxrefrac(struct transit *tr){
  static struct idxref st_idx;
  long r;            /* Radius index */
  PREC_ATM rho;      /* Density      */
  PREC_ATM nustp=0;  /* FINDME: Explain my name */
  /* TD: Allow for ray bending. Tau2 has to be enabled as well */

  /* Check radius array has been already sampled: */
  transitcheckcalled(tr->pi, "idxrefrac", 1, "makeradsample", TRPI_MAKERAD);

  /* Get struct objects: */
  tr->ds.ir     = &st_idx;
  prop_atm *atm = &tr->atm;

  /* Allocate space and initialize: */
  st_idx.n = (PREC_RES *)calloc(tr->rads.n, sizeof(PREC_RES));

  /* Calculate density at each radius: */
  for(r=0; r<tr->rads.n; r++){
    rho = stateeqnford(1, 1.0, atm->mm[r], 0, atm->p[r], atm->t[r]);
    st_idx.n[r] = 1 + rho*nustp/(LO*AMU*atm->mm[r]);
  }

  /* Set progress indicator and return success: */
  tr->pi |= TRPI_IDXREFRAC;
  return 0;
}


/* \fcnfh
   Free index of refraction array

   Return: 0 on success                             */
int
freemem_idexrefrac(struct idxref *ir, /* Index of refraction structure */
                   long *pi){
  /* Free arrays: */
  free(ir->n);

  /* Update progress indicator and return: */
  *pi &= ~(TRPI_IDXREFRAC|TRPI_TAU);
  return 0;
}


/* \fcnfh
   Restore hints structure, the structure needs to have been
   allocated before

   Return: 0 on success,
          -1 if not all the expected information is read
          -2 if info read is wrong
          -3 if cannot allocate memory
           1 if information read was suspicious                 */
int
restidxref(FILE *in,
           PREC_NREC nrad,
           struct idxref *ir){

  if(nrad<0)
    return -2;
  if(nrad>1000000)
    return 1;
  if((ir->n=(PREC_RES *)calloc(nrad, sizeof(PREC_RES)))==NULL)
    return -3;
  if(nrad==0)
    return 0;
  if(fread(ir->n, sizeof(PREC_RES), nrad, in) != nrad)
    return -1;

  return 0;
}


/* \fcnfh
   Write index of refraction values to file 'out'  */
void
saveidxref(FILE *out,
           PREC_NREC nrad,
           struct idxref *ir){
  if(nrad>0)
    fwrite(ir->n, sizeof(PREC_RES), nrad, out);
}
