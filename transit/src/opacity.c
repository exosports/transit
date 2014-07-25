/*
 * opacity.c - Calculate the opacity for a grid of temperature and pressures
 *             (or radius?) and store to file.
 *
 * Copyright (C) 2014 Patricio Cubillos
 *
 */


/* List of functions defined:
int opacity(struct transit *tr)
   Calculates the index of refraction.  Currently, it sets an index of
   refraction of 1.0 at all levels (no light bending).

int freemem_opacity(struct opacity *op, long *pi)
   Free index of refraction array

int restidxref(FILE *in, PREC_NREC nrad, struct idxref *ir)
   Restore hints structure, the structure needs to have been
   allocated before.

void saveidxref(FILE *out, PREC_NREC nrad, struct idxref *ir)
   Write index of refraction values to file pointed by out.
*/


#include <transit.h>

/* \fcnfh
   Calculate the opacity due to molecular transitions.

   Return: 0 on success                                                     */
int
opacity(struct transit *tr){
  //static struct opacity op;  /* The opacity struct  */
  long t, m, r, w;           /* Radius indices      */

  t = m = r = w = 1;
  w = t;
  /* Check radius array has been already sampled: */
  transitcheckcalled(tr->pi, "idxrefrac", 1, "makeradsample", TRPI_MAKERAD);

  transitprint(1, verblevel, "This is opacity grid talking to you!\n");
  /* Check if the opacity file exists:                                      */

  /* If the file exists, read the opacities from it:                        */

  /* If the file does not exists, calculate the opacities from hints:       */

  /* Get struct objects: */
  //tr->ds.ir     = &st_idx;
  //prop_atm *atm = &tr->atm;

  /* Allocate space and initialize: */
  //st_idx.n = (PREC_RES *)calloc(tr->rads.n, sizeof(PREC_RES));

  /* Calculate density at each radius: */
  //for(r=0; r<tr->rads.n; r++){
  //  rho = stateeqnford(1, 1.0, atm->mm[r], 0, atm->p[r], atm->t[r]);
  //  st_idx.n[r] = 1 + rho*nustp/(LO*AMU*atm->mm[r]);
  //}

  /* Set progress indicator and return success: */
  tr->pi |= TRPI_OPACITY;
  return 0;
}


/* \fcnfh
   Free opacity array

   Return: 0 on success                                                     */
int
freemem_opacity(struct opacity *op, /* Opacity structure                    */
                long *pi){          /* transit progress flag                */
  /* Free arrays: */
  free(op->opacity[0][0][0]);
  free(op->opacity[0][0]);
  free(op->opacity[0]);
  free(op->opacity);

  /* Update progress indicator and return:                                  */
  *pi &= ~(TRPI_OPACITY | TRPI_TAU);
  return 0;
}

