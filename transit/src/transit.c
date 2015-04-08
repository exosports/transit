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

/* Revision        January 23rd, 2014 Jasmina Blecic
                   implemented eclipse                                      */
/* Revision        March 19th,   2014 Jasmina Blecic
                   implemented switch eclipse/transit                       */
/* Revision        April 26th,   2014 Jasmina Blecic
                   implemented intensity grid and flux                      */

/* TBD: calloc checks */

#include <transit.h>

/* \fcnfh                                                                   */
int mytest(int argc,      /* Number of variables                              */
         char **argv){  /* Variables                                        */

  /* Initialization of data structure's pointers. Note that s_lt is not
     assigned because the array that is going to point has not been
     initialized yet.                                                       */
  struct transit transit;
  long itr=0;
  struct timeval tv;
  double t0=0.0;

  memset(&transit, 0, sizeof(struct transit));
  verblevel=2;

  /* Process the command line arguments:                                    */
  fw(processparameters, !=0, argc, argv, &transit);
  t0 = timecheck(verblevel, itr,  0, "processparameters", tv, t0);

  /* Accept all general hints:                                              */
  fw(acceptgenhints, !=0, &transit);
  /* Presentation:                                                          */
  printintro();

  /* No program warnings if verblevel is 0 or 1:                            */
  if(verblevel<2)
    transit_nowarn = 1;

  /* Make wavenumber binning:                                               */
  fw(makewnsample, <0, &transit);
  t0 = timecheck(verblevel, itr,  1, "makewnsample", tv, t0);
  if(fw_status>0)
    transitprint(7, verblevel,
                 "makewnsample() modified some of the hinted "
                 "parameters according to returned flag: 0x%lx.\n",
                 fw_status);
  transitprint(10, verblevel, "Wavenumber [%li]: dwn=%.6f\n",
                              transit.wns.n, transit.wns.d);
  transitprint(10, verblevel, "Oversampled wavenumber [%li]: dwn=%.6f\n",
                              transit.owns.n, (transit.owns.d/transit.owns.o));

  /* Read Atmosphere information:                                           */
  fw(getatm, !=0, &transit);
  t0 = timecheck(verblevel, itr,  2, "getatm", tv, t0);

  /* Read line info:                                                        */
  fw(readlineinfo, !=0, &transit);
  t0 = timecheck(verblevel, itr,  3, "readlineinfo", tv, t0);

  /* Make radius binning and interpolate data to new value:                 */
  fw(makeradsample, <0, &transit);
  t0 = timecheck(verblevel, itr,  4, "makeradsample", tv, t0);
  if(fw_status>0)
    transitprint(7, verblevel, "makeradsample() modified some of the hinted "
                               "parameters. Flag: 0x%lx.\n", fw_status);

  // // Compare atmopheric samplings:
  // // radius
  // transitprint(1, 2, "Radius sampling size: %li\n", transit.rads.n);
  // for (itr = 0; itr < transit.rads.n; itr++)
  //   transitprint(1, 2, "%.3e ", transit.rads.v[itr] - transit.ds.at->rads.v[itr]);
  // // pressure
  // transitprint(1,2,"\nPressure sampling:\n");
  // for (itr = 0; itr < transit.rads.n; itr++)
  //   transitprint(1, 2, "%.3e ", transit.atm.p[itr] - transit.ds.at->atm.p[itr]);
  // // temp
  // transitprint(1,2,"\nTemperature sampling:\n");
  // for (itr = 0; itr < transit.rads.n; itr++)
  //   transitprint(1, 2, "%.3e ", transit.atm.t[itr] - transit.ds.at->atm.t[itr]);
  // // mm
  // transitprint(1,2,"\nMmm sampling:\n");
  // for (itr = 0; itr < transit.rads.n; itr++)
  //   transitprint(1, 2, "%.3e ", transit.ds.at->mm[itr]);
    //transitprint(1, 2, "%.3e ", transit.atm.mm[itr] - transit.ds.at->mm[itr]);
  // // q
  // transitprint(1,2,"\nAbundance sampling:\n");
  // for (itr = 0; itr < transit.rads.n; itr++)
  //   transitprint(1, 2, "%.3e ", transit.ds.mol->molec[0].q[itr] - transit.ds.at->molec[0].q[itr]);
  // // d
  // transitprint(1,2,"\nDensity sampling:\n");
  // for (itr = 0; itr < transit.rads.n; itr++)
  //   transitprint(1, 2, "%.3e ", transit.ds.mol->molec[1].d[itr] - transit.ds.at->molec[1].d[itr]);
  // // Z
  // transitprint(1,2,"\nPartition-function sampling:\n");
  // for (itr = 0; itr < transit.rads.n; itr++)
  //   transitprint(1, 2, "%.3e ", transit.ds.iso->isov[0].z[itr] -
  //                               transit.ds.li->isov[0].z[itr]);
  // transitprint(1,2,"\n");

  // return 0;
  /* Calculate opacity grid:                                                */
  fw(opacity, <0, &transit);
  t0 = timecheck(verblevel, itr,  5, "opacity", tv, t0);

  if (transit.opabreak){
    /* FINDME: Free memory                                                  */
    return EXIT_SUCCESS;
  }

  /* Compute sampling of impact parameter:                                  */
  fw(makeipsample, <0, &transit);
  t0 = timecheck(verblevel, itr,  6, "makeipsample", tv, t0);
  if (fw_status>0)
    transitprint(7, verblevel, "makeipsample() modified some of the hinted "
                               "parameters. Flag: 0x%lx.\n", fw_status);

  /* Print sampling info:                                                   */
  fw(outsample, !=0, &transit);
  t0 = timecheck(verblevel, itr,  7, "outsample", tv, t0);

  /* Initialize CIA:                                                     */
  fw(readcia, !=0, &transit);
  t0 = timecheck(verblevel, itr,  8, "readcia", tv, t0);
  fw(interpolatecia, !=0, &transit);
  t0 = timecheck(verblevel, itr,  9, "interpolatecia", tv, t0);

  /* Compute index of refraction:                                        */
  fw(idxrefrac, !=0, &transit);
  t0 = timecheck(verblevel, itr,  10, "idxrefrac", tv, t0);
 
  /* Calculate extinction coefficient:                                   */
  fw(extwn, !=0, &transit);
  t0 = timecheck(verblevel, itr, 11, "extwn", tv, t0);
 
  /* Initialize structures for the optical-depth calculation:            */
  fw(init_optdepth, !=0, &transit);
  /* Calculate optical depth for eclipse:                                */
  if (strcmp(transit.sol->name, "eclipse") == 0){
    transitprint(1,verblevel, "\nCalculating eclipse:\n");

    /* Calculate intensity for each incident angle:                      */
    for(int i=0; i < transit.ann; i++){
      /* Fills out angle index                                           */
      transit.angleIndex = i;

      fw(tau, !=0, &transit);
      t0 = timecheck(verblevel, itr, 12, "tau eclipse", tv, t0);
  
      /* Calculates eclipse intensity:                                   */
      /* In cgs units erg/s/sr/cm                                        */
      fw(emergent_intens, !=0, &transit);
      t0 = timecheck(verblevel, itr, 13, "emergent intensity", tv, t0);
    }

    /* Calculates flux  erg/s/cm                                         */
    fw(flux, !=0, &transit);
    t0 = timecheck(verblevel, itr, 14, "flux", tv, t0);

    /* Free no longer needed memory                                      */
    freemem_intensityGrid(transit.ds.intens, &transit.pi);
  }

  /* Calculate optical depth for transit:                                */
  else if (strcmp(transit.sol->name, "transit") == 0){
    transitprint(1,verblevel, "\nCalculating transit:\n");
    fw(tau, !=0, &transit);
    t0 = timecheck(verblevel, itr, 12, "tau transit", tv, t0); 

   /* Calculates transit modulation:                                     */
    fw(modulation, !=0, &transit);
    t0 = timecheck(verblevel, itr, 13, "modulation", tv, t0);
  }
 
  /* Free no longer needed memory                                        */
  freemem_idexrefrac(transit.ds.ir,        &transit.pi);
  freemem_extinction(transit.ds.ex,        &transit.pi);
  freemem_tau(transit.ds.tau,              &transit.pi);

  free(transit.save.ext);
  freemem_cia      (transit.ds.cia, &transit.pi);
  freemem_outputray(transit.ds.out, &transit.pi);
  t0 = timecheck(verblevel, itr, 14, "THE END", tv, t0);
  transitprint(1, verblevel, "----------------------------\n");

  freemem_isotopes(     transit.ds.iso, &transit.pi);
  freemem_molecules(    transit.ds.mol, &transit.pi);
  freemem_atmosphere(   transit.ds.at,  &transit.pi);
  freemem_lineinfotrans(transit.ds.li,  &transit.pi);
  freemem_transit(&transit);

  return EXIT_SUCCESS;
}


/* \fcnfh
   Frees transit structure.                 */
void
freemem_transit(struct transit *tr){
  freemem_hints(tr->ds.th);

  freemem_samp(&tr->rads);
  freemem_samp(&tr->wns);
  freemem_samp(&tr->ips);
  free_atm(&tr->atm);

  free(tr->outpret);
  /* TBD: Free saves once it is enabled
  freemem_saves();                          */
}
