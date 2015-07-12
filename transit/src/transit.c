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

Copyright (C) 2015 University of Central Florida.  All rights reserved.

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

/* TBD: calloc checks */

#include <transit.h>
struct transit transit;
long itr=0;
struct timeval tv;
double t0=0.0;
int    init_run=0;

void transit_init(int argc, char **argv);
int  get_no_samples(void);
void get_waveno_arr(double *waveno_arr, int waveno);
void set_radius(double refradius);
void run_transit(double *re_input, int transint, double *transit_out,
                 int transit_out_size);
void do_transit(double *transit_out);


void transit_init(int argc, char **argv){
  /* The purpose of this function is to set up and initialize all the
     structures nessisary to run the transit code.                          */
  memset(&transit, 0, sizeof(struct transit));
  verblevel=2;

  /* Process the command line arguments:                                    */
  fw(processparameters, !=0, argc, argv, &transit);
  t0 = timecheck(verblevel, itr,  0, "processparameters", tv, t0);

  transitprint(1, verblevel, "verblevel: %d\n", verblevel);

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

  /* Calculate opacity grid:                                                */
  fw(opacity, <0, &transit);
  t0 = timecheck(verblevel, itr,  5, "opacity", tv, t0);

  /* Initialize CIA:                                                        */
  fw(readcia, !=0, &transit);
  t0 = timecheck(verblevel, itr,  6, "readcia", tv, t0);
  init_run = 1;
}


int get_no_samples(void){
  /* This function will return the size of the wave number array */
  return (int)transit.wns.n;
}

void get_waveno_arr(double * waveno_arr, int waveno){
  int i;
  if (init_run > 0){
    for(i=0; i < (int)transit.wns.n; i++){
      waveno_arr[i] = transit.wns.v[i];
    }
  }
  else{
    printf("Transit not initialized, please run init. Values set -1\n");
    for(i=0; i < (int)transit.wns.n; i++){
        waveno_arr[i] = -1;
    }
  }
}


void set_radius(double refradius){
  transit.r0 = refradius;
}


void run_transit(double *re_input, int transtint, double *transit_out,
                 int transit_out_size){
  fw(reloadatm, <0, &transit, re_input);
  do_transit(transit_out);
}


void do_transit(double * transit_out){
  int i;
  if (init_run > 0){
    fw(makeipsample, <0, &transit);
    t0 = timecheck(verblevel, itr,  6, "makeipsample", tv, t0);
    if(fw_status>0)
      transitprint(7, verblevel, "makeipsample() modified some of the hinted "
                                 "parameters. Flag: 0x%lx.\n", fw_status);

    /* Interpolate CIA:                                                     */
    fw(interpolatecia, !=0, &transit);
    t0 = timecheck(verblevel, itr,  9, "interpolatecia", tv, t0);

    /* Compute index of refraction:                                         */
    fw(idxrefrac, !=0, &transit);
    t0 = timecheck(verblevel, itr,  10, "idxrefrac", tv, t0);

    /* Calculate extinction coefficient:                                    */
    fw(extwn, !=0, &transit);
    t0 = timecheck(verblevel, itr, 11, "extwn", tv, t0);

    /* Initialize structures for the optical-depth calculation:            */
    fw(init_optdepth, !=0, &transit);

    /* Calculates optical depth for eclipse                                 */
    if(strcmp(transit.sol->name, "eclipse") == 0){
      transitprint(1, verblevel, "\nCalculating eclipse:\n");

      fw(tau, !=0, &transit);
      t0 = timecheck(verblevel, itr, 12, "tau eclipse", tv, t0);

      /* Calculate optical depth for eclipse:                               */
      for(i=0; i < transit.ann; i++){
        /* Fills out angle index                                            */
        transit.angleIndex = i;

        /* Calculates eclipse intensity:                                    */
        /* In cgs units erg/s/sr/cm                                         */
        fw(emergent_intens, !=0, &transit);
        t0 = timecheck(verblevel, itr, 13, "emergent intensity", tv, t0);
      }

      /* Calculates flux  erg/s/cm                                          */
      fw(flux, !=0, &transit);
      t0 = timecheck(verblevel, itr, 14, "flux", tv, t0);

      /* Free no longer needed memory                                       */
      freemem_intensityGrid(transit.ds.intens, &transit.pi);
    }

    /* Calculate optical depth for transit:                                 */
    else if (strcmp(transit.sol->name, "transit") == 0){
      transitprint(1, verblevel, "\nCalculating transit:\n");
      fw(tau, !=0, &transit);
      t0 = timecheck(verblevel, itr, 12, "tau transit", tv, t0);

      /* Calculates transit modulation:                                     */
      fw(modulation, !=0, &transit);
      t0 = timecheck(verblevel, itr, 13, "modulation", tv, t0);
   }

    for(int i=0; i < transit.wns.n; i++){
      transit_out[i] = transit.ds.out->o[i];
    }

    /* Free arrays allocated inside the individual call:                    */
    free(transit.save.ext);
    freemem_samp(&transit.ips);
    freemem_idexrefrac(transit.ds.ir,  &transit.pi);
    freemem_extinction(transit.ds.ex,  &transit.pi);
    freemem_tau(       transit.ds.tau, &transit.pi);
    freemem_outputray( transit.ds.out, &transit.pi);

    t0 = timecheck(verblevel, itr, 14, "THE END", tv, t0);
    transitprint(1, verblevel, "----------------------------\n");
    itr++;
  }
  else{
    printf("Transit init not run, please initialize transit\n");
  }
}

void free_memory(void){
  /* This function frees all the memory used in transit, and should be
     called at the end of the program. Check if all these data structures
     can be used when called from bart, In MPItransit not all the frees
     happen.                                                                */
  freemem_molecules( transit.ds.mol, &transit.pi);
  freemem_atmosphere(transit.ds.at,  &transit.pi);
  if (transit.fp_opa == NULL)
    freemem_linetransition(&transit.ds.li->lt,  &transit.pi);
  freemem_lineinfo(  transit.ds.li,  &transit.pi);
  freemem_cia(       transit.ds.cia, &transit.pi);
  freemem_transit(&transit);
  init_run = 0;
}


int main(int argc, char **argv){
  transit_init(argc, argv);
  int trans_size = get_no_samples();
  double tmp[trans_size];
  do_transit(tmp);
  free_memory();
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
