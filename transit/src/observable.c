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

/* Revision        April 26th,   2014 Jasmina Blecic
                   moved lines that free memory to transit.c                */

#include <transit.h>

/* \fcnfh
   Calculate the transit modulation at each wavenumber
   Return: 0 on success, else
          -1 if impact parameter sampling is not equispaced                 */
int
modulation(struct transit *tr){
  struct optdepth *tau = tr->ds.tau;
  struct geometry *sg  = tr->ds.sg;
  static struct outputray st_out;
  tr->ds.out = &st_out;

  long w;
  prop_samp *ip = &tr->ips;
  prop_samp *wn = &tr->wns;
  ray_solution *sol = tr->sol;

  /* Check that impact parameter and wavenumber samples exist:              */
  transitcheckcalled(tr->pi, "modulation", 3, "tau",          TRPI_TAU,
                                              "makeipsample", TRPI_MAKEIP,
                                              "makewnsample", TRPI_MAKEWN);

  /* Allocate the modulation array:                                         */
  PREC_RES *out = st_out.o = (PREC_RES *)calloc(wn->n, sizeof(PREC_RES));

  /* Set time to the user hinted default, and other user hints:             */
  setgeom(sg, HUGE_VAL, &tr->pi);

  /* Integrate for each wavelength:                                         */
  transitprint(1, verblevel, "Integrating over wavelength.\n");

  int nextw = wn->n/10;

  /* Calculate the modulation spectrum at each wavenumber:                  */
  for(w=0; w < wn->n; w++){
    out[w] = sol->spectrum(tr, tau->t[w], wn->v[w], tau->last[w],
                           tau->toomuch, ip);
    if (out[w] < 0){
      switch(-(int)out[w]){
      case 1:
        if(tr->modlevel == -1)
          transiterror(TERR_SERIOUS, "Optical depth didn't reach limiting "
                       "%g at wavenumber %g cm-1 (only reached %g).  Cannot "
                       "use critical radius technique (-1).\n", tau->toomuch,
                       tau->t[w][tau->last[w]], wn->v[w]*wn->fct);
      default:
        transiterror(TERR_SERIOUS, "There was a problem while calculating "
                     "modulation at wavenumber %g cm-1. Error code %i.\n",
                     wn->v[w]*wn->fct, (int)out[w]);
        break;
      }
      exit(EXIT_FAILURE);
    }

    /* Print to screen the progress status:                                 */
    if(w==nextw){
      nextw += wn->n/10;
      transitprint(2, verblevel, "%i%% ", (10*(int)(10*w/wn->n+0.9999999999)));
    }
  }
  transitprint(1, verblevel, "\nDone.\n");

  /* Set progress indicator, and print output:                              */
  tr->pi |= TRPI_MODULATION;
  printmod(tr);
  return 0;
}


/* \fcnfh
   Print (to file or stdout) the modulation as function of wavelength */
void
printmod(struct transit *tr){
  FILE *outf = stdout;
  struct outputray *outray = tr->ds.out;
  int rn;

  /* Open file: */
  if(tr->f_outmod && tr->f_outmod[0] != '-')
    outf = fopen(tr->f_outmod, "w");

  transitprint(1, verblevel,
               "\nPrinting in-transit/out-transit modulation in '%s'.\n",
               tr->f_outmod?tr->f_outmod:"standard output");

  /* Print: */
  char wlu[20], /* Wavelength units name */
       wnu[20]; /* Wavenumber units name (the inverse, actually) */
  long nsd = (long)(1e6);

  /* Get wavenumber units name: */
  if((long)(nsd*tr->wns.fct)==nsd) strcpy(wnu, "cm");
  else if((long)(nsd*1e-1*tr->wns.fct)==nsd) strcpy(wnu, "mm");
  else if((long)(nsd*1e-4*tr->wns.fct)==nsd) strcpy(wnu, "um");
  else if((long)(nsd*1e-7*tr->wns.fct)==nsd) strcpy(wnu, "nm");
  else if((long)(nsd*1e-8*tr->wns.fct)==nsd) strcpy(wnu, "A ");
  else sprintf(wnu, "%6.1g cm", 1/tr->wns.fct);

  /* Get wavelength units name: */
  if((long)(nsd*tr->wavs.fct)==nsd) strcpy(wlu, "cm");
  else if((long)(nsd*1e1*tr->wavs.fct)==nsd) strcpy(wlu, "mm");
  else if((long)(nsd*1e4*tr->wavs.fct)==nsd) strcpy(wlu, "um");
  else if((long)(nsd*1e7*tr->wavs.fct)==nsd) strcpy(wlu, "nm");
  else if((long)(nsd*1e8*tr->wavs.fct)==nsd) strcpy(wlu, "A ");
  else sprintf(wlu, "%8.1g cm", tr->wavs.fct);

  /* Print header: */
  fprintf(outf, "#wvl [um]        modulation\n");

  /* Print wavelength (in microns) and modulation at each wavenumber:       */
  for(rn=0; rn<tr->wns.n; rn++)
    fprintf(outf, "%-17.9g%-18.9g\n",
                  1/(tr->wns.v[rn]/tr->wns.fct*1e-4),
                  outray->o[rn]);

  fclose(outf);
  return;
}


/*\fcnfh
  Free the transit modulation array

  Return: 0 on success                          */
int
freemem_outputray(struct outputray *out,
                  long *pi){
  /* Free arrays: */
  free(out->o);

  /* Clear PI and return: */
  *pi &= ~(TRPI_MODULATION);
  return 0;
}
