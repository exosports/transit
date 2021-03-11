// Copyright (C) 2015-2016 University of Central Florida. All rights reserved.
// Transit is under an open-source, reproducible-research license (see LICENSE).

#include <transit.h>


/* Initial version January 23rd, 2014 Jasmina Blecic
                   implemented eclipse                                     */
/* Revision        March 19th,   2014 Jasmina Blecic
                   implemented switch eclipse/transit                      */
/* Revision        April 26th,   2014 Jasmina Blecic
                   implemented intensity grid and flux                     */


/* Defines static variables                                                 */
static PREC_RES *area_grid;

/* #########################################################
    CALCULATES OPTICAL DEPTH AT VARIOUS POINTS ON THE PLANET
   ######################################################### */


/* FUNCTION
   Computes optical depth for eclipse geometry for one ray, one wn,
   at various incident angles on the planet surface,
   between a certain layer in the atmosphere up to the top layer.
   Returns: Optical depth divided by rad.fct:  \frac{tau}{units_{rad}}      */
static PREC_RES
eclipsetau(struct transit *tr,
           PREC_RES height,    /* Altitude down to where calculate tau      */
           PREC_RES *ex){      /* Extinction per layer [rad]                */

  /* Layers radius array:                                                   */
  prop_samp *rads = &tr->rads;  /* Radius sampling                          */
  PREC_RES *rad  = rads->v;     /* Radius array                             */

  /* Get the index rs, of the sampled radius immediately below or equal
     to height (i.e. rad[rs] <= height < rad[rs+1]):                        */
  int rs = binsearchapprox(rad, height, 0, tr->rads.n-1);

  /* Auxiliary variables for Simson integration:                            */
  double *hsum, *hratio, *hfactor, *h;

  /* Returns 0 if this is the top layer (no distance travelled):            */
  if (rs == tr->rads.n-1)
    return 0.0;

  /* Move pointers to the location of height:                               */
  rad += rs;
  ex  += rs;

  /* Number of layers beween height and the top layer:                      */
  int nrad = tr->rads.n - rs;

  PREC_RES res;          /* Optical depth divided by units of radius        */
  PREC_RES x3[3], r3[3]; /* Interpolation variables                         */

  /* Distance along the path:                                               */
  PREC_RES s[nrad];

  /* Providing three necessary points for spline integration:               */
  const PREC_RES tmpex  = *ex;
  const PREC_RES tmprad = *rad;

  if(nrad==2) *ex = interp_parab(rad-1, ex-1, rad[0]);
  else        *ex = interp_parab(rad,   ex,   rad[0]);

  if(nrad==2){
    x3[0] = ex[0];
    x3[2] = ex[1];
    x3[1] = (ex[1]+ex[0])/2.0;
    r3[0] = rad[0];
    r3[2] = rad[1];
    r3[1] = (rad[0]+rad[1])/2.0;
    *rad = tmprad;
    *ex  = tmpex;
    rad  = r3;
    ex   = x3;
    nrad++;
  }

  /* Distance along the path:                                               */
  s[0] = 0.0;
  for(int i=1; i < nrad; i++){
    s[i] = s[i-1] + (rad[i] - rad[i-1]);
  }

  hsum    = calloc((nrad-1)/2, sizeof(double));
  hratio  = calloc((nrad-1)/2, sizeof(double));
  hfactor = calloc((nrad-1)/2, sizeof(double));
  h       = calloc(nrad-1, sizeof(double));

  /* Integrate extinction along the path:                                   */
  makeh(s, h, nrad);
  geth(h, hsum, hratio, hfactor, nrad);
  res = simps(ex, h, hsum, hratio, hfactor, nrad);

  free(hsum);
  free(hratio);
  free(hfactor);
  free(h);

  /* Optical depth divided by units of radius:                              */
  return res;
}


/* #################################################
    CALCULATES EMERGENT INTENSITY FOR ONE WAVENUMBER
   ################################################# */

/* \fcnfh
   Calculates emergent intensity.
   Return: emergent intensity for one wavenumber                            */

/* DEF */
static PREC_RES
eclipse_intens(struct transit *tr,  /* Transit structure                    */
               PREC_RES *tau,       /* Optical depth array                  */
               PREC_RES w,          /* Current wavenumber value             */
               long last,           /* Index where tau == toomuch           */
               double toomuch,      /* Maximum optical depth calculated     */
               prop_samp *rad){     /* Radius array                         */
  /* FINDME: toomuch is not needed as a parameter                           */

  /* General variables:                                                     */
  PREC_RES res;                  /* Result                                  */
  PREC_ATM *temp = tr->atm.t;    /* Temperatures                            */

  PREC_RES angle = tr->angles[tr->angleIndex] * DEGREES;

  /* Takes sampling properties for wavenumber from tr:                      */
  prop_samp *wn = &tr->wns;
  /* Wavenumber units factor to cgs:                                        */
  double wfct  = wn->fct;

  /* Radius parameter variables:                                            */
  long rnn  = rad->n;
  long i;

  /* Blackbody function at each layer:                                      */
  PREC_RES B[rnn];
  /* Integration parts:                                                     */
  PREC_RES dtau[rnn];      /* Tau integration variable                      */

  /* Integrate for each of the planet's layer starting from the
     outermost until the closest layer.                                     */

  /* Planck function (erg/s/sr/cm) for wavenumbers:
        B_\nu = 2 h {\bar\nu}^3 c^2 \frac{1}
                {\exp(\frac{h \bar \nu c}{k_B T})-1}                        */
  for(i=0; i <= last; i++){
    dtau[i] = exp(-tau[i]/cos(angle));
    B[i] = (2.0 * H * pow(w*wfct,3.0) * LS * LS) /
            (exp(H * w*wfct * LS / (KB * temp[rnn-1-i])) - 1.0);
  }
  /*    Background emission, medium emission                                */
  res = B[last]*dtau[last] - integ_trapz(dtau, B, last+1);
  return res;
}


/* ###############################################################
    CALCULATES EMERGENT INTENSITY AT VARIOUS POINTS ON THE PLANET
   ############################################################### */

/* \fcnfh
   Calculates the emergent intensity (ergs/s/sr/cm) for the whole range
   of wavenumbers at the various points on the planet
   Returns: emergent intensity for the whole wavenumber range               */
/* DEF */
int
emergent_intens(struct transit *tr){  /* Transit structure                  */
  static struct outputray st_out;     /* Output structure                   */
  tr->ds.out = &st_out;

  /* Initial variables:                                                     */
  long w;
  prop_samp *rad = &tr->rads;          /* Radius array pointer              */
  prop_samp *wn  = &tr->wns;           /* Wavenumber array pointer          */
  long int wnn   = wn->n;              /* Wavenumbers                       */
  ray_solution *sol = tr->sol;         /* Eclipse ray solution pointer      */

  /* Reads angle index from transit structure                               */
  long int angleIndex = tr->angleIndex;
  /* Intensity for all angles and all wn                                    */
  PREC_RES **intens_grid = tr->ds.intens->a;

  /* Intensity array for one angle all wn:                                  */
  PREC_RES *out = intens_grid[angleIndex];

  /* Reads the tau array from transit structure                             */
  struct optdepth *tau = tr->ds.tau;

  /* Integrate for each wavelength:                                         */
  tr_output(TOUT_RESULT, "Integrating over wavelength.\n");

  /* Printing process variable:                                             */
  int nextw = wn->n/10;

  /* Calculates the intensity integral at each wavenumber:                  */
  for(w=0; w<wnn; w++){
    //transitprint(1, 2, "[%li]", w);
    //if (w == 1612 || w == 1607){
    //  transitprint(1, 2, "\nTau (%.3f) [%li,%li]= np.array([", wn->v[w],
    //                      rad->n, tau->last[w]);
    //  //for (int ii=rad->n-1; ii>tau->last[w]; ii--)
    //  for (int ii=0; ii<rad->n; ii++)
    //    transitprint(1, 2, "%.4e, ", tau->t[w][ii]);
    //  transitprint(1, 2, "])\n");
    //}
    //if (fabs(wn->v[w] - 1844.59) < 0.005)
    //  transitprint(1, verblevel, "\nWavenumber index is: %li\n", w);

    /* Calculate the intensity spectrum (call to eclipse_intens):           */
    out[w] = sol->spectrum(tr, tau->t[w], wn->v[w], tau->last[w],
                           tau->toomuch, rad);

    /* Prints to screen the progress status:                                */
    if(w == nextw){
      nextw += wn->n/10;
      tr_output(TOUT_DEBUG, "%i%% ", (10*(int)(10*w/wn->n+0.9999999999)));
    }
  }
  tr_output(TOUT_DEBUG, "\n");
  tr_output(TOUT_RESULT, "Done.\n");

  /* Sets progress indicator, and prints output:                             */
  tr->pi |= TRPI_MODULATION; /* FINDME: this is not a modulation calculation */
  if (tr->angleIndex == tr->ann-1)
    printintens(tr);
  return 0;
}


/* FUNCTION
   Calculate the flux spectrum
   Formula:
   Flux = pi * SUMM_i [I_i * (sin(theta_fin)^2 - sin(theta_in)^2)]
   I_i are calculated for each angle defined in the configuration file
   Returns: zero on success                                                 */
int
flux(struct transit *tr){  /* Transit structure                             */
  static struct outputray st_out;
  tr->ds.out = &st_out;

  /* Get angles and number of angles from transithint:                      */
  PREC_RES *angles = tr->angles;      /* Angles                             */
  long int an = tr->ann;              /* Number of angles                   */

  /* Intensity for all angles and all wn                                    */
  PREC_RES **intens_grid = tr->ds.intens->a;

  long int i, w;  /* for-loop indices                                       */
  PREC_RES area,  /* Projected area                                         */
           *out;  /* Output flux array (per wavenumber)                     */

  prop_samp *wn  = &tr->wns; /* Wavenumber sample                           */
  long int wnn = wn->n;      /* Number of wavenumbers                       */

  /* Allocate area grid and set its first and last value:                   */
  area_grid = (PREC_RES *)calloc(an+1, sizeof(PREC_RES));
  area_grid[ 0] =  0.0 * DEGREES;
  area_grid[an] = 90.0 * DEGREES;

  /* Fills out area grid array. Converts to radians.
     Limits of each area defined in the middle of the angles given:         */
  for(i = 1; i < an; i++)
    area_grid[i] = (angles[i-1] + angles[i]) * DEGREES / 2.0;

  /* Allocates array for the emergent flux:                                 */
  out = st_out.o = (PREC_RES *)calloc(wnn, sizeof(PREC_RES));

  /* Add weighted Intensity to get the flux:                                */
  for(i = 0; i < an; i++){
    area = pow(sin(area_grid[i+1]), 2.0) - pow(sin(area_grid[i]), 2.0);
    for(w=0; w < wnn; w++)
      out[w] += PI * intens_grid[i][w] * area;
  }

  /* Free memory that is no longer needed                                   */
  freemem_localeclipse();

  /* prints output                                                          */
  printflux(tr);
  return 0;
}


/* \fcnfh
   Print (to file or stdout) the emergent intensities as function of wavelength
   for each angle)                                                          */
void
printintens(struct transit *tr){
  long int i, w;                  /* Auxilliary for-loop indices            */
  FILE *outf = stdout;            /* Output file pointer                    */

  PREC_RES *angles = tr->angles;  /* Array of incident angles               */
  int an = tr->ann;               /* Number of incident angles              */

  prop_samp *wn = &tr->wns;       /* Wavenumber sample                      */
  long int wnn = wn->n;           /* Number of wavenumber samples           */

  /* Intensity for all angles and all wn:                                   */
  PREC_RES **intens_grid = tr->ds.intens->a;

  /* Adds string to the output files to differentiate between outputs        */
  char our_fileName[512];

  /* Open file:                                                             */
  if (tr->f_outintens && tr->f_outintens[0] != '-'){
    strncpy(our_fileName, tr->f_outintens, 512);
    outf = fopen(our_fileName, "w");
  }
  else{
    tr_output(TOUT_INFO, "No intensity file requested.\n");
    return;
  }
  tr_output(TOUT_INFO, "\nPrinting intensity in '%s'\n",
                             tr->f_outintens ? our_fileName:"standard output");

  /* Print the header:                                                      */
  //fprintf(outf, "#wvl [um]%*s", 6, " ");
  //for(i=0; i < an; i++)
  //    fprintf(outf, "I[%4.1lf deg]%*s", angles[i], 7, " ");
  //fprintf(outf, "[erg/s/cm/sr]\n");

  /* Print the header:                                                      */
  fprintf(outf, "#wvl %*s", 10, " ");
  for(i=0; i < an; i++)
      fprintf(outf, "I[%4.1lf deg]%*s", angles[i], 7, " ");
  fprintf(outf, "\n#[um]%*s", 10, " ");
  for(i=0; i < an; i++)
      fprintf(outf, "[erg/s/cm/sr]%*s", 5, " ");
  fprintf(outf, "\n");


  /* Fills out each column with the correct output intensity                 */
  for(w=0; w<wnn; w++){
    fprintf(outf, "%-15.10g", 1e4/(tr->wns.v[w]/tr->wns.fct));
    for(i=0; i < an; i++)
      fprintf(outf, "%-18.9g", intens_grid[i][w]);
    fprintf(outf,"\n");
  }

  /* Closes the file:                                                        */
  fclose(outf);
  return;
}


/* \fcnfh
   Print (to file or stdout) the emergent intensity as function of wavenumber
   (and wavelength)                                                         */
void
printflux(struct transit *tr){
  FILE *outf=stdout;
  /* The flux per wavenumber array:                                         */
  PREC_RES *Flux = tr->ds.out->o;
  int rn;

  /* Open file:                                                             */
  if(tr->f_outspec && tr->f_outspec[0] != '-')
    outf = fopen(tr->f_outspec, "w");

  tr_output(TOUT_INFO, "\nPrinting day-side emission flux spectrum in '%s'\n",
            tr->f_outspec? tr->f_outspec:"standard output");

  /* Print the header:                                                      */
  fprintf(outf, "#wvl [um]%*sFlux [erg/s/cm]\n", 6, " ");

  /* Print wavelength and flux:                                             */
  for(rn=0; rn < tr->wns.n; rn++)
    fprintf(outf, "%-15.10g%-18.9g\n", 1e4/(tr->wns.v[rn]/tr->wns.fct),
            Flux[rn]);

  /* Closes the file:                                                       */
  fclose(outf);
  return;
}


/* \fcnfh
   Frees eclipse pointer arrays. Data array should already be free          */
void
freemem_localeclipse(){
  /* Free auxiliar variables:                                               */
  free(area_grid);
}


/* \fcnfh
   Free intensity grid structure arrays
   Return 0 on success                                                      */
int
freemem_intensityGrid(struct grid *intens,   /* grid structure              */
                      long *pi){             /* progress indicator flag     */

  /* Free arrays:                                                           */
  free(intens->a[0]);
  free(intens->a);

  /* Update indicator and return:                                           */
  *pi &= ~(TRPI_GRID);
  return 0;
}

const ray_solution eclipsepath = {
  "eclipse",         /* Name of the solution                                */
  "eclipse.c",       /* Source code file name                               */
  0,                 /* Request equispaced layer sampling                   */
  &eclipsetau,       /* Optical-depth calculation function                  */
  &eclipse_intens,   /* Intensity calculation function                      */
};
