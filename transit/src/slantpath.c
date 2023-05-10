// Copyright (C) 2015-2016 University of Central Florida. All rights reserved.
// Transit is under an open-source, reproducible-research license (see LICENSE).

#include <transit.h>


/* FUNCTION
   Compute the light path and optical depth at a given impact parameter
   and wavenumber, for a medium with constant index of refraction

   The impact parameter, b, and radii rad should have the same units.
   The result needs to be multiplied by the units of 'rad' to be real.
   It can take nrad values of 1 or bigger. However, if 2 is given then
   'ex' and 'rad' need to have a referenceable element at position
   -1; i.e. rad[-1] and ex[-1] need to exist.

   Return: Optical depth divided by rad.fct:  \frac{tau}{units_{rad}}       */
static PREC_RES
totaltau1(PREC_RES b,    /* Impact parameter                                */
          PREC_RES *rad, /* Layers radius array                             */
          PREC_RES refr, /* Refractivity index                              */
          PREC_RES *ex,  /* Extinction[rad]                                 */
          long nrad){    /* Number of radii elements                        */

  int rs, i;             /* Auxilliary radius and for-loop indices          */
  PREC_RES res;          /* Half-path optical depth                         */
  PREC_RES x3[3], r3[3]; /* Auxiliary interpolation variables               */

  double *hsum, *hratio, *hfactor, *h;

  /* Closest approach radius:                                               */
  PREC_RES r0 = b/refr;

  /* Get the index rs, of the sampled radius immediately below or equal
     to r0 (i.e. rad[rs] <= r0 < rad[rs+1]):                                */
  rs = binsearch(rad, 0, nrad-1, r0);
  if ((rs == -5) || (rs == -2))
    return 0;  /* If it is the outmost layer                                */
  /* If some other error occurred:                                          */
  else if(rs < 0){
    tr_output(TOUT_ERROR, "Closest approach (%g) is larger than the "
      "top layer of the atmosphere (%g).\n", r0, rad[nrad-1]);
    exit(EXIT_FAILURE);
  }

  /* Move extinction and radius pointers to the rs-th element:              */
  rad  += rs;
  ex   += rs;
  nrad -= rs;

  /* Calculate the extinction coefficient at r0 using a parabolic
     interpolation.  Temporarily store extinction(r0) and r0 at the rs-th
     position of ex and rad arrays:                                         */
  const PREC_RES tmpex  = *ex;  /* Get a constant copy of the arrays        */
  const PREC_RES tmprad = *rad;
  /* Interpolate and modify arrays values:                                  */
  if(nrad==2) *ex = interp_parab(rad-1, ex-1, r0);
  else        *ex = interp_parab(rad,   ex,   r0);
  *rad = r0;
  /* Handle special case of only two elements in array, make an
     intermediate point in between values:                                  */
  if(nrad == 2){
    x3[0] = ex[0];
    x3[2] = ex[1];
    x3[1] = (ex[0] + ex[1])/2.0;
    r3[0] = rad[0];
    r3[2] = rad[1];
    r3[1] = (rad[0] + rad[1])/2.0;
    *rad = tmprad;
    *ex  = tmpex;
    rad  = r3;
    ex   = x3;
    nrad++;
  }

  /* Convert to s spacing, i.e. distance along the path (between the
     closest-approach point and the i-th layer).
     Radius needs to be equispaced:                                         */
  PREC_RES s[nrad];  /* Distance along ray path                             */
  s[0] = 0;

  /* Only valid for refraction index = 1.0                                  */
  for(i=1; i<nrad; i++){
    s[i] = sqrt(rad[i]*rad[i] - r0*r0);
  }

  /* Integrate Extinction along ray path:                                   */
  hsum    = calloc((nrad-1)/2, sizeof(double));
  hratio  = calloc((nrad-1)/2, sizeof(double));
  hfactor = calloc((nrad-1)/2, sizeof(double));
  h       = calloc(nrad-1, sizeof(double));

  makeh(s, h, nrad);
  geth(h, hsum, hratio, hfactor, nrad);
  res = simps(ex, h, hsum, hratio, hfactor, nrad);

  /* Reset original values of arrays:                                       */
  *ex  = tmpex;
  *rad = tmprad;

  /* Return:                                                                */
  free(hsum);
  free(hratio);
  free(hfactor);
  free(h);

  return 2*res;
}


/* FUNCTION
   Compute the light path and optical depth at a given impact parameter
   and wavenumber, for a medium with constant index of refraction

   The impact parameter b needs to be given in units of 'rad' and the
   result needs to be multiplied by the units 'rad' to be real.

   Return: $\frac{tau}{units_{rad}}$ returns optical depth divided by
           units of 'rad'                                                   */
static PREC_RES
totaltau2(PREC_RES b,      /* Impact parameter                              */
          PREC_RES *rad,   /* Radius array                                  */
          PREC_RES *refr,  /* Refractivity index                            */
          PREC_RES *ex,    /* Extinction[rad]                               */
          long nrad){      /* Number of radii elements                      */
  PREC_RES dt[nrad];
  PREC_RES r0a = b;
  PREC_RES r0 = 0;
  int rs, i=0;
  const int maxiterations=50;

  /* Auxiliary variables for Simson integration:                            */
  double *hsum, *hratio, *hfactor, *h;

  tr_output(TOUT_ERROR, "This routine has not been implemented yet.\n");
  exit(EXIT_FAILURE);

  /* Look for closest approach radius:                                      */
  while(1){
    r0 = b/lineinterp(r0a, rad, refr, nrad);
    if(r0==r0a)
      break;
    if(i++ > maxiterations) {
      tr_output(TOUT_ERROR, "Maximum iterations (%i) reached while "
        "looking for r0. Convergence not reached (%.6g!=%.6g).\n",
        maxiterations, r0, r0a);
      exit(EXIT_FAILURE);
    }
    r0a = r0;
  }

  /* Get bin value 'rs' such that r0 is between rad[rs-1] inclusive
     and rad[rs] exclusive:                                                 */
  /* If we are looking at the outmost layer, then return:                   */
  if((rs=binsearch(rad, 0, nrad-1, r0))==-5)
    return 0;

  /* If some other error occurred:                                          */
  else if(rs<0) {
    tr_output(TOUT_ERROR,
      "Closest approach value(%g) is outside sampled radius "
      "range(%g - %g).\n", r0, rad[0], rad[nrad-1]);
    exit(EXIT_FAILURE);
  }
  /* Move pointer to rs-th element:                                         */
  rs++;
  /* FINDME: This will break                                                */

  /* A fraction 'analiticfrac' of the integration near the closest
     appraoach is calcualated analitically, otherwise, I get a division
     by zero in formula:
     \tau_{\nu}(\rho) =
     \underbrace{\frac{2e_{\nu}\rho}{n}
                  \left(\sqrt{\left(\frac{nr_1}{\rho}\right)^2-1} -
                        \sqrt{\left(\frac{nr_0}{\rho}\right)^2-1}\right)
     }_{\rm{analitic}} +
     \underbrace{2\int_{r_1=r_0+\delta r}^{\infty}
                 \frac{e_{\nu}~n~r}{\sqrt{n^2r^2-\rho^2}}{\rm d} r
     }_{\rm{numerical}}                                                     */

  /* First for the analitical part of the integral:                         */
  PREC_RES res;
  if(ex[rs-1]==ex[rs])
    res = ex[rs] * r0 * (sqrt( rad[rs]*rad[rs] / (r0*r0) - 1));
  else{
    PREC_RES alpha = ( ex[rs] - ex[rs-1] ) / ( rad[rs] - rad[rs-1] );
    PREC_RES rm    = rad[rs];
    if(alpha<0)
      res = - alpha * (rm * sqrt( rm*rm -  r0*r0) - r0*r0 *
                        log(sqrt( rm*rm / (r0*r0) - 1) + rm/r0 )) / 2.0;
    else
      res =   alpha * (rm * sqrt( rm*rm -  r0*r0) + r0*r0 *
                        log(sqrt( rm*rm / (r0*r0) - 1) + rm/r0 ))  / 2.0;
  }

  /* And now for the numerical integration.  Set the variables:             */
  for(i=rs; i<nrad; i++){
    r0a = b / (refr[i]*rad[i]);
    transitASSERT(r0a>1, "Condition could not be asserted, b/(nr)=%g > 1.\n",
                  r0a);

    dt[i] = ex[i]/sqrt(1-r0a*r0a);
  }

  /* Allocate auxillary integration arrays                                  */
  hsum    = calloc((nrad-rs-1)/2, sizeof(double));
  hratio  = calloc((nrad-rs-1)/2, sizeof(double));
  hfactor = calloc((nrad-rs-1)/2, sizeof(double));
  h       = calloc(nrad-rs-1,   sizeof(double));

  /* Allocate array to hold integration range                               */
  double *radtemp;
  double *dttemp;
  radtemp = calloc(nrad-rs, sizeof(double));
  dttemp  = calloc(nrad-rs, sizeof(double));

  for(i=rs;i<nrad;i++){
    radtemp[i-rs] = rad[i];
    dttemp[i-rs]  = dt[i];
  }

  /* Integrate extinction along the path:                                    */
  makeh(radtemp, h, nrad-rs);
  geth(h, hsum, hratio, hfactor, nrad);
  res += simps(dttemp, h, hsum, hratio, hfactor, nrad-rs);

  free(hsum);
  free(hratio);
  free(hfactor);
  free(h);

  return 2*res;
}


/* FUNCTION
 Wrapper function to calculate the optical depth along the path for a
 given impact parameter at a specific wavenumber.

 Return: $\frac{tau}{units_{rad}}$ returns optical depth divided by units
                                    of 'rad'                           */
static inline PREC_RES
transittau(struct transit *tr,
           PREC_RES b,      /* Impact paramer                               */
           PREC_RES *ex){   /* Extinction array [rad]                       */

  /* Optical depth calculation level:                                       */
  tr->taulevel = tr->ds.th->taulevel;
  PREC_RES *refr = tr->ds.ir->n;     /* Index-of-refraction array           */
  PREC_RES *rad  = tr->rads.v;       /* Layer radius array                  */
  PREC_RES nrad  = tr->rads.n;       /* Number of layers                    */

  switch(tr->taulevel){
  case 1: /* Constant index of refraction:                                  */
    return totaltau1(b, rad, *refr, ex, nrad);
    break;
  case 2: /* Variable index of refraction:                                  */
    return totaltau2(b, rad,  refr, ex, nrad);
    break;
  default:
    tr_output(TOUT_ERROR,
      "slantpath:: totaltau:: Level %i of detail has not been "
      "implemented to compute optical depth.\n", tr->taulevel);
    exit(EXIT_FAILURE);
    return 0;
  }
}


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
  tr_output(TOUT_RESULT, "Integrating over wavelength.\n");

  int nextw = wn->n/10;

  /* Calculate the modulation spectrum at each wavenumber:                  */
  for(w=0; w < wn->n; w++){
    out[w] = sol->spectrum(tr, tau->t[w], wn->v[w], tau->last[w],
                           tau->toomuch, ip);
    if (out[w] < 0){
      switch(-(int)out[w]){
      case 1:
        if(tr->modlevel == -1) {
          tr_output(TOUT_ERROR, "Optical depth didn't reach limiting "
            "%g at wavenumber %g cm-1 (only reached %g).  Cannot "
            "use critical radius technique (-1).\n", tau->toomuch,
            tau->t[w][tau->last[w]], wn->v[w]*wn->fct);
          exit(EXIT_FAILURE);
        }
      default:
        tr_output(TOUT_ERROR, "There was a problem while calculating "
          "modulation at wavenumber %g cm-1. Error code %i.\n",
          wn->v[w]*wn->fct, (int)out[w]);
        exit(EXIT_FAILURE);
        break;
      }
      exit(EXIT_FAILURE);
    }

    /* Print to screen the progress status:                                 */
    if(w==nextw){
      nextw += wn->n/10;
      tr_output(TOUT_DEBUG, "%i%% ", (10*(int)(10*w/wn->n+0.9999999999)));
    }
  }
  tr_output(TOUT_DEBUG, "\n");
  tr_output(TOUT_RESULT, "Done.\n");

  /* Set progress indicator, and print output:                              */
  tr->pi |= TRPI_MODULATION;
  printmod(tr);
  return 0;
}


/* FUNCTION
   Calculate the transit modulation at a given wavenumber for no-limb
   darkening nor emitted flux.

   Return: the transit's modulation:
   1 - in-transit/out-of-transit flux ratio (Equation 3.12):
    M_{\lambda} = \frac{1}{R_\star^2}\left(R^2 -
                   2\int_{0}^{R} \exp^{-\tau_\lambda(r)} r\,{\rm d}r\right) */
static PREC_RES
modulation1(PREC_RES *tau,        /* Optical depth array                    */
            long last,            /* Index where tau = toomuch              */
            double toomuch,       /* Maximum optical depth calculated       */
            prop_samp *ip,        /* Impact parameter                       */
            struct geometry *sg){ /* Geometry struct                        */
  PREC_RES res;
  /* Stellar radius:                                                        */
  double srad = sg->starrad * sg->starradfct;

  /* Impact parameter variables:                                            */
  long ipn  = ip->n;
  long ipn1 = ipn-1;
  long i;
  double *hsum, *hratio, *hfactor, *h;

  /* Max overall tau, for the tr.ds.sg.transparent=True case:               */
  const PREC_RES maxtau = tau[last] > toomuch? tau[last]:toomuch;

  PREC_RES rinteg[ipn],  /* Integrand                                       */
           ipv[ipn];     /* Impact parameter where to integrate             */

  /* Integrate for each of the planet's layer starting from the
     outermost until the closest layer:                                     */
  for(i=0; i<=last; i++){
    ipv   [ipn1-i] = ip->v[i] * ip->fct;
    rinteg[ipn1-i] = exp(-tau[i]) * ipv[ipn1-i];
  }
  /* Add one more layer with 0. Only two to have a nice ending
    spline and not unnecessary values:                                      */
  last += 1;
  if (last > ipn1)
    last = ipn1;
  for(; i <= last; i++){
    ipv   [ipn1-i] = ip->v[i] * ip->fct;
    rinteg[ipn1-i] = 0;
  }
  /* FINDME: there is no need to use a for-loop here, the first two lines
     are confusing also. This could have been written much simpler          */

  /* Increment last to represent the number of elements, check that we
     have enough:                                                           */
  last++;
  if(last < 3) {
    tr_output(TOUT_ERROR, "Condition failed, less than 3 items "
      "(only %i) for radial integration.\n", last);
    exit(EXIT_FAILURE);
  }

  /* Integrate along radius:                                                */
  hsum    = calloc((last-1)/2, sizeof(double));
  hratio  = calloc((last-1)/2, sizeof(double));
  hfactor = calloc((last-1)/2, sizeof(double));
  h       = calloc(last-1, sizeof(double));


  makeh(ipv+ipn-last, h, last);
  geth(h, hsum, hratio, hfactor, last);
  res = simps(rinteg+ipn-last, h, hsum, hratio, hfactor, last);

  /* Substract the total area blocked by the planet. This is from the
     following:
     \begin{eqnarray}
     1-M = &1-\frac{\int_0^{r_p}\int e^{-\tau}r{\rm d}\theta {\rm d}r
                  \ +\ \int_{r_p}^{R_s}{\rm d}A} {\pi R_s^2}       \\
         = & -\frac{\int_0^{r_p}\int e^{-\tau}r{\rm d}\theta {\rm d}r
                  \ +\ Area_{p}} {\pi R_s^2}                       \\
         = & -\frac{2\int_0^{r_p} e^{-\tau}r{\rm d}r
                \ +\ r_p^2} {\pi R_s^2}
     \end{eqnarray}                                                         */
  res = ipv[ipn1]*ipv[ipn1] - 2.0*res;

  /* If the planet is going to be transparent with its maximum optical
     depth given by toomuch then:                                           */
  if(sg->transpplanet)
    res -= exp(-maxtau) * ipv[ipn-last] * ipv[ipn-last];

  /* Normalize by the stellar radius:                                       */
  res *= 1.0 / (srad*srad);

  free(hsum);
  free(hratio);
  free(hfactor);
  free(h);

  return res;
}


/* FUNCTION
   Calculate the modulation at a given wavenumber, considering the planet
   as an opaque disc of radius r = r(tau=toomuch), for no-limb darkening
   nor planet emission.

   Return: transit's modulation, else
           -1 if toomuch was not reached                               */
static inline PREC_RES
modulationm1(PREC_RES *tau,        /* Optical depth array              */
             long last,            /* Index where tau = toomuch        */
             double toomuch,       /* Maximum optical depth calculated */
             prop_samp *ip,        /* Impact parameter                 */
             struct geometry *sg){ /* Geometry struct                  */
  long i,         /* Auxiliary for-loop index         */
       ini;       /* Impact parameter index           */
  double ipv[2];  /* Impact parameter around toomuch  */
  double muchrad; /* Radius where tau reaches toomuch */
  double srad = sg->starrad * sg->starradfct; /* Stellar radius */

  /* toomuch was not reached: */
  if(tau[last] < toomuch)
    return -1;

  /* Get impact parameter before and after tau=toomuch: */
  ini = ++last-2;
  if(ini<0) ini = 0;
  for(i=ini; i<last; i++)
    ipv[i-ini] = ip->v[i]*ip->fct;

  /* Find the planet radius through a linear interpolation: */
  muchrad = interp_line(tau+ini, ipv, toomuch);

  /* Return the modulation: */
  return muchrad * muchrad / (srad*srad);
}


/* FUNCTION
   Wrapper function to calculate the modulation in/out-of-transit ratio.
   Return: modulation                                      */
static PREC_RES
modulationperwn(struct transit *tr,  /* Transit structure                   */
                PREC_RES *tau,       /* Optical depth                       */
                PREC_RES w,          /* Wavenumber value                    */
                long last,           /* Radius index where tau = toomuch    */
                double toomuch,      /* Maximum optical depth calculated    */
                prop_samp *ip){      /* Impact parameter array              */

  tr->modlevel = tr->ds.th->modlevel;
  struct geometry *sg  = tr->ds.sg;

  switch(tr->modlevel){
  case 1:
    return modulation1(tau,  last, toomuch, ip, sg);
    break;
  case -1:
    return modulationm1(tau, last, toomuch, ip, sg);
    break;
  default:
    tr_output(TOUT_ERROR, "slantpath:: modulationperwn:: Level %i of "
      "detail has not been implemented to compute modulation.\n",
      tr->modlevel);
    exit(EXIT_FAILURE);
    return 0;
  }
}



/* \fcnfh
   Print (to file or stdout) the modulation as function of wavelength */
void
printmod(struct transit *tr){
  FILE *outf = stdout;
  struct outputray *outray = tr->ds.out;
  int rn;

  /* Open file: */
  if(tr->f_outspec && tr->f_outspec[0] != '-')
    outf = fopen(tr->f_outspec, "w");

  tr_output(TOUT_INFO, "\nPrinting modulation spectrum in '%s'.\n",
            tr->f_outspec?tr->f_outspec:"standard output");

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

const ray_solution slantpath = {
  "transit",        /* Name of the solution                     */
  "slantpath.c",    /* Source code file name                    */
  0,                /* Request equispaced impact parameter      */
  &transittau,      /* Optical-depth calculation function       */
  &modulationperwn, /* Modulation calculation function          */
};
