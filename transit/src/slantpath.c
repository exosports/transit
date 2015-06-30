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

  /* Closest approach radius:                                               */
  PREC_RES r0 = b/refr;

  /* Get the index rs, of the sampled radius immediately below or equal
     to r0 (i.e. rad[rs] <= r0 < rad[rs+1]):                                */
  rs = binsearch(rad, 0, nrad-1, r0);
  if ((rs == -5) || (rs == -2))
    return 0;  /* If it is the outmost layer                                */
  /* If some other error occurred:                                          */
  else if(rs < 0){
    transiterror(TERR_CRITICAL, "Closest approach (%g) is larger than the "
                 "top layer of the atmosphere (%g).\n", r0, rad[nrad-1]);
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
  //transitprint(100, verblevel, "Delta-s = [");
  for(i=1; i<nrad; i++){
    s[i] = sqrt(rad[i]*rad[i] - r0*r0);
    //transitprint(100, verblevel, "%.6f, ", s[i]);
  }
  //transitprint(100, verblevel, "]\n");

  /* Integrate Extinction along ray path:                                   */
  /* Use spline if GSL is available along with at least 3 points:           */


  double *hsum;
  double *hratio;
  double *hfactor;
  double *h;

  hsum    = calloc(nrad/2, sizeof(double));
  hratio  = calloc(nrad/2, sizeof(double));
  hfactor = calloc(nrad/2, sizeof(double));
  h       = calloc(nrad-1, sizeof(double));
  
  makeh(s, h, nrad);

  geth(h, hsum, hratio, hfactor, nrad);

  res = simps(ex, h, hsum, hratio, hfactor, nrad);

  /*
#ifdef _USE_GSL
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_interp       *spl = gsl_interp_alloc(gsl_interp_cspline, nrad);
  gsl_interp_init(spl, s, ex, nrad);
  res = gsl_interp_eval_integ(spl, s, ex, 0, s[nrad-1], acc);
  gsl_interp_free(spl);
  gsl_interp_accel_free(acc);
#else
#error non equispaced integration is not implemented without GSL
#endif */ /* _USE_GSL                                                          */

  /* Reset original values of arrays:                                       */
  *ex  = tmpex;
  *rad = tmprad;

  /* Return:                                                                */
  return 2*(res);
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

  transiterror(TERR_CRITICAL|TERR_ALLOWCONT, "This routine has not been "
               "tested yet.  You have been warned.\n");

  /* Look for closest approach radius:                                      */
  while(1){
    r0 = b/lineinterp(r0a, rad, refr, nrad);
    if(r0==r0a)
      break;
    if(i++ > maxiterations)
      transiterror(TERR_CRITICAL, "Maximum iterations (%i) reached while "
                   "looking for r0. Convergence not reached (%.6g!=%.6g).\n",
                   maxiterations, r0, r0a);
    r0a = r0;
  }

  /* Get bin value 'rs' such that r0 is between rad[rs-1] inclusive
     and rad[rs] exclusive:                                                 */
  /* If we are looking at the outmost layer, then return:                   */
  if((rs=binsearch(rad, 0, nrad-1, r0))==-5)
    return 0;

  /* If some other error occurred:                                          */
  else if(rs<0)
    transiterror(TERR_CRITICAL,
                 "Closest approach value(%g) is outside sampled radius "
                 "range(%g - %g).\n", r0, rad[0], rad[nrad-1]);
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

  /* Integrate: Use spline if GSL is available along with at least 3 points: */
#ifdef _USE_GSL
  if(nrad-rs > 2){
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *spl = gsl_spline_alloc(gsl_interp_cspline, nrad-rs);
    gsl_spline_init(spl, rad+rs, dt+rs, nrad-rs);
    res += gsl_spline_eval_integ(spl, rad[rs], rad[nrad-1], acc);
    gsl_spline_free(spl);
    gsl_interp_accel_free(acc);
  }
  else /* Else, use trapezoidal integration when there are only two points: */
#endif /* _USE_GSL */
  if(nrad-rs > 1)
    res += integ_trasim(rad[1]-rad[0], dt+rs, nrad-rs);

  return 2*(res);
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
    transiterror(TERR_CRITICAL,
                 "slantpath:: totaltau:: Level %i of detail has not been "
                 "implemented to compute optical depth.\n", tr->taulevel);
    return 0;
  }
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
  if(last<3)
    transiterror(TERR_CRITICAL, "Condition failed, less than 3 items "
                                "(only %i) for radial integration.\n", last);

  /* Integrate along radius:                                                */
#ifdef _USE_GSL
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_interp *spl       = gsl_interp_alloc(gsl_interp_cspline, last);
  gsl_interp_init(spl, ipv+ipn-last, rinteg+ipn-last, last);
  res = gsl_interp_eval_integ(spl, ipv+ipn-last, rinteg+ipn-last,
                               ipv[ipn-last], ipv[ipn1], acc);
  gsl_interp_free(spl);
  gsl_interp_accel_free (acc);
#else
# error computation of modulation() without GSL is not implemented
#endif

  /* TD: Add real unblocked area of the star, considering geometry          */
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
    transiterror(TERR_CRITICAL, "slantpath:: modulationperwn:: Level %i of "
                 "detail has not been implemented to compute modulation.\n",
                 tr->modlevel);
    return 0;
  }
}

const ray_solution slantpath = {
  "transit",        /* Name of the solution                     */
  "slantpath.c",    /* Source code file name                    */
  0,                /* Request equispaced impact parameter      */
  &transittau,      /* Optical-depth calculation function       */
  &modulationperwn, /* Modulation calculation function          */
};
