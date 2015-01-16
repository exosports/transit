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
#include <extraext.h>


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

/* \fcnfh
   Computes optical depth for eclipse geometry for one ray, one wn, 
   at various incident angles on the planet surface, 
   between a certain layer in the atmosphere up to the top layer. 
   Returns: Optical depth divided by rad.fct:  \frac{tau}{units_{rad}}      */
/* DEF */
static PREC_RES
totaltau_eclipse(PREC_RES *rad,  /* Equispaced radius array: From current   */
                                 /* to outmost layer                        */
                 PREC_RES *ex,   /* Extinction[rad], all radii, one wn      */
                                 /* Often used as ex+ri than &ex[ri]        */
                 PREC_RES angle, /* Incident ray-path angle (degrees)       */
                 long nrad){     /* Number of radii between current and     */
                                 /* outmost layer                           */

  PREC_RES res;          /* Optical depth divided by units of radius        */
  PREC_RES x3[3], r3[3]; /* Interpolation variables                         */

  /* Conversion to radian:                                                  */
  PREC_RES angle_rad = angle * DEGREES;   

  /* Distance along the path:                                               */
  PREC_RES s[nrad];

  /* Returns 0 if at outmost layer. No distance travelled.                  */
  if(nrad == 1)  
    return 0.0;

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
    s[i] = s[i-1] + (rad[i] - rad[i-1])/cos(angle_rad);
  }

  /* Integrate extinction along the path:                                   */
  /* Use spline if GSL is available along with at least 3 points:           */
//#ifdef _USE_GSL
//  gsl_interp_accel *acc = gsl_interp_accel_alloc();
//  gsl_interp *spl = gsl_interp_alloc(gsl_interp_cspline, nrad);
//  gsl_interp_init(spl, s, ex, nrad);
//  res = gsl_interp_eval_integ(spl, s, ex, 0, s[nrad-1], acc);
//  gsl_interp_free(spl);
//  gsl_interp_accel_free(acc);
//#else
//#error non equispaced integration is not implemented without GSL
//#endif /* _USE_GSL */

  /* Safety mode: GSL is acting up sometimes                                */
  res = integ_trapz(s, ex, nrad);

  /* Optical depth divided by units of radius:                              */
  return res;
}


/* ################################### 
    GRID FOR INTENSITY CALCULATION
   ################################### */

/* \fcnfh
   Allocates intensity grid. 2D array of angles and wn to pack intensities
   Return: 0 on success  */

/* DEF */
int
intens_grid(struct transit *tr){ 
  prop_samp *wn = &tr->wns;                /* Wavenumber sample             */
  long int wnn = wn->n;                    /* Wavenumbers                   */

  /* Declaring indices:                                                     */
  long int i;                              /* For counting angles           */

  /* Reads number of angles from transithint structure                      */
  struct transithint *trh=tr->ds.th;
  long int an = trh->ann;                   

  /* Initialize tau structure:                                              */
  static struct optdepth tau;
  tr->ds.tau = &tau;

  /* Set maximum optical depth:                                             */
  if(tr->ds.th->toomuch > 0)
    tau.toomuch = trh->toomuch;

  /* Allocates array to store radii indices where tau reaches toomuch:      */
  tau.last = (long      *)calloc(wnn,        sizeof(long));
  /* Allocates 2D array [rad][wn]:                                          */
  tau.t    = (PREC_RES **)calloc(wnn,        sizeof(PREC_RES *));
  tau.t[0] = (PREC_RES  *)calloc(wnn*tr->rads.n, sizeof(PREC_RES  ));
  for(i=1; i < wnn; i++)
    tau.t[i] = tau.t[0] + i*tr->rads.n;

  /* Allocates intensity grid structure                                     */
  static struct grid intens; 
  memset(&intens, 0, sizeof(struct grid));  
  /* Connects the intensity gird structure to the transit structure         */
  tr->ds.intens = &intens;                 

  /* Allocates 2D array of intensities [angle][wn]:                         */
  intens.a    = (PREC_RES **)calloc(an   ,  sizeof(PREC_RES *));
  intens.a[0] = (PREC_RES  *)calloc(wnn*an, sizeof(PREC_RES  ));

  /* Connects intensity array with correct an and wn:                       */
  for(i = 1; i < an; i++)
    intens.a[i] = intens.a[0] + i * wnn;
  return 0;
}


/* ################################### 
    FILLS OUT OPTICAL DEPTH 2D ARRAY
   ################################### */

#define CIA_DOFLOAT  2
#define CIA_RADFIRST 1

/* \fcnfh
   Calculates optical depth as a function of radii and wavenumber.
   Return: 0 on success                                                     */

/* DEF */
int
tau_eclipse(struct transit *tr){           /* Transit structure             */

  struct optdepth *tau=tr->ds.tau;         /* Def optical depth structure   */
  //tr->ds.tau = &tau;                     /* Called from transit structure */
  prop_samp *rad = &tr->rads;              /* Radius sample                 */
  prop_samp *wn  = &tr->wns;               /* Wavenumber sample             */
  struct extinction *ex = tr->ds.ex;       /* Def extinction structure      */
  PREC_RES **e = ex->e;                    /* Extinction array [rad][wn]    */
  PREC_RES (*fcn)() = tr->ecl->tauEclipse; /* Function for totaltau_eclipse */
                                           /* See structure_tr.h line ~91   */

  long wi, ri; /* Indices for wavenumber, and radius                       */
  int rn;      /* Returns error code from computeextradiaus=>exstradis     */

  PREC_RES *r  = rad->v;         /* Values of radii                         */
  PREC_RES *tau_wn;              /* Array of optical depths [rad], one wn   */
  PREC_ATM *temp = tr->atm.t,    /* Temp array from atmospheric file        */
           tfct  = tr->atm.tfct; /* Temp units from atmospheric file        */

  /* Number of elements:                                                    */
  long int wnn = wn->n;        /* Wavenumbers                              */
  long int rnn = rad->n;       /* Radii                                    */
  double wfct  = wn->fct;      /* Wavenumber units factor to cgs           */

  struct transithint *trh = tr->ds.th;
  /* Pato Rojo explanation: 
  Consider:
  TRU_*BITS
  Those mark the bits that are carrying flags for each aspect. 
  i.e. TRU_EXTBITS is 1 for all the flag positions that are relevant
  for the extinction computation (only TRU_OUTTAU in that case); 
  TRU_ATMBITS is 1 for all the flag positions that are relevant for 
  the atmospheric computation; and so  on...therefore the line 64 
  passes all the TAU relevant flags from the user hint to 
  the main structure.  In particular, it only passes whether TRU_OUTTAU 
  was 1 or 0 as specified by the user. */
  transitacceptflag(tr->fl, tr->ds.th->fl, TRU_TAUBITS);

  ///* Set maximum optical depth:                                             */
  //if(tr->ds.th->toomuch > 0)   
  //  tau.toomuch = trh->toomuch; 

  ///* Allocates array to store radii indices where tau reaches toomuch:      */
  //tau.last = (long      *)calloc(wnn,        sizeof(long));
  ///* Allocates 2D array [rad][wn]:                                          */
  //tau.t    = (PREC_RES **)calloc(wnn,        sizeof(PREC_RES *));
  //tau.t[0] = (PREC_RES  *)calloc(wnn*rad->n, sizeof(PREC_RES  ));
  ///* Connects tau.t with correct wn and rad:                                */
  //for(wi=1; wi < wnn; wi++)
  //  tau.t[wi] = tau.t[0] + wi*rad->n;

  /* Set cloud structure:                                                   */
  /* Linear for rini to rfin, grey opacity from 0 to maxe.                  */
  static struct extcloud cl;
  cl.maxe = tr->ds.th->cl.maxe; /* Maximum opacity                          */
  cl.rini = tr->ds.th->cl.rini; /* Top layer radius                         */
  cl.rfin = tr->ds.th->cl.rfin; /* Radius of maxe                           */
  if(tr->ds.th->cl.rfct==0)
    cl.rfct = rad->fct;
  else
    cl.rfct = tr->ds.th->cl.rfct;
  tr->ds.cl = &cl;

  /* Tests if ex is calculated for particular radius, if true/False:        */
  _Bool *comp = ex->computed;

  /* Restors extinction from the savefile if given:                         */
  if(tr->save.ext)
    restfile_extinct(tr->save.ext, e, comp, rnn, wnn);

  /* Computes extinction at the outermost layer:                            */
  /* Pato Rojo explanation:                                                 */
  /* If outmost layer not computed r[lastr] will segfault otherwise.        */
  /* Last computed layer in that case is unnecessary.                       */
  /* Compute extinction when radius is below the last computed layer:       */
  if(!comp[rnn-1]){
    transitprint(1, verblevel, "Computing extinction in the outermost "
                               "layer.\n");
    if (tr->fp_opa != NULL)
      rn = interpolmolext(tr, rnn-1, ex->e);
    else
      if((rn=computemolext(tr, rnn-1, ex->e))!=0)
      transiterror(TERR_CRITICAL,
                   "computeextradius() returned error code %i.\n", rn);
  }

  /* Gets a copy of the radius units factor:                                */
  double rfct = rad->fct;

  /* Reads angle index from transit structure                               */
  long int angleIndex = tr->angleIndex;

  /* Reads the current angle for optical depth calculation                  */
  PREC_RES angle = trh->angles[angleIndex];

  /* Requests at least four radii to calculate a spline interpolation:      */
  if(rnn < 4)
    transiterror(TERR_SERIOUS, "tau(): At least four radius points (%d given) "
                 "are required. (three for spline and one for the analytical "
                 "part)\n", rnn);

  transitprint(1, verblevel, "Calculating optical depth at various radii for "
                             "angle %2.1lf degrees.\n", angle);

  PREC_RES er[rnn];        /* Array of extinction per radius                */

  int lastr = rnn-1;       /* Radius index of last computed extinction      */
                           /* from the outmost layer inwards                */

  /* Tenth of wavenumber sample size, used for progress printing:           */
  int wnextout = (long)(wnn/10.0); 

  /* Extinction from scattering and clouds:                                 */
  double e_s[rnn], e_c[rnn];
  /* Extinction from CIA:                                                   */
  PREC_CIA **e_cia = tr->ds.cia->e;
  /* Extinction from scattering:                                            */
  /* This is a hook for future handling of scattering.                      */
  struct extscat *sc = tr->ds.sc;

  /* Flow of this part:
     For each wavenumber of the selected scale, the optical depth along 
     the path with a given radius is computed. 
     The computation starts with a radius r[rnn], equal to the radius
     of the outmost layer and descends into the atmosphere until the 
     optical depth rises above a certain user-defined limit (--toomuch). 
     Each time this procedure reaches a radius whose extinction has not 
     been computed, it calls the extinction computing routine, 
     to compute the extinction over the whole wavenumber range.             */
  
  /* For each wavenumber:                                                   */
  for(wi=0; wi < wnn; wi++){

    /* Pointing to the declared tau array. All radii for one wn.            */
    tau_wn = tau->t[wi];

    /* Print output every 10\% that is ready:                               */
    if(wi > wnextout){
      transitprint(2, verblevel, "%i%%\n", (int)(100*(float)wi/wnn+0.5));
      wnextout += (long)(wnn/10.0);
    }

    /* Calculate extinction from scattering, clouds, and CIA.               */
    /* At each level just for the wn of interest:                           */
    computeextscat(e_s,  rnn, sc, rad->v, rad->fct, temp, tfct, wn->v[wi]*wfct);
    computeextcloud(e_c, rnn, &cl, rad, temp, tfct, wn->v[wi]*wfct);

    /* Calculate extinction for all radii if extinction exists:             */
    for(ri=0; ri < rnn; ri++)
      er[ri] = e[ri][wi] + e_cia[wi][ri] + e_s[ri] + e_c[ri];

    /* For each radii:                                                      */
    for(ri=rnn-1; ri > -1; ri--){
      transitprint(30, verblevel, "Radius r[%li] = %9.4g\n", ri, r[ri]);

      /* Compute extinction only if radius is smaller than the radius of
         the last calculated extinction:                                    */
      if(r[ri] < r[lastr]){
        transitprint(30, verblevel, "Last Tau (r=%9.4g, wn=%9.4g): %10.4g.\n",
                                   r[ri-1], wn->v[wi], tau_wn[ri-1]);
        if(ri)
          transitprint(30, verblevel, "Last Tau (r=%9.4g, wn=%9.4g): %10.4g.\n",
                                     r[ri-1], wn->v[wi], tau_wn[ri-1]);

        /* If extinction was not computed, compute it using computextradius. 
           The function calls extradius for one radius, but for all wn.
           Then, the code updates extinction for the particular wn. 
           Computation starts from the layer below outermost, since the 
           outmost layer is computed above.                                 */ 
                                          
        do{
          if(!comp[--lastr]){
            /* Compute a new extinction at given radius printing error if
               something happen:                                            */
            transitprint(2, verblevel, "Radius %i: %.9g cm ...\n",
                                       lastr+1, r[lastr]);
            if (tr->fp_opa != NULL)
              rn = interpolmolext(tr, lastr, ex->e);
            else
              if((rn=computemolext(tr, lastr, ex->e))!=0)
                transiterror(TERR_CRITICAL,
                           "computeextradius() return error code %i while "
                           "computing radius #%i: %g\n", rn, r[lastr]*rfct);
            /* Update the value of the extinction at the right place:       */
            er[lastr] = e[lastr][wi] + e_s[lastr] + e_c[lastr] +
                        e_cia[wi][lastr];
          }
        }while(r[ri]*rad->fct < r[lastr]*rfct);
      }

      //if (wi == 1612 || wi == 1611){
      //if (wi == 1611){
      //  transitprint(1,2, "ext(%3li) = tot:%.4e -- mol:%.4e -- cia:%.4e\n"
      //                    "    rad: %.2f   angle: %.2f  rnn-ri: %li\n",
      //                    ri, er[ri], e[ri][wi], e_cia[wi][ri],
      //                    *(r+ri), angle, rnn-ri);
      //}
      /* Compute tau_wn[ri] array until tau reaches toomuch
         first tau[0], starts from the outmost layer.                         
         Also compute tau.last (radii where tau reaches toomuch at each wn): */

      if( (tau_wn[rnn-ri-1] = rfct * fcn(r+ri, er+ri, angle, rnn-ri)) > tau->toomuch){
        tau->last[wi] = rnn-ri-1;

        if (ri < 3){
          transitprint(1, verblevel,
                       "WARNING: At wavenumber %g (cm-1), the critical TAU "
                       "value (%g) was exceeded with tau=%g at the radius "
                       "level %li (%g km), this should have "
                       "happened in a deeper layer (check IP sampling or ATM "
                       "file).\n", wn->v[wi], tau->toomuch, tau_wn[ri],
                       ri, r[ri]*rfct/1e5);
        }
        /* Exit tau-calculation loop if it reached toomuch:                 */
        break;
      }

      /* Sets tau of the outermost layer to zero:                           */
      tau_wn[0] = 0;

      transitDEBUG(22, verblevel, "Tau(lambda %li=%9.07g, r=%9.4g) : %g "
                                  "(toomuch: %g)\n", wi, wn->v[wi], r[ri],
                                  tau_wn[ri], tau->toomuch);
    }

    if(ri==rnn){
      transitprint(1, verblevel,
                   "WARNING: At wavenumber %g cm-1, the bottom of the "
                   "atmosphere was reached before obtaining the critical TAU "
                   "value of %g.\nMaximum TAU reached: %g.\n",
                   wn->v[wi], tau->toomuch, tau_wn[ri]);
      tau->last[wi] = ri-1;
    }
  }

  transitprint(1, verblevel, " Done.\nOptical depth calculated up to %g.\n",
               tr->ds.tau->toomuch);

  /* Print detailed output if requested:                                    */
  if(tr->ds.det->tau.n)
    detailout(&tr->wns, &tr->rads,  &tr->ds.det->tau, tau->t, 0);
  if(tr->ds.det->ext.n)
    detailout(&tr->wns, &tr->rads, &tr->ds.det->ext, e, CIA_RADFIRST);
  if(tr->ds.det->cia.n)
    detailout(&tr->wns, &tr->rads, &tr->ds.det->cia, (double **)e_cia,
              CIA_DOFLOAT);

  if(tr->save.ext)
    savefile_extinct(tr->save.ext, e, comp, rnn, wnn);

  /* Print lowest radius before optical depth gets too big:                 */
  if(tr->f_toomuch)
    printtoomuch(tr->f_toomuch, tr->ds.tau, &tr->wns, &tr->rads);

  /* Set progress indicator and output tau if requested, 
     otherwise return success:                                              */
  tr->pi |= TRPI_TAU;
  if(tr->fl & TRU_OUTTAU)
    printtau(tr);
  return 0;
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
               PREC_RES *tau,          /* Optical depth array tau.c==>tau   */
               PREC_RES w,             /* Current wavenumber value          */
               long last,              /* Index where tau = toomuch         */
               double toomuch,         /* Maximum optical depth calculated  */
               prop_samp *rad){        /* Radii array                       */
  /* FINDME: toomuch is not needed as a parameter        */

  /* General variables:                                                     */
  PREC_RES res;                  /* Result                                  */
  PREC_ATM *temp = tr->atm.t;    /* Temperatures                            */

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
  PREC_RES tauInteg[rnn],  /* Integrand function                            */
           tauIV[rnn];     /* Tau integration variable                      */

  /* Integrate for each of the planet's layer starting from the           
     outermost until the closest layer. 
     The order is opposite since tau starts from the top and 
     radius array starts from the bottom.                                   */

  /* Planck function (erg/s/sr/cm) for wavenumbers:
        B_\nu = 2 h {\bar\nu}^3 c^2 \frac{1}
                {\exp(\frac{h \bar \nu c}{k_B T})-1}                        */
  for(i=0; i <= last; i++){
    tauIV[i] = tau[i];
    B[i] =  (2.0 * H * w * w * w * wfct * wfct * wfct * LS * LS)
          / (exp(H * w * wfct * LS / (KB * temp[rnn-1-i])) - 1.0);
    tauInteg[i] = B[i] * exp(-tau[i]);
  }

  /* Added 0 at the end when tau reach toomuch, so the spline looks nice    */
  /* Add all other layers to be 0.                                          */
  for(; i<rnn; i++){
    tauInteg[i] = 0;
    /* Geometric progression is used to provide enough elements 
       for integral to work. It does not change the final outcome/result.   */
    tauIV[i] = tauIV[i-1] + 1;
   }

  /* Adding additional 0 layer, plus the last represent number of elements 
     is -1, so we need to add one more. 2 in total.                         */
  last += 2;

  /* If atmosphere is transparent, and at last level tau has not reached 
     tau.toomuch, last is set to max number of layers (rnn, instead of rnn-1
     because we added 2 on the previous step). The code requests never
     to go over it.                                                         */
  if(last > rnn)    
    last = rnn;

  /* Checks if we have enough radii to do spline, at least 3:               */
  if(last < 3)
    transiterror(TERR_CRITICAL, "Less than 3 items (%i given) for radial "
                                "integration.\n", last);

  /* Integrate along tau up to tau = toomuch:                               */
#ifdef _USE_GSL
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  gsl_interp *spl       = gsl_interp_alloc(gsl_interp_cspline, last);
  gsl_interp_init(spl, tauIV, tauInteg, last);
  res = gsl_interp_eval_integ(spl, tauIV, tauInteg,
                               tauIV[0], tauIV[last-1], acc);
  gsl_interp_free(spl);
  gsl_interp_accel_free (acc);
#else
# error computation of modulation() without GSL is not implemented
#endif

  /* GSL is stupid. I will use a trapezoidal rule for integration instead
     (when I want to run the code in safety mode):                          */
  //res = integ_trapz(tauIV, tauInteg, last);

  //if (fabs(w-2877.00) <= 0.5){
  //  transitprint(1, 2, "\nI(w=%.10g) = %.10g\n", w, res);
  //  for (i=0; i<last; i++)
  //    transitprint(1, 2, "  tau: %.10e   int:%.10e\n", tauIV[i], tauInteg[i]);
  //  double res2 = integ_trapz(tauIV, tauInteg, last-1);
  //  transitprint(1,2, "Trapezoidal integration: %.10e\n", res2);
  //}

  //if (res < 0){
  //if (fabs(w-1844.59) <= 0.005){
  //  double res2 = integ_trapz(tauIV, tauInteg, last-1);
  //  transitprint(1,2, "Trapezoidal integration: %.10e\n", res2);
  //  transitprint(1, 2, "\nI(w=%.10g) = %.10g\n", w, res);
  //  for (i=0; i<last; i++)
  //    transitprint(1, 2, "  tau: %.10e   int:%.10g\n", tauIV[i], tauInteg[i]);
  //}
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
  static struct outputray st_out;     /* Defines output structure           */
  tr->ds.out = &st_out;

  /* Initial variables:                                                     */
  long w;
  prop_samp *rad = &tr->rads;          /* Radius array pointer              */
  prop_samp *wn  = &tr->wns;           /* Wavenumber array pointer          */
  long int wnn   = wn->n;              /* Wavenumbers                       */
  eclipse_ray_solution *ecl = tr->ecl; /* Eclipse ray solution pointer      */

  /* Reads angle index from transit structure                               */
  long int angleIndex = tr->angleIndex;

  /* Intensity for all angles and all wn                                    */
  PREC_RES **intens_grid = tr->ds.intens->a;

  /* Intensity array for one angle all wn:                                  */
  PREC_RES *out = intens_grid[angleIndex];

  /* Reads the tau array from transit structure                             */
  struct optdepth *tau = tr->ds.tau;

  /* Integrate for each wavelength:                                         */
  transitprint(1, verblevel, "\nIntegrating over wavelength.\n");

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
    out[w] = ecl->eclIntenWn(tr, tau->t[w], wn->v[w], tau->last[w], 
                             tau->toomuch, rad);

    /* Prints to screen the progress status:                                 */
    if(w==nextw){
      nextw += wn->n/10;
      transitprint(2, verblevel, "%i%% ", (10*(int)(10*w/wn->n+0.9999999999)));
    }
  }

  transitprint(1, verblevel, "\nDone.\n");

  /* Sets progress indicator, and prints output:                             */
  tr->pi |= TRPI_MODULATION; /* FINDME: this is not a modulation calculation */
  printintens(tr);
  return 0;
}


/* #################
    CALCULATES FLUX
   ################# */

/* \fcnfh
   Calculates flux
   Formula: 
   Flux = pi * SUMM_i [I_i * (sin(theta_fin)^2 - sin(theta_in)^2)]
   I_i are calculated for each angle defined in the configuration file
   Returns: zero on success                                                */

/* DEF */
int
flux(struct transit *tr){  /* Transit structure                             */
  struct transithint *trh=tr->ds.th;
  static struct outputray st_out;
  tr->ds.out = &st_out;

  /* Get angles and number of angles from transithint:                      */
  PREC_RES *angles = trh->angles;     /* Angles                             */
  long int an = trh->ann;             /* Number of angles                   */

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

  /* Control code, to check areas:                                          */
  //for(i=0; i<an+1; i++)
  //  printf("Area grid [%ld] = %lf\n", i, area_grid[i]/DEGREES);  

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
  /* HACK: */
  if (verblevel > 0)
    printflux(tr);
  return 0;
}


/* ############################
    OUTPUT FILE FOR INTENSITIES
   ############################ */

/* \fcnfh
   Print (to file or stdout) the emergent intensities as function of wavelength
   for each angle)                                                          */

void
printintens(struct transit *tr)
{

  /* Intensity for all angles and all wn                                    */
  PREC_RES **intens_grid = tr->ds.intens->a;   

  /* Takes number of angles from the transithint structure                  */
  struct transithint *trh=tr->ds.th;
  long int an = trh->ann;                   /* Number of angles            */

  prop_samp *wn = &tr->wns;                  /* Wavenumber sample           */
  long int wnn = wn->n;                     /* Wavenumbers                 */

  /* Declaring indexes:                                                     */
  long int i;                               /* For counting angles         */
  long int w;                               /* For counting wn             */

  /* Reads angles from transithint structure                                */
  PREC_RES *angles = trh->angles;      

  /* Declares stdout                                                        */
  FILE *outf=stdout;

  /* Adds string to the output files to differentiate between outputs        */
  char our_fileName[512];
  strncpy(our_fileName, tr->f_out, 512);
  strcat(our_fileName, ".-Intens");

  /* Opens a file:                                                          */
  if(tr->f_out&&tr->f_out[0]!='-')
    outf=fopen(our_fileName,"w");

  transitprint(1, verblevel, "\nPrinting intensity in '%s'\n",
               tr->f_out ? tr->f_out:"standard output");

  /* Prints:                                                                */
  char wlu[20],wnu[20];                       /* Wavelength units name     */
  long nsd=(long)(1e6);                      /* Wavenumber units name      */

  /* Get wavenumber units name:                                             */
  if((long)(nsd*tr->wns.fct)==nsd) strcpy(wnu,"cm");
  else if((long)(nsd*1e-1*tr->wns.fct)==nsd) strcpy(wnu,"mm");
  else if((long)(nsd*1e-4*tr->wns.fct)==nsd) strcpy(wnu,"um");
  else if((long)(nsd*1e-7*tr->wns.fct)==nsd) strcpy(wnu,"nm");
  else if((long)(nsd*1e-8*tr->wns.fct)==nsd) strcpy(wnu,"a");
  else sprintf(wnu,"%6.1g cm",1/tr->wns.fct);

  /* Get wavenumber units name:                                              */
  if((long)(nsd*tr->wavs.fct)==nsd) strcpy(wlu,"cm");
  else if((long)(nsd*1e1*tr->wavs.fct)==nsd) strcpy(wlu,"mm");
  else if((long)(nsd*1e4*tr->wavs.fct)==nsd) strcpy(wlu,"um");
  else if((long)(nsd*1e7*tr->wavs.fct)==nsd) strcpy(wlu,"nm");
  else if((long)(nsd*1e8*tr->wavs.fct)==nsd)  strcpy(wlu,"a");
  else sprintf(wlu,"%8.1g cm",tr->wavs.fct);

  /* Prints the header:                                                      */
  fprintf(outf, "#wvl [um]%*s", 6, " ");
  for(i=0; i < an; i++)
      fprintf(outf, "I[%4.1lf deg]%*s", angles[i], 7, " ");
  fprintf(outf, "[erg/s/cm/sr]\n");

  /* Fills out each column with the correct output intensity                 */
  for(w=0; w<wnn; w++){
    fprintf(outf, "%-15.10g", 1e4/(tr->wns.v[w]/tr->wns.fct));
    for(i = 0; i < an; i++)
      fprintf(outf, "%-18.9g", intens_grid[i][w]);
    fprintf(outf,"\n");
  }

  /* Closes the file:                                                        */
  fclose(outf);
  return;
}


/* ######################
    OUTPUT FILE FOR FLUX
   ###################### */

/* \fcnfh
   Print (to file or stdout) the emergent intensity as function of wavenumber 
   (and wavelength)                                                         */

void
printflux(struct transit *tr){
  FILE *outf=stdout;
  /* The flux per wavenumber array:                                         */
  PREC_RES *Flux = tr->ds.out->o; 
  int rn;

  /* Adds string to the output files to differentiate between outputs:      */
  char our_fileName[512];
  strncpy(our_fileName, tr->f_out, 512);
  strcat(our_fileName, ".-Flux");

  /* Opens a file:                                                          */
  if(tr->f_out && tr->f_out[0] != '-')
    outf=fopen(our_fileName, "w");

  transitprint(1, verblevel, "\nPrinting flux in '%s'\n",
               tr->f_out ? our_fileName:"standard output");

  /* Prints:                                                                */
  char wlu[20], wnu[20];         /* Wavelength units name                   */
  long nsd=(long)(1e6);        /* Wavenumber units name                     */

  /* Get wavenumber units name:                                             */
  if((long)(nsd*tr->wns.fct)==nsd) strcpy(wnu,"cm");
  else if((long)(nsd*1e-1*tr->wns.fct)==nsd) strcpy(wnu,"mm");
  else if((long)(nsd*1e-4*tr->wns.fct)==nsd) strcpy(wnu,"um");
  else if((long)(nsd*1e-7*tr->wns.fct)==nsd) strcpy(wnu,"nm");
  else if((long)(nsd*1e-8*tr->wns.fct)==nsd) strcpy(wnu,"a");
  else sprintf(wnu, "%6.1g cm", 1/tr->wns.fct);

  /* Get wavenumber units name:                                             */
  if((long)(nsd*tr->wavs.fct)==nsd) strcpy(wlu,"cm");
  else if((long)(nsd*1e1*tr->wavs.fct)==nsd) strcpy(wlu,"mm");
  else if((long)(nsd*1e4*tr->wavs.fct)==nsd) strcpy(wlu,"um");
  else if((long)(nsd*1e7*tr->wavs.fct)==nsd) strcpy(wlu,"nm");
  else if((long)(nsd*1e8*tr->wavs.fct)==nsd)  strcpy(wlu,"a");
  else sprintf(wlu,"%8.1g cm",tr->wavs.fct);

  /* Prints the header:                                                     */
  fprintf(outf, "#wvl [um]%*sFlux [erg/s/cm]\n", 6, " ");

  /* Prints wavelength and flux:                                            */
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


/* ############################
    CALLS ECLIPSE RAY SOLUTION
   ############################ */

  /* Executes eclipse module:                                               */
const eclipse_ray_solution eclipsepath = {
       "Eclipse Path",   /* Name of the solution                            */
       "eclipse.c",      /* Source code file name                           */
       &totaltau_eclipse,/* Per per wavenumber value computation            */
       &eclipse_intens,  /* Per wavenumber value computation                */
       };


