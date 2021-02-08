// Copyright (C) 2015-2016 University of Central Florida. All rights reserved.
// Transit is under an open-source, reproducible-research license (see LICENSE).

#include <transit.h>

/* \fcnfh DEF
   Create a sampling array. Take values from hint or else from a
   reference sampling.

   Get units factor.  Get inital and final values.
   Check that there is a non-zero interval.
   Check that only the spacing or the number of points have been defined.
   If numpoints was provided, use the given array of values and ommit
   oversampling.  If spacing was provided, calculate numpoints, oversample,
   and fill in values.

   Return: 1 for modified initial value
           2 for modified final   value
           0 if nothing was changed but there is a sampled array
          -1 if hinted initial is bigger than maximum allowed
          -2 if hinted final is smaller than minimum allowed
          -3 if accepted initial value is greater or equal to final one
          -4 if both spacing and number of elements were hinted
          -5 if none or both of spacing or number of elements were in
             the referenced
          -6 Not valid oversampling was given by ref when requested      */
int
makesample1(prop_samp *samp,       /* transit sampling    */
            prop_samp *ref,        /* reference sampling  */
            const long fl){        /* Sampling flag       */

  int res=0;              /* Return output             */
  PREC_RES *v;            /* The sampling values       */
  double osd,             /* Oversampled delta         */
         si;              /* Sample initial value      */
  _Bool dhint=ref->d!=0;  /* True if hint.d is defined */
  /* Acceptable ratio exceeding final value without truncating the last bin: */
  double okfinalexcess = 1e-8;

  /* Get units factor: */
  samp->fct = ref->fct;

  /* Get initial value: */
  si = samp->i = ref->i;

  /* Get final value: */
  samp->f = ref->f;
  /* Check non-zero interval: */
  if (samp->f < samp->i){
    tr_output(TOUT_ERROR,
      "Hinted final value for %s sampling (%g) is smaller than "
      "hinted initial value %.8g.\n", TRH_NAME(fl), samp->f, samp->i);
    return -3;
  }

  /* Progress report: print flag, spacing, and number of points from hint: */
  tr_output(TOUT_DEBUG, "Flags: 0x%lx    hint.d: %g   hint.n: %li\n",
    fl, ref->d, ref->n);

  if (!dhint){
    /* If none of the ref exists then throw error: */
    tr_output(TOUT_ERROR,
      "Spacing (%g) was not hinted in %s sampling.\n", ref->d, TRH_NAME(fl));
    return -5;
  }
  /* If spacing exists then trust it: */
  else if(ref->d!=0){
      samp->d = ref->d;
  }
  /* Else: */
  else{
    tr_output(TOUT_ERROR,
      "Invalid spacing (%g) in %s sampling.\n", samp->d, TRH_NAME(fl));
    exit(EXIT_FAILURE);
  }

  /* Make sampling based on spacing:       */
  if(samp->d<0) okfinalexcess = -okfinalexcess;
  /* Define number of points:              */
  samp->n = ((1.0+okfinalexcess)*samp->f - samp->i)/samp->d+1;
  if(samp->n<0)
    samp->n = -samp->n;

  /* Check for hinted oversampling:        */
  if(ref->o <= 0){
    tr_output(TOUT_ERROR,
      "Invalid hinted oversampling for %s sampling.\n", TRH_NAME(fl));
    return -6;
  }
  samp->o = ref->o;

  /* Oversampled number of points:      */
  samp->n = (samp->n-1)*samp->o+1;
  /* Oversampled delta:                 */
  osd = samp->d/(double)samp->o;

  /* Allocate and fill sampling values: */
  v = samp->v = (PREC_RES *)calloc(samp->n, sizeof(PREC_RES));
  /* Fill-in values:                    */
  PREC_NREC n = samp->n;
  *v = si;
  v += --n;
  while(n)
    *v-- = si + n--*osd;
  /* FINDME: Why so complicated? This is far easier:
  for (i=0; i<samp->n; i++)
    *v+i = samp->i + i*osd               */

  /* Check the final point: */
  if(si!=0 && samp->v[samp->n-1]!=samp->f && verblevel>2)
    /* FINDME: Consider removig the verblevel condition */
    tr_output(TOUT_WARN,
      "Final sampled value (%g) of the %li points doesn't coincide "
      "exactly with required value (%g). %s sampling with "
      "pre-oversampling spacing of %g.\n", samp->v[samp->n-1],
      samp->n, samp->f, TRH_NAME(fl), samp->d);

  /* Return the flags of accepted values: */
  return res;
}


/* \fcnfh  DEF
   Create a sampling array. Take values from hint or else from a
   reference sampling.

   Get units factor.  Get inital and final values.
   Check that there is a non-zero interval.
   Check that only the spacing or the number of points have been defined.
   If numpoints was provided, use the given array of values and ommit
   oversampling.  If spacing was provided, calculate numpoints, oversample,
   and fill in values.

   Return: 1 for modified initial value
           2 for modified final   value
           0 if nothing was changed but there is a sampled array
          -1 if hinted initial is bigger than maximum allowed
          -2 if hinted final is smaller than minimum allowed
          -3 if accepted initial value is greater or equal to final one
          -4 if both spacing and number of elements were hinted
          -5 if none or both of spacing or number of elements were in
             the referenced
          -6 Not valid oversampling was given by ref when requested      */
int
makesample(prop_samp *samp,  /* transit sampling                            */
           prop_samp *hint,  /* hint    sampling                            */
           prop_samp *ref,   /* reference sampling                          */
           const long fl){   /* Sampling flag                               */

  int res=0;              /* Return output                                  */
  PREC_RES *v;            /* The sampling values                            */
  double osd,             /* Oversampled delta                              */
         si;              /* Sample initial value                           */
  _Bool dhint=hint->d!=0; /* True if hint.d is defined                      */
  /* Acceptable ratio exceeding final value without truncating the last bin: */
  double okfinalexcess = 1e-8;

  /* Get units factor:                                                      */
  if(hint->fct<=0)
    samp->fct = ref->fct;
  else
    samp->fct = hint->fct;

  /* Set initial value:                                                     */
  /* If hint unset, use reference:                                          */
  if(hint->i <= 0){
    samp->i = ref->i;
    tr_output(TOUT_INFO,
      "Using ref sampling %g [cgs] for initial value of %s.\n",
      samp->i*samp->fct, TRH_NAME(fl));
    res |= 0x1;  /* Turn on modification flag for return output             */
  }
  else
    samp->i = hint->i;

  /* Set final value (if hint unset, use reference):                        */
  if (hint->f <= 0){
    samp->f = ref->f;
    tr_output(TOUT_INFO,
      "Using ref sampling %g [cgs] for final value of %s.\n",
      samp->f*samp->fct, TRH_NAME(fl));
    res |= 0x2; /* Turn on modification flag for return output              */
  } else
    samp->f = hint->f;

  /* Print flag, spacing, and number of points from hint:                   */
  tr_output(TOUT_DEBUG, "Flags: 0x%lx    hint.d: %g   hint.n: %li\n",
    fl, hint->d, hint->n);

  /* If dhint has not been hinted then use ref's:                           */
  if(!dhint){
    /* If none of the ref exists then throw error:                          */
    if((ref->d==0 && ref->n<=0)){
      tr_output(TOUT_ERROR,
        "Spacing (%g) and number of elements (%i) were either both "
        "or none in the reference for %s sampling. And yes, none "
        "were hinted.\n", ref->d, ref->n, TRH_NAME(fl));
      return -5;
    }
    /* If spacing exists then trust it:                                     */
    if(ref->d!=0){
      samp->d = ref->d;
    }
    /* Else, use set array: */
    else{
      /* If initial or final value were modified, then print out a warning: */
      if(res){
        tr_output(TOUT_WARN,
          "Array of length %i was given as reference for %s "
          "sampling, but the initial (%g -> %g) or final "
          "(%g -> %g) values MIGHT have been modified.\n",
          ref->n, TRH_NAME(fl), ref->i, samp->i, ref->f, samp->f);
      }
      samp->n = ref->n;
      samp->d = 0;
      samp->v = (PREC_RES *)calloc(samp->n, sizeof(PREC_RES));
      memcpy(samp->v, ref->v, samp->n*sizeof(PREC_RES));
      if(ref->o != 0)
        tr_output(TOUT_WARN,
          "Fixed sampling array of length %i was referenced. "
          "But also oversampling was given (%li), ignoring it "
          "in %s sampling.\n", samp->n, ref->o, TRH_NAME(fl));
      samp->o = 0;
      /* Return any possible modification: */
      return res;
    }
  }
  /* Else if spacing was hinted, then it has to be positive at this point: */
  else if(dhint){
    transitASSERT(hint->d <= 0, "Error: Logic test 1 failed in %s's "
                                "makesample()\n", TRH_NAME(fl));
    samp->d = hint->d;
  }
  else{
    /* n can't be hinted: */
    tr_output(TOUT_ERROR, "Invalid sampling inputs.\n");
    exit(EXIT_FAILURE);
  }

  /* Check non-zero interval: */
  if((samp->f <= samp->i && samp->d > 0) ||
     (samp->f >= samp->i && samp->d < 0)){
    tr_output(TOUT_ERROR,
      "Initial accepted sampling value (%g) is greater or equal "
      "than final accepted sample value (%g). %s was being "
      "hinted.\n", samp->i, samp->f, TRH_NAME(fl));
    return -3;
  }

  /* Make sampling based on spacing:       */
  if(samp->d<0) okfinalexcess = -okfinalexcess;
  /* Define number of points:              */
  samp->n = ((1.0+okfinalexcess)*samp->f - samp->i)/samp->d + 1;
  if(samp->n<0)  /* FINDME: Explain how can this happen */
    samp->n = -samp->n;

  /* Check for hinted oversampling:        */
  if(hint->o<=0){
    /* If not, check for ref oversampling: */
    if(ref->o<=0){  /* If not, throw error */
      tr_output(TOUT_ERROR,
        "Not valid oversampling in the reference for "
        "%s sampling.\n", TRH_NAME(fl));
      return -6;
    }
    samp->o = ref->o;
  }
  else
    samp->o = hint->o;

  /* Oversampled number of points:      */
  samp->n = (samp->n - 1)*samp->o + 1;
  /* Oversampled delta:                 */
  osd = samp->d/(double)samp->o;

  /* Allocate and fill sampling values: */
  si = samp->i;  /* FINDME: delete si and keep using samp->i */
  v = samp->v = (PREC_RES *)calloc(samp->n, sizeof(PREC_RES));
  /* Fill-in values:                    */
  PREC_NREC n = samp->n;
  *v = si;
  v += --n;
  while(n)
    *v-- = si + n--*osd;
  /* FINDME: Why so complicated? This is far easier:
  for (i=0; i<samp->n; i++)
    *v+i = samp->i + i*osd               */

  /* Check the final point: */
  if(si != 0  &&  samp->v[samp->n-1] != samp->f  &&  verblevel > 2)
    /* FINDME: Consider removig the verblevel condition */
    tr_output(TOUT_WARN,
      "Final sampled value (%g) of the %li points doesn't coincide "
      "exactly with required value (%g). %s sampling with "
      "pre-oversampling spacing of %g.\n", samp->v[samp->n-1],
      samp->n, samp->f, TRH_NAME(fl), samp->d);

  /* Return the flags of accepted values: */
  return res;
}


/* \fcnfh  DEF
   Call makesample with the appropiate parameters and set the flags
   for a wavenumber sampling

   Return: makesample() output                                              */
int
makewnsample(struct transit *tr){
  struct transithint *th = tr->ds.th;
  prop_samp *hsamp  = &th->wns;          /* Hinted wavenumber sampling      */
  prop_samp *wlsamp = &th->wavs;         /* Hinted wavelength sampling      */
  prop_samp rsamp;                       /* Reference wavenumber sampling   */
  memset(&rsamp, 0, sizeof(prop_samp));
  int res;                               /* Result status                   */

  /* Get initial wavenumber value:                                          */
  if (hsamp->i > 0){
    if (hsamp->fct <= 0) {
      tr_output(TOUT_ERROR, "User specified wavenumber factor is "
        "negative (%g).\n", hsamp->fct);
      exit(EXIT_FAILURE);
    }
    rsamp.i = hsamp->i*hsamp->fct;
    tr_output(TOUT_RESULT, "wave i1: %.3f = %.2f * %.2f\n", rsamp.i,
      hsamp->i, hsamp->fct);
  }
  else if (wlsamp->f > 0){
    if (wlsamp->fct <= 0) {
      tr_output(TOUT_ERROR, "User specified wavelength factor is "
        "negative (%g).\n", wlsamp->fct);
      exit(EXIT_FAILURE);
    }
    rsamp.i = 1.0/(wlsamp->f*wlsamp->fct);
  }
  else {
    tr_output(TOUT_ERROR, "Initial wavenumber (nor final wavelength) "
      "were correctly provided by the user.\n");
    exit(EXIT_FAILURE);
  }

  /* Get final wavenumber value:                                            */
  if (hsamp->f > 0){
    if (hsamp->fct < 0) {
      tr_output(TOUT_ERROR, "User specified wavenumber factor is "
        "negative (%g).\n", hsamp->fct);
      exit(EXIT_FAILURE);
    }
    rsamp.f = hsamp->f*hsamp->fct;
  }
  else if (wlsamp->i > 0){
    if (wlsamp->fct < 0) {
      tr_output(TOUT_ERROR, "User specified wavelength factor is "
        "negative (%g).\n", wlsamp->fct);
      exit(EXIT_FAILURE);
    }
    rsamp.f = 1.0/(wlsamp->i*wlsamp->fct);
  }
  else {
    tr_output(TOUT_ERROR, "Final wavenumber (nor initial wavelength) "
      "were correctly provided by the user.\n");
    exit(EXIT_FAILURE);
  }

  /* Set up reference wavenumber sampling:                                  */
  rsamp.o = hsamp->o;

  /* Reference's wavenumber units are cm-1:                                 */
  rsamp.fct = 1;

  /* Don't set the number of elements:                                      */
  rsamp.n = 0;

  /* Set spacing such that the wavenumber grid has the same
     number of points as the wavelength grid:                */
  if (hsamp->d <= 0) {
    tr_output(TOUT_ERROR,
      "Incorrect wavenumber spacing (%g), it must be positive.\n", hsamp->d);
    exit(EXIT_FAILURE);
  }
  rsamp.d = hsamp->d;

  /* Make the oversampled wavenumber sampling:                              */
  res = makesample1(&tr->owns, &rsamp, TRH_WN);
  /* Make the wavenumber sampling:                                          */
  rsamp.o = 1;
  res = makesample1(&tr->wns,  &rsamp, TRH_WN);
  /* Get the exact divisors of the oversampling factor:                     */
  tr->odivs = divisors(tr->owns.o, &tr->ndivs);
  tr_output(TOUT_DEBUG, "There are %i divisors of the oversampling "
    "factor (%i):\n", tr->ndivs, tr->owns.o);
  for (int i=0; i<tr->ndivs; i++)
    tr_output(TOUT_DEBUG, "%5i", tr->odivs[i]);
  tr_output(TOUT_DEBUG, "\n");

  /* Set progress indicator if sampling was successful and return status:   */
  if(res >= 0)
    tr->pi |= TRPI_MAKEWN;
  return res;
}


/* \fcnfh  DEF
 Calls makesample with the appropiate parameters and set the flags

 @returns makesample() output
          1 Only one point value was requested                              */
int
makeradsample(struct transit *tr){
  int res,    /* Return output                                              */
      i, j,   /* Auxiliary for-loop indices                                 */
      iso1db; /* First isotope's index in a database                        */
  struct lineinfo *li=tr->ds.li;    /* Lineinfo struct                      */
  struct atm_data *atms=tr->ds.at;  /* Atmosphere data struct               */
  struct isotopes  *iso=tr->ds.iso; /* Isotopes struct                      */
  struct molecules *mol=tr->ds.mol; /* Molecules struct                     */

  int nrad,            /* Number of radii to sample                         */
      niso=iso->n_i,   /* Number of isotopes                                */
      ndb =iso->n_db,  /* Number of data bases                              */
      nmol=mol->nmol;  /* Number of molecules                               */

  long flag;                      /* Flag for interpolation type            */
  prop_isov *isovs;               /* lineinfo's isotopes information        */
  prop_samp *rsamp = &atms->rads; /* Atmosphere's radii sampling            */
  prop_samp *rad   = &tr->rads;   /* Output radius sampling                 */
  prop_atm  *atmt  = &tr->atm;    /* Array to store p, t, and mm sampling   */
  prop_mol  *molec = mol->molec;  /* Molecular variable information         */


  /* Check that getatm() and readinfo_tli() have been already called:       */
  transitcheckcalled(tr->pi, "makeradsample", 2, "getatm",       TRPI_GETATM,
                                                 "readinfo_tli", TRPI_READINFO);
  /* Exception for re-runs:                                                 */
  if (tr->pi & TRPI_MAKERAD){
    /* Free memory before re-allocating them:                               */
    free_atm(atmt);
    for (i=0; i<nmol; i++)
      free_mol(mol->molec+i);
    for (i=0; i<niso; i++)
      free_isov(iso->isov+i);
    freemem_samp(rad);
    tr->pi &= ~(TRPI_MAKERAD);
  }
  /* Check that variables are not NULL:                                     */
  transitASSERT(atms->rads.n<1 || !nmol,
  //transitASSERT(atms->rads.n<1 || !ndb || !niso || !nmol,
                "makeradsample():: called but essential variables are "
                "missing (%d, %d).\n", atms->rads.n, nmol);
                //"missing (%d, %d, %d, %d).\n", atms->rads.n, ndb, niso, nmol);

  /* Set interpolation function flag:                                       */
  flag = tr->interpflag;
  tr_output(TOUT_DEBUG, "transit interpolation flag: %li.\n", flag);

  /* We need to set-up limit so that the hinted values are compatible
     with the atmosphere                                                    */

  /* If there is only one atmospheric point, don't do makesample:           */
  if(rsamp->n==1){
    rad->n    = 1;
    rad->i    = rsamp->i;
    rad->f    = rsamp->f;
    rad->fct  = rsamp->fct;
    rad->d    = 0;
    rad->v    = (PREC_RES *)calloc(1, sizeof(PREC_RES));
    rad->v[0] = rsamp->v[0];
    res       = 0;   /* makesample()-like output                            */
    /* FINDME: warn that hinted values are going to be useless              */
  }
  /* Keep atmospheric-file sampling:                                        */
  else if (tr->ds.th->rads.d == -1){
    rad->n    = rsamp->n;
    rad->i    = rsamp->i;
    rad->f    = rsamp->f;
    rad->fct  = rsamp->fct;
    rad->d    = 0;
    rad->v    = (PREC_RES *)calloc(rad->n, sizeof(PREC_RES));
    for (i=0; i < rad->n; i++)
      rad->v[i] = rsamp->v[i];
    res       = 0;
  }
  /* Resample to equidistant radius array:                                  */
  else{
    res = makesample(rad, &tr->ds.th->rads, rsamp, TRH_RAD);
  }
  nrad = rad->n;

  /* Allocate arrays that will receive the interpolated data:               */
  for(i=0; i<nmol; i++){
    molec[i].d = (PREC_ATM *)calloc(nrad, sizeof(PREC_ATM));
    molec[i].q = (PREC_ATM *)calloc(nrad, sizeof(PREC_ATM));
    molec[i].n = nrad;
  }
  for (i=0; i<niso; i++){
    iso->isov[i].z = (PREC_ZREC *)calloc(nrad, sizeof(PREC_ZREC));
    iso->isov[i].n = nrad;
  }

  atmt->tfct = atms->atm.tfct;
  atmt->pfct = atms->atm.pfct;
  atmt->t  = (PREC_ATM *)calloc(nrad, sizeof(PREC_ATM));
  atmt->p  = (PREC_ATM *)calloc(nrad, sizeof(PREC_ATM));
  atmt->mm = (double   *)calloc(nrad, sizeof(double));

  /* Interpolate the atm. temperature, pressure, and mean molecular mass:   */
  splinterp(rsamp->n, rsamp->v, atms->atm.t, nrad, rad->v, atmt->t);
  splinterp(rsamp->n, rsamp->v, atms->atm.p, nrad, rad->v, atmt->p);
  splinterp(rsamp->n, rsamp->v, atms->mm,    nrad, rad->v, atmt->mm);

  /* Temperature boundary check:                                            */
  for (i=0; i<nrad; i++){
    if (atmt->t[i] < li->tmin) {
      tr_output(TOUT_ERROR, "The layer %d in the atmospheric model has "
        "a lower temperature (%.1f K) than the lowest allowed "
        "TLI temperature (%.1f K).\n", i, atmt->t[i], li->tmin);
      exit(EXIT_FAILURE);
    }
    if (atmt->t[i] > li->tmax) {
      tr_output(TOUT_ERROR, "The layer %d in the atmospheric model has "
        "a higher temperature (%.1f K) than the highest allowed "
        "TLI temperature (%.1f K).\n", i, atmt->t[i], li->tmax);
      exit(EXIT_FAILURE);
    }
  }

  /* Interpolate molecular density and abundance:                           */
  for(i=0; i<nmol; i++){
    splinterp(rsamp->n, rsamp->v, atms->molec[i].d,  nrad, rad->v, molec[i].d);
    splinterp(rsamp->n, rsamp->v, atms->molec[i].q,  nrad, rad->v, molec[i].q);
  }

  /* Interpolate isotopic partition function and cross section:             */
  for(i=0; i<ndb; i++){       /* For each database separately:              */
    iso1db = iso->db[i].s;    /* Index of first isotope in current DB       */
    isovs  = li->isov + iso1db;
    for(j=0; j < iso->db[i].i; j++){
      transitASSERT(iso1db + j > niso-1,
                    "Trying to reference an isotope (%i) outside the extended "
                    "limit (%i).\n", iso1db+j, niso-1);
      splinterp(li->db[i].t, li->db[i].T, isovs[j].z,
                nrad,        atmt->t,     iso->isov[iso1db+j].z);
    }
  }
  /* Set progress indicator and return:                                     */
  if(res>=0)
    tr->pi |= TRPI_MAKERAD;
  return res;
}


/* \fcnfh  DEF
   Calls makesample with the appropiate parameters and set the flags for
   an impact parameter sampling

   @returns makesample() output                                             */
int
makeipsample(struct transit *tr){
  struct transithint *th = tr->ds.th;
  int res=0; /* Return output */
  int i;

  /* Case when I used the atmospheric radius sampling:                      */
  if (th->rads.d == -1){
    tr->ips.n = tr->rads.n;
    tr->ips.d = 0;
    tr->ips.i = tr->rads.f;
    tr->ips.f = tr->rads.i;
    tr->ips.v = (PREC_RES *)calloc(tr->ips.n, sizeof(PREC_RES));
    for(i=0; i < tr->ips.n; i++)
      tr->ips.v[i] = tr->rads.v[tr->ips.n-i-1];
    tr->ips.o = 0;
    tr->ips.fct = tr->rads.fct;
  }
  /* FINDME: This is not even an option                                     */
  else{
    /* Hinted ip sample:                                                    */
    prop_samp usamp = {0, -th->ips.d,    th->ips.f,              th->ips.i,
                       th->ips.o,        NULL,                   th->ips.fct};
    /* Reference ip sample (inverse radius sample):                         */
    prop_samp rsamp = {0, -tr->rads.d,  tr->rads.v[tr->rads.n-1],
                       tr->rads.v[0],   tr->rads.o,   NULL,     tr->rads.fct};

    if(usamp.f < usamp.i) {
      tr_output(TOUT_ERROR,
        "Wrong specification of impact parameter, final value (%g) "
        "has to be bigger than initial (%g).\n", usamp.f, usamp.i);
      exit(EXIT_FAILURE);
    }

    transitcheckcalled(tr->pi, "makeipsample", 1, "makeradsample",TRPI_MAKERAD);

    /* Make the sampling taking as reference the radius sampling:           */
    res = makesample(&tr->ips, &usamp, &rsamp, TRH_IPRM);
  }

  /* Print sample information to file                                       */
  if (th->savefiles)
    outsample(tr);

  /* Set progress indicator:                                                */
  if(res >= 0)
    tr->pi |= TRPI_MAKEIP;
  return res;
}


/* \fcnfh  DEF
   Calls makesample with the appropiate parameters and set the flags for
   a temperature sampling

   Return: makesample() output                                              */
int
maketempsample(struct transit *tr){
  int res; /* Return output                                                 */
  struct transithint *th = tr->ds.th; /* transithint                        */

  /* Hinted temperature sample:                                             */
  prop_samp usamp = {0, th->temp.d, th->temp.i, th->temp.f, 1, NULL, 1};
  /* Reference temperature sample:                                          */
  prop_samp rsamp = {0, 0,          0,          0,          1, NULL, 1};

  if (usamp.f < usamp.i) {
    tr_output(TOUT_ERROR, "Wrong specification of temperature, final "
      "value (%g) has to be bigger than initial (%g).\n", usamp.f, usamp.i);
    exit(EXIT_FAILURE);
  }

  /* Make the sampling:                                                     */
  res = makesample(&tr->temp, &usamp, &rsamp, TRH_TEMP);

  /* Set progress indicator:                                                */
  if (res >= 0)
    tr->pi |= TRPI_MAKETEMP;
  return res;
}


/* FUNCTION
   Print sample info for a structure                                        */
static void
printsample(FILE *out,        /* File pointer to write out                  */
            prop_samp *samp,  /* Sample strucuture to print                 */
            char *desc,       /* Description                                */
            long fl){         /* Various flag                               */
  int i;

  /* Print header:                                                          */
  fprintf(out, "############################\n"
               "   %-12s Sampling\n"
               "----------------------------\n", desc);

  /* Print units, sample limits, and spacing:                               */
  fprintf(out, "Factor to cgs units: %g\n",            samp->fct);
  fprintf(out, "Initial value: %g\nFinal value: %g\n", samp->i, samp->f);
  fprintf(out, "Spacing: %g\n",                        samp->d);

  /* Print oversampling if necessary:                                       */
  if(!(fl&TRF_NOOVERSAMP))
    fprintf(out, "Oversample: %i\n", samp->o);

  /* Print sample size:                                                     */
  fprintf(out, "Number of elements: %lli\n", samp->n);

  /* Print sample array:                                                    */
  if(!(fl&TRF_NOVALUE)){
    fprintf(out, "Values: ");
    for(i=0; i<samp->n; i++)
      fprintf(out, " %12.8g", samp->v[i]);
    fprintf(out, "\n");
  }
}


/* \fcnfh  DEF
   Saves in binary the sample structure               */
void
savesample(FILE *out,        /* File pointer to write */
           prop_samp *samp){ /* Sampling strucuture   */
  fwrite(samp, sizeof(prop_samp), 1, out);
  savesample_arr(out, samp);
}


/* \fcnfh  DEF
   Saves in binary the sample structure's arrays */
void
savesample_arr(FILE *out,        /* File pointer to write */
               prop_samp *samp){ /* Sampling strucuture   */
  if(samp->n>0)
    fwrite(samp->v, sizeof(PREC_RES), samp->n, out);
}


/* \fcnfh  DEF
   Restore a binary sample structure

   Return: 0 on success, else
          -1 if not all the expected information is read
          -2 if info read is wrong
          -3 if cannot allocate memory
           1 if information read was suspicious      */
int
restsample(FILE *in,         /* File pointer to read */
           prop_samp *samp){ /* Sampling strucuture  */
  if(fread(samp, sizeof(prop_samp), 1, in) != 1)
    return -1;
  return restsample_arr(in, samp);
}


/* \fcnfh  DEF
   Restore a binary sample structure

   Return: 0 on success, else
          -1 if not all the expected information is read
          -2 if info read is wrong
          -3 if cannot allocate memory
           1 if information read was suspicious             */
int
restsample_arr(FILE *in,         /* File pointer to read */
               prop_samp *samp){ /* Sampling strucuture  */
  if(samp->n<0)
    return -2;
  if(samp->n>1000000)
    return 1;
  if((samp->v = (PREC_RES *)calloc(samp->n, sizeof(PREC_RES)))==NULL)
    return -3;
  if(samp->n==0)
    return 0;
  if(fread(samp->v, sizeof(PREC_RES), samp->n, in) != samp->n)
    return -1;

  return 0;
}


/* \fcnfh  DEF
   Print the sample data to file.

   Return: 0 on success, else
           1 if cannot write to file  */
int
outsample(struct transit *tr){
  FILE *out = stdout;
  char *filename = tr->f_outsample; /* Sample filename */

  /* Check output filename exists: */
  if(!filename)
    return 0;

  /* If default f_outsample, print to screen: */
  if(strcmp(filename, "-") != 0)
    if((out=fopen(filename, "w"))==NULL){
      tr_output(TOUT_WARN, "Cannot open file '%s' for writing sampling "
        "data.\n", filename);
      return 1;
    }

  tr_output(TOUT_INFO, "Printing sampling information in '%s'.\n\n",
    filename);

  /* Print each sample: */
  printsample(out, &tr->wns,  "Wavenumber",       TRF_NOVALUE);
  printsample(out, &tr->wavs, "Wavelength",       TRF_NOVALUE);
  printsample(out, &tr->rads, "Radius",           TRF_NOOVERSAMP);
  printsample(out, &tr->ips,  "Impact parameter", 0);

  fclose(out);

  return 0;
}


/* \fcnfh  DEF
 Frees the sampling structure */
void
freemem_samp(prop_samp *samp){
  if(samp->n)
    free(samp->v);
}


#ifdef DBGSAMPLE
int
main(int argc, char *argv[]){

  prop_samp lim, hint={0,0,0,0,0}, res;
  float mar=0;
  int i;

  if(argc<5){
    fprintf(stderr, "Syntax:\n"
            "    dbgsample <ini> <fin> <delt> <oversampling> [<margin>]\n");
    exit(0);
  }

  lim.i = atof(argv[1]);
  lim.f = atof(argv[2]);
  lim.d = atof(argv[3]);
  lim.o = atof(argv[4]);
  lim.n = 0;
  if(argc==6)
    mar = atof(argv[5]);

  i = makesample(&res, &hint, &lim, 0, 0, mar);

  fprintf(stderr, "Makesample returned %i\n"
          "Initial %g, final %g, delta %g, oversamp %i, number %li, "
          "margin %g\n", i, res.i, res.f, res.d, res.o, res.n, mar);

  if(i<0)
    exit(1);

  for(i=0; i<res.n; i++)
    fprintf(stderr," rad(%i): %g\n", i, res.v[i]);
}

#endif /* debug */
