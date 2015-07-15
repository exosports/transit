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

#include <transit.h>

/* FUNCTION: Wrapper to calculate a Voigt profile
   Return: 1/2 of the number of points in the profile                       */
int
getprofile(PREC_VOIGT **pr,  /* Pointer to 1D profile                       */
           PREC_RES dwn,     /* wavenumber spacing                          */
           PREC_VOIGT dop,   /* Doppler width                               */
           PREC_VOIGT lor,   /* Lorentz width                               */
           float ta,         /* times of alpha                              */
           int nwave){       /* Maximum half-size of profile                */

  PREC_VOIGTP bigalpha, /* Largest width (Doppler or Lorentz)               */
              wvgt;     /* Calculated half-width of profile                 */
  int nvgt,             /* Number of points in profile                      */
      j;                /* Auxiliary index                                  */

  /* Get the largest width (alpha Doppler or Lorentz):                      */
  bigalpha = dop;
  if(bigalpha<lor)
    bigalpha = lor;

  /* Half width from line center where the profile will be computed:        */
  wvgt = bigalpha*ta;
  /* Number of points in profile:                                           */
  nvgt = 2*(long)(wvgt/dwn+0.5) + 1;
  /* The profile needs to contain at least three wavenumber elements:       */
  if (nvgt < 2)
    nvgt = 3;
  /* Profile does not need to be larger than the wavenumber range:          */
  if (nvgt > 2*nwave)
    nvgt = 2*nwave + 1;

  /* Basic check that 'lor' or 'dop' are within sense:                      */
  if(nvgt < 0)
    transiterror(TERR_CRITICAL, "Number of Voigt bins (%d) are not positive.  "
                 "Doppler width: %g, Lorentz width: %g.\n", nvgt, dop, lor);

  /* Allocate profile array:                                                */
  *pr = (PREC_VOIGT *)calloc(nvgt, sizeof(PREC_VOIGT));

  /* Calculate voigt using a width that gives an integer number of 'dwn'
     spaced bins:                                                           */
  if((j=voigtn(nvgt, dwn*(long)(nvgt/2), lor, dop, pr, -1,
               nvgt > _voigt_maxelements ? VOIGT_QUICK:0)) != 1)
    transiterror(TERR_CRITICAL, "voigtn2() returned error code %i.\n", j);

  return nvgt/2;
}


/* FUNCTION:
   Saving extinction for a possible next run                                */
void
savefile_extinct(char *filename,
                 PREC_RES **e,
                 _Bool *c,
                 long nrad,
                 long nwav){
  FILE *fp;

  if((fp=fopen(filename, "w")) == NULL){
    transiterror(TERR_WARNING,
                 "Extinction savefile '%s' cannot be opened for writing.\n"
                 " Continuing without saving\n"
                 ,filename);
    return;
  }

  transitprint(2, verblevel, "Saving extinction file '%s'", filename);

  const char mn[] = "@E@S@";
  fwrite(mn, sizeof(char), 5, fp);
  fwrite(e[0], sizeof(PREC_RES), nrad*nwav, fp);
  fwrite(c, sizeof(_Bool), nrad, fp);

  fclose(fp);

  int i;
  for (i=0 ; i<nrad ; i++)
    if (c[i]) break;

  transitprint(2, verblevel, " done (%li/%li radii computed)\n", nrad-i, nrad);
}


/* FUNCTION:
   Restoring extinction for a possible next run                             */
void
restfile_extinct(char *filename,
                 PREC_RES **e,
                 _Bool *c,
                 long nrad,
                 long nwav){

  FILE *fp;

  if((fp=fopen(filename, "r")) == NULL){
    transiterror(TERR_WARNING,
                 "Extinction savefile '%s' cannot be opened for reading.\n"
                 "Continuing without restoring. You can safely ignore "
                 "this warning if this the first time you run for this "
                 "extinction savefile.\n", filename);
    return;
  }

  char mn[5];
  if(fread(mn, sizeof(char), 5, fp) != 5 || strncmp(mn,"@E@S@",5) != 0){
     transiterror(TERR_WARNING,
                  "Given filename for extinction savefile '%s' exists\n"
                  "and is not a valid extinction file. Remove it\n"
                  "before trying to use extinction savefile\n", filename);
     fclose(fp);
     return;
  }

  transitprint(2, verblevel, "Restoring extinction file '%s'", filename);

  fread(e[0], sizeof(PREC_RES), nrad*nwav, fp);
  fread(c,    sizeof(_Bool),    nrad,      fp);

  long i;
  for (i=0; i<nrad; i++)
    if (c[i]) break;

  transitprint(2, verblevel, " done (From the %lith radii)\n", i);

  fclose(fp);

}


/* FUNCTION:
   Fill up the extinction information in tr->ds.ex
   TD: Scattering parameters should be added at some point here.
  Return: 0 on success, else
          computeextradius()                                                */
int
extwn(struct transit *tr){
  struct transithint *th=tr->ds.th;
  static struct extinction st_ex;
  tr->ds.ex = &st_ex;
  struct extinction *ex = &st_ex;
  int i;

  /* Check these routines have been called: */
  transitcheckcalled(tr->pi, "extwn", 4,
                             "readinfo_tli",  TRPI_READINFO,
                             "readdatarng",   TRPI_READDATA,
                             "makewnsample",  TRPI_MAKEWN,
                             "makeradsample", TRPI_MAKERAD);
  transitacceptflag(tr->fl, tr->ds.th->fl, TRU_EXTBITS);

  struct isotopes *iso = tr->ds.iso;
  int niso = iso->n_i;
  int nrad = tr->rads.n;
  int nwn  = tr->wns.n;

  /* Check there is at least one atmospheric layer:                         */
  /* FINDME: Move to readatm                                                */
  if(tr->rads.n < 1){
    transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
                 "There are no atmospheric parameters specified. I need at "
                 "least one atmospheric point to calculate a spectra.\n");
    return -2;
  }
  /* Check there are at least two wavenumber sample points:                 */
  /* FINDME: Move to makewnsample                                           */
  if(tr->wns.n < 2){
    transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
                 "I need at least 2 wavenumber points to compute "
                 "anything; I need resolution.\n");
    return -3;
  }
  /* Check there is at least one isotope linelist                           */
  /* FINDME: This should not be a condition, we may want
     to calculate an atmosphere with only CIA for example.                  */
  if(niso < 1){
    transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
                 "You are requiring a spectra of zero isotopes!.\n");
    return -5;
  }

  /* Get the extinction coefficient threshold:                              */
  ex->ethresh = th->ethresh;

  /* Declare extinction-coefficient array:                                  */
  ex->e        = (PREC_RES **)calloc(nrad,     sizeof(PREC_RES *));
  if((ex->e[0] = (PREC_RES  *)calloc(nrad*nwn, sizeof(PREC_RES)))==NULL)
    transiterror(TERR_CRITICAL|TERR_ALLOC, "Unable to allocate %li = %li*%li "
                 "for the extinction coefficient.\n", nrad*nwn, nrad, nwn);

  for(i=1; i<nrad; i++){
    ex->e[i] = ex->e[0] + i*nwn;
  }

  /* Has the extinction been computed at given radius boolean:              */
  ex->computed = (_Bool *)calloc(nrad, sizeof(_Bool));

  /* Set progress indicator, and print and output extinction if one P,T
     was desired, otherwise return success:                                 */
  tr->pi |= TRPI_EXTWN;
  return 0;
}


/* FUNCTION:
   Free extinction coefficient structure arrays
   Return 0 on success */
int
freemem_extinction(struct extinction *ex, /* Extinciton struct       */
                   long *pi){             /* progress indicator flag */
  /* Free arrays: */
  free(ex->e[0]);
  free(ex->e);
  free(ex->computed);

  /* Update indicator and return: */
  *pi &= ~(TRPI_EXTWN);
  return 0;
}


/* FUNCTION:
   Restore hints structure, the structure needs to have been allocated
   before
   Returns 0 on success
          -1 if not all the expected information is read
          -2 if info read is wrong
          -3 if cannot allocate memory
           1 if information read was suspicious                             */
int
restextinct(FILE *in,
            PREC_NREC nrad,
            short niso,
            PREC_NREC nwn,
            struct extinction *ex){
  PREC_NREC nr,nw;
  short i;
  size_t rn;

  rn=fread(&nr,sizeof(PREC_NREC),1,in);
  if(rn!=1) return -1;
  rn=fread(&i,sizeof(short),1,in);
  if(rn!=1) return -1;
  rn=fread(&nw,sizeof(PREC_NREC),1,in);
  if(rn!=1) return -1;

  //no more than 10 000 isotopes, or 10 000 000 of radii or wavelenth
  if(nr!=nrad||nw!=nwn||i!=niso||
     niso>10000||nrad>10000000||nwn>10000000)
    return -2;

  if((ex->e    = (PREC_RES **)calloc(nrad,     sizeof(PREC_RES * ))) == NULL)
    return -3;
  if((ex->e[0] = (PREC_RES  *)calloc(nrad*nwn, sizeof(PREC_RES   ))) == NULL)
    return -3;

  rn = fread(ex->e[0], sizeof(PREC_RES), nwn*nrad,in);
  if(rn!=nwn*nrad) return -1;
  rn = fread(ex->computed, sizeof(PREC_RES), nrad, in);
  if(rn!=nrad) return -1;

  for(i=1; i<nrad; i++){
    ex->e[i] = ex->e[0] + i*nwn;
  }
  return 0;
}


/* FUNCTION: Compute the molecular extinction.
   Store results in kiso.  If permol is true, calculate extinction per
   molecule separately; else, collapse all extinction into kiso[0].         */
int
computemolext(struct transit *tr, /* transit struct                         */
              PREC_RES **kiso,    /* Extinction coefficient array [mol][wn] */
              PREC_ATM temp,      /* Temperature                            */
              PREC_ATM *density,  /* Density per species                    */
              double *Z,          /* Partition Function per isotope         */
              int permol){        /* Calculate the extinction per molecule  */

  /* Transit structures:                                                    */
  struct opacity    *op =tr->ds.op;
  struct isotopes   *iso=tr->ds.iso;
  struct molecules  *mol=tr->ds.mol;
  struct line_transition *lt=&(tr->ds.li->lt);

  PREC_NREC ln;
  int i, m=0,
      *idop, *ilor;
  long j, maxj, minj, offset;

  /* Voigt profile variables:                                               */
  PREC_VOIGT ***profile=op->profile;  /* Voigt profile                      */
  PREC_NREC **profsize=op->profsize;  /* Voigt-profile half-size            */
  double *aDop=op->aDop,          /* Doppler-width sample                   */
         *aLor=op->aLor;          /* Lorentz-width sample                   */
  int nDop=op->nDop,              /* Number of Doppler samples              */
      nLor=op->nLor;              /* Number of Lorentz samples              */

  PREC_NREC subw,
            nlines=tr->ds.li->n_l; /* Number of line transitions            */
  PREC_RES wavn, next_wn;
  double fdoppler, florentz, /* Doppler and Lorentz-broadening factors      */
         csdiameter;         /* Collision diameter                          */
  double propto_k;
  double *kmax, *kmin,       /* Maximum and minimum values of propto_k      */
         **ktmp;

  PREC_VOIGTP *alphal, *alphad;

  int niso = iso->n_i,        /* Number of isotopes in atmosphere           */
      nmol = mol->nmol,       /* Number of species in atmosphere            */
      Nmol;                   /* Number of species with line-transitions    */

  int iown, idwn;             /* Line-center indices                        */

  double maxwidth=0,   /* Maximum width between Lorentz and Doppler         */
         minwidth=1e5; /* Minimum width among isotopes in a Layer           */

  int ofactor;  /* Dynamic oversampling factor                              */

  long nadd  = 0, /* Number of co-added lines                               */
       nskip = 0, /* Number of skipped lines                                */
       neval = 0; /* Number of evaluated profiles                           */

  /* Wavenumber array variables:                                            */
  PREC_RES  *wn = tr->wns.v;
  PREC_NREC onwn = tr->owns.n,
            dnwn;

  /* Wavenumber sampling intervals:                                         */
  PREC_RES  dwn = tr->wns.d /tr->wns.o,   /* Output array                   */
           odwn = tr->owns.d/tr->owns.o,  /* Oversampling array             */
           ddwn;                          /* Dynamic sampling array         */

  /* Allocate alpha Lorentz and Doppler arrays:                             */
  alphal = (PREC_VOIGTP *)calloc(niso, sizeof(PREC_VOIGTP));
  alphad = (PREC_VOIGTP *)calloc(niso, sizeof(PREC_VOIGTP));

  /* Allocate width indices array:                                          */
  idop = (int *)calloc(niso, sizeof(int));
  ilor = (int *)calloc(niso, sizeof(int));

  kmax = (double *)calloc(op->Nmol, sizeof(double));
  kmin = (double *)calloc(op->Nmol, sizeof(double));

  /* Number of species in output array:                                     */
  if (permol)
    Nmol = op->Nmol;
  else
    Nmol = 1;

  /* Temporary extinction array:                                            */
  ktmp    = (double **)malloc(Nmol            * sizeof(double *));
  ktmp[0] = (double  *)malloc(Nmol*tr->owns.n * sizeof(double  ));
  for (i=1; i<Nmol; i++)
    ktmp[i] = ktmp[0] + tr->owns.n * i;

  /* Constant factors for line widths:                                      */
  fdoppler = sqrt(2*KB*temp/AMU) * SQRTLN2 / LS;
  florentz = sqrt(2*KB*temp/PI/AMU) / (AMU*LS);

  /* Calculate the isotope's widths for this layer:                         */
  for(i=0; i<niso; i++){
    /* Lorentz profile width:                                               */
    alphal[i] = 0.0;
    for(j=0; j<nmol; j++){
      /* Isotope's collision diameter:                                      */
      csdiameter = (mol->radius[j] + mol->radius[iso->imol[i]]);
      /* Line width:                                                        */
      alphal[i] += density[j]/mol->mass[j] * csdiameter * csdiameter *
                   sqrt(1/iso->isof[i].m + 1/mol->mass[j]);
    }
    alphal[i] *= florentz;

    /* Doppler profile width (divided by central wavenumber):               */
    alphad[i] = fdoppler / sqrt(iso->isof[i].m);

    /* Print Lorentz and Doppler broadening widths:                         */
    if(i <= 0)
      transitprint(1, verblevel, "Broadening (cm-1): Lorentz: %.5e, Doppler: "
              "%.5e (T=%.2f).\n", alphal[i], alphad[i]*wn[0], temp);

    maxwidth = fmax(alphal[i], alphad[i]*wn[0]); /* Max between Dop and Lor */
    minwidth = fmin(minwidth, maxwidth);

    /* Search for aDop and aLor indices for alphal[i] and alphad[i]:        */
    idop[i] = binsearchapprox(aDop, alphad[i]*wn[0], 0, nDop);
    ilor[i] = binsearchapprox(aLor, alphal[i],       0, nLor);
  }

  transitprint(10, verblevel, "Minimum width in layer: %.9f\n", minwidth);
  /* Set oversampling resolution:                                           */
  for (i=1; i < tr->ndivs; i++)
    if (tr->odivs[i]*(dwn/tr->owns.o) >= 0.5 * minwidth){
      break;
    }
  ofactor = tr->odivs[i-1];         /* Dynamic-sampling oversampling factor */
  ddwn    = odwn * ofactor;         /* Dynamic-sampling grid interval       */
  dnwn    = 1 + (onwn-1) / ofactor; /* Number of dynamic-sampling values    */
  transitprint(10, verblevel, "Dynamic-sampling grid interval: %.9f  "
               "(scale factor:%i)\n", ddwn, ofactor);
  transitprint(10, verblevel, "Number of dynamic-sampling values:%li\n",
                                dnwn);

  /* Determine the maximum and minimum line-strength per isotope:           */
  for(ln=0; ln<nlines; ln++){
    /* Wavenumber of line transition:                                       */
    wavn = 1.0 / (lt->wl[ln] * lt->wfct);
    /* Isotope ID of line:                                                  */
    i = lt->isoid[ln];
    /* Species index in output array:                                       */
    if (permol)
      m = valueinarray(op->molID, mol->ID[iso->imol[i]], op->Nmol);

    /* If it is beyond the lower limit, skip to next line transition:       */
    if ((wavn < tr->wns.i) || (wavn > tr->owns.v[onwn-1]))
      continue;

    /* Calculate the extinction coefficient except the broadening factor:   */
    propto_k = iso->isoratio[i]               *       /* Density            */
            SIGCTE     * lt->gf[ln]           *       /* Constant * gf      */
            exp(-EXPCTE*lt->efct*lt->elow[ln]/temp) * /* Level population   */
            (1-exp(-EXPCTE*wavn/temp))        /       /* Induced emission   */
            iso->isof[i].m                    /       /* Isotope mass       */
            Z[i];                                     /* Partition function */
    /* Maximum line strength among all transitions for each species:        */
    if (kmax[m] == 0){
      kmax[m] = kmin[m] = propto_k;
    } else{
      kmax[m] = fmax(kmax[m], propto_k);
      kmin[m] = fmin(kmin[m], propto_k);
    }
  }

  /* Compute the spectra, proceed for every line:                           */
  for(ln=0; ln<nlines; ln++){
    wavn = 1.0/(lt->wl[ln]*lt->wfct);
    i    = lt->isoid[ln];
    if (permol)
      m = valueinarray(op->molID, mol->ID[iso->imol[i]], op->Nmol);

    if ((wavn < tr->wns.i) || (wavn > tr->owns.v[onwn-1]))
      continue;

    /* Extinction coefficient (factors depending on the line transition):   */
    propto_k = lt->gf[ln]                              *
               exp(-EXPCTE*lt->efct*lt->elow[ln]/temp) *
               (1-exp(-EXPCTE*wavn/temp));

    /* Index of closest oversampled wavenumber:                             */
    iown = (wavn - tr->wns.i)/odwn;
    if (fabs(wavn - tr->owns.v[iown+1]) < fabs(wavn - tr->owns.v[iown]))
      iown++;

    /* Check if the next line falls on the same sampling index:             */
    while (ln != nlines-1 && lt->isoid[ln+1] == i){
      next_wn = 1.0/(lt->wl[ln+1]*lt->wfct);
      if (fabs(next_wn - tr->owns.v[iown]) < odwn){
        nadd++;
        ln++;
        /* Add the contribution from this line into the opacity:            */
        propto_k += lt->gf[ln]                                    *
                    exp(-EXPCTE * lt->efct * lt->elow[ln] / temp) *
                    (1-exp(-EXPCTE*next_wn/temp));
      }
      else
        break;
    }
    /* The rest of the factors:                                             */
    propto_k *= SIGCTE*iso->isoratio[i] / (iso->isof[i].m * Z[i]);

    /* If line is too weak, skip it:                                        */
    if (propto_k < tr->ds.th->ethresh * kmax[m]){
      nskip++;
      continue;
    }
    /* Multiply by the species density:                                     */
    if (permol == 0)
      propto_k *= density[iso->imol[i]];

    /* Index of closest (but not larger than) dynamic-sampling wavenumber:  */
    idwn = (wavn - tr->wns.i)/ddwn;

    transitprint(1000, 2, "own[nown:%li]=%.3f  (wf=%.3f)\n",
                          onwn, tr->owns.v[onwn-1], tr->wns.f);
    transitprint(1000, 2, "wavn=%.3f   own[%i]=%.3f\n",
                          wavn, iown, tr->owns.v[iown]);

    /* FINDME: de-hard code this threshold                                  */
    /* Update Doppler width according to the current wavenumber:            */
    if (alphad[i]*wavn/alphal[i] >= 1e-1){
      /* Recalculate index for Doppler width:                               */
      idop[i] = binsearchapprox(aDop, alphad[i]*wavn, 0, nDop);
    }

    /* Sub-sampling offset between center of line and dyn-sampled wn:       */
    subw = iown - idwn*ofactor;
    /* Offset between the profile and the wavenumber-array indices:         */
    offset = ofactor*idwn - profsize[idop[i]][ilor[i]] + subw;
    /* Range that contributes to the opacity:                               */
    /* Set the lower and upper indices of the profile to be used:           */
    minj = idwn - (profsize[idop[i]][ilor[i]] - subw) / ofactor;
    maxj = idwn + (profsize[idop[i]][ilor[i]] + subw) / ofactor;
    if (minj < 0)
      minj = 0;
    if (maxj > dnwn)
      maxj = dnwn;

    transitprint(1000, verblevel, "minj:%li  maxj:%li  subw:%li  offset:%li  "
                       "index1:%li\nf=np.array([",
                 minj, maxj, subw, offset, ofactor*minj - offset);
    /* Add the contribution from this line to the opacity spectrum:         */
    /* Adding in more complex but faster array indexing based on simpler
     * pointer arrithmatic                                                  */
    PREC_VOIGT * tmp_point = profile[idop[i]][ilor[i]];
    int beg_j = ofactor*minj - offset;
    for(j=minj; j<maxj; ++j){
      ktmp[m][j] += propto_k * tmp_point[beg_j];
      beg_j += ofactor;
    }
    neval++;
  }
  /* Downsample ktmp to the final sampling size:                            */
  for (m=0; m < Nmol; m++)
    downsample(ktmp[m], kiso[m], dnwn, tr->owns.o/ofactor);

  transitprint(9, verblevel, "Number of co-added lines:     %8li  (%5.2f%%)\n",
                             nadd,  nadd*100.0/nlines);
  transitprint(9, verblevel, "Number of skipped profiles:   %8li  (%5.2f%%)\n",
                             nskip, nskip*100.0/nlines);
  transitprint(9, verblevel, "Number of evaluated profiles: %8li  (%5.2f%%)\n",
                             neval, neval*100.0/nlines);

  /* Free allocated memory:                                                 */
  free(ktmp[0]);
  free(ktmp);
  free(alphal);
  free(alphad);
  free(idop);
  free(ilor);
  free(kmax);
  free(kmin);

  return 0;
}


/* Obtain the molecular extinction by interpolating the opacity grid at
   the specified atmospheric layer:                                         */
int
interpolmolext(struct transit *tr, /* transit struct                        */
               PREC_NREC r,        /* Radius index                          */
               PREC_RES **kiso){   /* Extinction coefficient array          */

  struct opacity    *op=tr->ds.op;  /* Opacity struct                       */
  struct molecules *mol=tr->ds.mol;

  long Nmol, Ntemp, Nwave;
  PREC_RES *gtemp;
  int       *gmol;
  int itemp, imol,
      i, m;   /* for-loop indices                                           */
  double ext; /* Interpolated extinction coefficient                        */

  /* Layer temperature:                                                     */
  PREC_ATM temp = tr->atm.t[r] * tr->atm.tfct;
  /* Gridded temperatures:                                                  */
  gtemp = op->temp;
  Ntemp = op->Ntemp;
  /* Gridded molecules list:                                                */
  gmol = op->molID;
  Nmol = op->Nmol;
  /* Wavenumber array size:                                                 */
  Nwave = op->Nwave;

  /* Interpolate:                                                           */
  /* Find index of grid-temperature immediately lower than temp:            */
  itemp = binsearchapprox(gtemp, temp, 0, Ntemp);
  if (temp < gtemp[itemp])
    itemp--;
  transitprint(30, verblevel, "Temperature: T[%i]=%.0f < %.2f < T[%.i]=%.0f\n",
               itemp, gtemp[itemp], temp, itemp+1, gtemp[itemp+1]);

  for (i=0; i < Nwave; i++){
    /* Add contribution from each molecule:                                 */
    for (m=0; m < Nmol; m++){
      /* Linear interpolation of the extinction coefficient:                */
      ext = (op->o[r][itemp  ][m][i] * (gtemp[itemp+1]-temp) +
             op->o[r][itemp+1][m][i] * (temp - gtemp[itemp]) ) /
                                                 (gtemp[itemp+1]-gtemp[itemp]);
      imol = valueinarray(mol->ID, gmol[m], mol->nmol);
      kiso[r][i] += mol->molec[imol].d[r] * ext;
    }
  }

  return 0;
}

/* \fcnfh
   Compute scatering contribution to extinction
*/
void
computeextscat(double *e,
               long n,
               struct extscat *sc,
               double *rad,
               double trad,
               double *temp,
               double tcft,
               double wn){
  long i;

  for(i=0; i<n; i++)
    e[i] = 0;
}


/* \fcnfh
   Compute cloud contribution to extinction                                 */
void
computeextcloud(double *e,
               long n,
               struct extcloud *cl,
               prop_samp *rad,
               double *temp,
               double tcft,
               double wn){
  long i;
  double *radius  = rad->v,  /* Atmospheric-model layer's radius            */
       rfct       = rad->fct,
       extinction = cl->cloudext,
       top        = cl->cloudtop * rfct,
       bottom     = cl->cloudbot * rfct,
       slope;

  /* If there are no clouds, set array to zero:                             */
  if(extinction == 0){
    memset(e, 0, n*sizeof(double));
    return;
  }

  /* Slope of extinction per altitude:                                      */
  slope = extinction / (bottom - top);

  /* Find radius index right below the cloud top layer:                     */
  for(i=n-1; i>=0; i--){
    if(radius[i]*rfct <= top)
      break;
    e[i] = 0;  /* The extinction is zero above the cloud                    */
  }

  /* Set cloud extinction between cloud top and bottom:                     */
  for(; i>=0; i--){
    if(radius[i] * rfct <= bottom)
      break;
    e[i] = slope * (radius[i]*rfct - top);
  }

  /* Keep constant extinction until the bottom:                             */
  for(; i>=0; i--)
      e[i] = extinction;
}
