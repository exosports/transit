/*
 * extinction.c   - Computes extinction coefficient. Component of the
 *                  Transit program.
 *
 * Copyright (C) 2004-2006 Patricio Rojo (pato@astro.cornell.edu)
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of version 2 of the GNU General 
 * Public License as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 * 02111-1307, USA.
 */

#include <transit.h>

/* \fcnfh
   Wrapper to calculate a Voigt profile 'pr'

   Return: 1/2 of the number of points in the non-oversampled profile       */
inline int
newprofile(PREC_VOIGT **pr,  /* Output 2D profile                           */
           int vf,           /* Number of fine resolution bins              */
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

  //transitprint(1, verblevel, "aDop=%.3g,  aLor=%.3g\nwvgt=%.4g,  "
  //                           "alpha=%.4g.\n", dop, lor, wvgt, bigalpha);

  /* Basic check that 'lor' or 'dop' are within sense:                      */
  if(nvgt < 0)
    transiterror(TERR_CRITICAL, "Number of Voigt bins (%d) are not positive.  "
                 "Doppler width: %g, Lorentz width: %g.\n", nvgt, dop, lor);

  /* Allocate profile array:                                                */
  //transitprint(1, verblevel, "size = %d.\n", nvgt*vf);
  *pr = (PREC_VOIGT *)calloc(nvgt*vf, sizeof(PREC_VOIGT));
  //transitprint(1, verblevel, "FLAG alloc.\n");
  for(j=0; j<vf; j++)
    pr[j] = pr[0] + j*nvgt;

  /* Calculate voigt using a width that gives an integer number of 'dwn'
     spaced bins:                                                           */
  if((j=voigtn(vf, nvgt, dwn*(long)(nvgt/2), lor, dop, pr, -1, 
               nvgt>_voigt_maxelements?VOIGT_QUICK:0)) != 1)
    transiterror(TERR_CRITICAL, "voigtn() returned error code %i.\n", j);

  return nvgt/2;
}


/* Re-worked version of newprofile that calculates a profile in a 1D array
   (voigtfine = 1 case)                                                     */
/* FINDME: Investigate why it doesn't work.                                 */
inline int
getprofile(PREC_VOIGT *pr,   /* Output 1D profile                           */
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

  //transitprint(1, verblevel, "size = %d.\n", nvgt*vf);
  /* Allocate profile array:                                                */
  pr = (PREC_VOIGT *)calloc(nvgt, sizeof(PREC_VOIGT));

  /* Calculate voigt using a width that gives an integer number of 'dwn'
     spaced bins:                                                           */
  if((j=voigtn2(nvgt, dwn*(long)(nvgt/2), lor, dop, pr, -1,
               nvgt > _voigt_maxelements ? VOIGT_QUICK:0)) != 1)
    transiterror(TERR_CRITICAL, "voigtn2() returned error code %i.\n", j);

  //for (j=0; j<nvgt; j++){
  //  transitprint(1, 2, "%.5e, ", pr[j]);
  //}
  //transitprint(1, 2, "])\n");

  return nvgt/2;
}


/* \fcnfh
 Saving extinction for a possible next run                                  */
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


/* \fcnfh
   Restoring extinction for a possible next run
*/
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


/*\fcnfh
  Fill up the extinction information in tr->ds.ex  
  TD: Scattering parameters should be added at some point here.

  Return: 0 on success, else
          computeextradius()  */
int
extwn(struct transit *tr){
  struct transithint *th=tr->ds.th;
  static struct extinction st_ex;
  tr->ds.ex = &st_ex;
  struct extinction *ex = &st_ex;
  int i, j;

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

  /* Set extinction-per-isotope boolean: */
  ex->periso = ((tr->fl&TRU_EXTINPERISO)==TRU_EXTINPERISO);

  /* Arrange the extinctions so that the order is [iso][rad][wn] */
  int nni = ex->periso ? niso:1;
  int nnr = ex->periso ? nrad:0;

  /* FINDME: Should this be nni instead of niso? */
  ex->e           = (PREC_RES ***)calloc(niso,         sizeof(PREC_RES **));
  ex->e[0]        = (PREC_RES  **)calloc(nni*nrad,     sizeof(PREC_RES *));
  if((ex->e[0][0] = (PREC_RES   *)calloc(nni*nrad*nwn, sizeof(PREC_RES)))==NULL)
    transiterror(TERR_CRITICAL|TERR_ALLOC,
                 "Unable to allocate %li = %li*%li*%li to calculate "
                 "extinction for every radii, %stry to shorten the  "
                 "wavenumber range.\n", nrad*nni*nwn, nrad, nni, nwn,
                 ex->periso?"try disabling exctinction per isotope "
                               "(option --no-per-iso), or ":"");
  for(i=0; i<niso; i++){
    ex->e[i] = ex->e[0] + i*nnr;
    if(!i || ex->periso)
      for(j=0; j<nrad; j++)
        ex->e[i][j] = ex->e[0][0] + nwn*(j + nnr*i);
  }

  /* Has the extinction been computed at given radius boolean: */
  ex->computed = (_Bool *)calloc(nrad, sizeof(_Bool));

  /* For each radius (index 'r'):                                           */
  transitprint(1, verblevel, "\nThere are %d radii samples.\n", nrad);

  /* Set progress indicator, and print and output extinction if one P,T
     was desired, otherwise return success:                                 */
  tr->pi |= TRPI_EXTWN;
  if(tr->rads.n == 1)
    printone(tr);
  return 0;
}


/* \fcnfh
   Printout for one P,T conditions                                          */
void
printone(struct transit *tr){
  int rn;
  FILE *out=stdout;

  /* Open file:                                                             */
  if(tr->f_out && tr->f_out[0] != '-')
    out = fopen(tr->f_out, "w");

  transitprint(1, verblevel, "\nPrinting extinction for one radius (at %gcm) "
                             "in '%s'\n", tr->rads.v[0],
                             tr->f_out?tr->f_out:"standard output");

  /* Print:                                                                 */
  fprintf(out, "#wavenumber[cm-1]   wavelength[nm]   extinction[cm-1]   "
               "cross-section[cm2]\n");
  for(rn=0; rn < tr->wns.n; rn++)
    fprintf(out, "%12.6f%14.6f%17.7g%17.7g\n", tr->wns.fct*tr->wns.v[rn],
            1/(tr->wavs.fct * tr->wns.v[rn] * tr->wns.fct),
            tr->ds.ex->e[0][0][rn],
            AMU*tr->ds.ex->e[0][0][rn] *
            tr->ds.iso->isof[0].m/tr->ds.mol->molec[tr->ds.iso->imol[0]].d[0]);

  exit(EXIT_SUCCESS);
}


/* \fcnfh
   Free extinction coefficient structure arrays

   Return 0 on success */
int
freemem_extinction(struct extinction *ex, /* Extinciton struct       */
                   long *pi){             /* progress indicator flag */
  /* Free arrays: */
  free(ex->e[0][0]);
  free(ex->e[0]);
  free(ex->e);
  free(ex->computed);

  /* Update indicator and return: */
  *pi &= ~(TRPI_EXTWN);
  return 0;
}


/* \fcnfh
   Restore hints structure, the structure needs to have been allocated
   before

   @returns 0 on success
            -1 if not all the expected information is read
            -2 if info read is wrong
            -3 if cannot allocate memory
            1 if information read was suspicious
*/
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

  int nni=ex->periso?niso:1;

  if((ex->e      =(PREC_RES ***)calloc(niso        ,sizeof(PREC_RES **)))
     ==NULL)
    return -3;
  if((ex->e[0]   =(PREC_RES **) calloc(nni*nrad    ,sizeof(PREC_RES * )))
     ==NULL)
    return -3;
  if((ex->e[0][0]=(PREC_RES *)  calloc(nrad*nni*nwn,sizeof(PREC_RES   )))
     ==NULL)
    return -3;

  rn=fread(ex->e[0][0],sizeof(PREC_RES),nwn*nni*nrad,in);
  if(rn!=nwn*nni*nrad) return -1;
  rn=fread(ex->computed,sizeof(PREC_RES),nrad,in);
  if(rn!=nrad) return -1;

  int nnr = ex->periso?nrad:0;
  for(i=0;i<niso;i++){
    ex->e[i]=ex->e[0]+i*nnr;
    if(!i||ex->periso)
      for(nr=0;nr<nrad;nr++)
        ex->e[i][nr]=ex->e[0][0] + nwn*( nr + nni*i );
  }

  return 0;
}


/* \fcnfh
   Frees voigt profile pointer arrays. Data array should already be free  */
void
freemem_localextinction(){
  /* profile[0][0] will always be freed after a succesfull completion
     of extradius: */
  //free(profile[0]);
  //free(profile);

  /* Free auxiliar variables that were allocated to number of isotopes: */
  //free(wa);
  //free(alphal);
  //free(ziso);
  //free(densiso);
}


/* Compute the molecular extinction:                                        */
int
computemolext(struct transit *tr, /* transit struct                         */
              PREC_NREC r,        /* Radius index                           */
              PREC_RES ***kiso){  /* Extinction coefficient array           */

  struct opacity *op=tr->ds.op;
  struct isotopes  *iso=tr->ds.iso; /* Isotopes struct                      */
  struct molecules *mol=tr->ds.mol;
  struct extinction *ex=tr->ds.ex;
  struct line_transition *lt=&(tr->ds.li->lt); /* Line transition struct    */

  PREC_NREC ln;
  int i, *idop, *ilor;
  long j, maxj, minj, offset;

  /* Voigt profile variables:                                               */
  PREC_VOIGT ****profile=op->profile; /* Voigt profile                      */
  //PREC_VOIGT ***vprofile=op->vprofile; /* Voigt profile                     */
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
  double kmax=0, kmin;       /* Maximum and minimum values of propto_k      */

  PREC_VOIGTP *alphal, *alphad;

  int niso = iso->n_i,        /* Number of isotopes                         */
      nmol = mol->nmol;       /* Number of species                          */

  int iown, idwn;             /* Line-center indices                        */
  double minwidth, maxwidth;  /* FINDME: For-test only */

  /* Temporal extinction array:                                             */
  double *ktmp = (double *)calloc(tr->owns.n, sizeof(double));
  int ofactor;  /* Dynamic oversampling factor                              */

  long nadd  = 0, /* Number of co-added lines                               */
       nskip = 0, /* Number of skipped lines                                */
       neval = 0; /* Number of evaluated profiles                           */

  /* Wavenumber array variables:                                            */
  PREC_RES  *wn = tr->wns.v;
  PREC_NREC // nwn = tr->wns.n,
            onwn = tr->owns.n,
            dnwn;

  /* Wavenumber sampling intervals:                                         */
  PREC_RES  dwn = tr->wns.d /tr->wns.o,   /* Output array                   */
           odwn = tr->owns.d/tr->owns.o,  /* Oversampling array             */
           ddwn;                          /* Dynamic sampling array         */
  /* Layer temperature:                                                     */
  PREC_ATM temp = tr->atm.t[r]*tr->atm.tfct;

  /* Constant factors for line widths:                                      */
  fdoppler = sqrt(2*KB*temp/AMU) * SQRTLN2 / LS;
  florentz = sqrt(2*KB*temp/PI/AMU) / (AMU*LS);

  /* Allocate alpha Lorentz and Doppler arrays:                             */
  alphal = (PREC_VOIGTP *)calloc(niso, sizeof(PREC_VOIGTP));
  alphad = (PREC_VOIGTP *)calloc(niso, sizeof(PREC_VOIGTP));

  /* Allocate width indices array:                                          */
  idop = (int *)calloc(niso, sizeof(int));
  ilor = (int *)calloc(niso, sizeof(int));

  maxwidth = 0.0;
  minwidth = 1e5;
  /* Calculate the isotope's widths for this layer:                         */
  for(i=0; i<niso; i++){
    /* Lorentz profile width:                                               */
    alphal[i] = 0.0;
    for(j=0; j<nmol; j++){
      /* Isotope's collision diameter:                                      */
      csdiameter = (mol->radius[j] + mol->radius[iso->imol[i]]);
      /* Line width:                                                        */
      alphal[i] += mol->molec[j].d[r]/mol->mass[j] * csdiameter * csdiameter *
                   sqrt(1/iso->isof[i].m + 1/mol->mass[j]);
    }
    alphal[i] *= florentz;

    /* Doppler profile width (divided by central wavenumber):               */
    alphad[i] = fdoppler/sqrt(iso->isof[i].m);

    /* Print Lorentz and Doppler broadening widths:                         */
    if(i <= 0)
      transitprint(1, verblevel, "Lorentz: %.9f, Doppler: %.9f broadening "
              "(T=%d, r=%li).\n", alphal[i], alphad[i]*wn[0], (int)temp, r);

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
  transitprint(100, verblevel, "Dynamic-sampling grid interval: %.9f  "
               "(scale factor:%i)\n", ddwn, ofactor);
  transitprint(100, verblevel, "Number of dynamic-sampling values:%li\n",
                                dnwn);

  /* Determine the maximum and minimum extinction-coefficient in this layer */
  for(ln=0; ln<nlines; ln++){
    /* Skip lines with strength lower than minelow:                         */
    if(lt->elow[ln] < tr->minelow)
      continue;

    /* Wavenumber of line transition:                                       */
    wavn = 1.0/(lt->wl[ln]*lt->wfct);
    /* Isotope ID of line:                                                  */
    i = lt->isoid[ln];

    /* If it is beyond the lower limit, skip to next line transition:       */
    if(wavn < tr->wns.i)
      continue;
    if(wavn > tr->owns.v[onwn-1])
      continue;

    /* Calculate the extinction coefficient except the broadening factor:   */
    propto_k = mol->molec[iso->imol[i]].d[r]*iso->isoratio[i] *  /* Density */
              SIGCTE     * lt->gf[ln]           *       /* Constant * gf    */
              exp(-EXPCTE*lt->efct*lt->elow[ln]/temp) * /* Level population */
              (1-exp(-EXPCTE*wavn/temp))        /       /* Induced emission */
              iso->isof[i].m                    /       /* Isotope mass     */
              iso->isov[i].z[r];                      /* Partition function */
    if (kmax == 0){
      kmax = kmin = propto_k;
    } else{
      kmax = fmax(kmax, propto_k);
      kmin = fmin(kmin, propto_k);
    }
  }

  /* Compute the spectra, proceed for every line:                           */
  for(ln=0; ln<nlines; ln++){
    wavn = 1.0/(lt->wl[ln]*lt->wfct); /* Wavenumber */
    i    = lt->isoid[ln];             /* Isotope ID */

    if(wavn < tr->wns.i)
      continue;
    if(wavn > tr->owns.v[onwn-1])
      continue;

    /* Extinction coefficient:                                              */
    propto_k = mol->molec[iso->imol[i]].d[r] * iso->isoratio[i] *
               SIGCTE * lt->gf[ln]                              *
               exp(-EXPCTE*lt->efct*lt->elow[ln]/temp)          *
               (1-exp(-EXPCTE*wavn/temp)) / iso->isof[i].m      /
               iso->isov[i].z[r];

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
        propto_k += mol->molec[iso->imol[i]].d[r] * iso->isoratio[i] *
                    SIGCTE * lt->gf[ln]                              *
                    exp(-EXPCTE * lt->efct * lt->elow[ln] / temp)    *
                    (1-exp(-EXPCTE*next_wn/temp)) / iso->isof[i].m   /
                    iso->isov[i].z[r];
      }
      else
        break;
    }

    /* If line is too weak, skip it:                                        */
    if (propto_k < ex->ethresh * kmax){
      nskip++;
      continue;
    }
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

    if (r== 100 && ln >= 1 && ln <= 19){
      transitprint(100, verblevel, "k=%.10e, d=%.4e, rat=%.4e, gf=%.4e, "
                   "elow=%.4e, T=%.4e, w=%.4e, m=%.4e, z=%.4e\n", propto_k,
                   mol->molec[iso->imol[i]].d[r], iso->isoratio[i], lt->gf[ln],
                   lt->elow[ln], temp, wavn, iso->isof[i].m, iso->isov[i].z[r]);
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
    for(j=minj; j<maxj; j++){
      ktmp[j] += propto_k * profile[idop[i]][ilor[i]][0][ofactor*j - offset];
      transitprint(1000, verblevel, "%.4e, ",
                            profile[idop[i]][ilor[i]][0][ofactor*j - offset]);
      //ktmp[j] += propto_k * vprofile [idop[i]] [ilor[i]] [j-offset];
    }
    neval++;
  }
  transitprint(10, verblevel, "Kmin: %.5e   Kmax: %.5e\n", kmin, kmax);
  /* Downsample ktmp to the final sampling size:                          */
  downsample(ktmp, kiso[i][r], dnwn, tr->owns.o/ofactor);
  transitprint(9, verblevel, "Number of co-added lines:     %8li  (%5.2f%%)\n",
                             nadd,  nadd*100.0/nlines);
  transitprint(9, verblevel, "Number of skipped profiles:   %8li  (%5.2f%%)\n",
                             nskip, nskip*100.0/nlines);
  transitprint(9, verblevel, "Number of evaluated profiles: %8li  (%5.2f%%)\n",
                             neval, neval*100.0/nlines);

  ex->computed[r] = 1;
  return 0;
}


/* Obtain the molecular extinction by interpolating the opacity grid at
   the specified atmospheric layer:                                         */
int
interpolmolext(struct transit *tr, /* transit struct                        */
               PREC_NREC r,        /* Radius index                          */
               PREC_RES ***kiso){  /* Extinction coefficient array          */

  struct opacity    *op=tr->ds.op;  /* Opacity struct                       */
  struct molecules *mol=tr->ds.mol;
  struct extinction *ex=tr->ds.ex;

  long Nmol, Ntemp, Nwave;
  PREC_RES *gtemp;
  int       *gmol;
  int itemp,
      i, m;   /* for-loop indices                                           */
  double ext; /* Interpolated extinction coefficient                        */

  double *ktmp = (double *)calloc(tr->owns.n, sizeof(double));
  /* Layer temperature:                                                     */
  //transitprint(1,2,"\n Radius index: %li\n", r);
  PREC_ATM temp = tr->atm.t[r] * tr->atm.tfct;
  /* Gridded temperatures:                                                  */  
  gtemp = op->temp;
  Ntemp = op->Ntemp;
  /* Gridded molecules list:                                                */
  gmol = op->molID;
  Nmol = op->Nmol;
  /* Wavenumber array size:                                                 */
  Nwave = op->Nwave;

  /* Interpolate:         */
  /* Find index of grid-temperature immediately lower than temp:            */
  itemp = binsearchapprox(gtemp, temp, 0, Ntemp);
  if (gtemp[itemp] < temp)
    itemp--;
  //transitprint(1, 2, "Temperature: T[%i]=%.0f < %.2f < T[%.i]=%.0f\n",
  //                   itemp, gtemp[itemp], temp, itemp+1, gtemp[itemp+1]);

  for (i=0; i < Nwave; i++){
    /* Add contribution from each molecule:                                 */
    for (m=0; m < Nmol; m++){
      /* Linear interpolation of the extinction coefficient:                */
      ext = (op->o[m][itemp  ][r][i] * (gtemp[itemp+1]-temp) + 
             op->o[m][itemp+1][r][i] * (temp - gtemp[itemp]) ) /
                                                 (gtemp[itemp+1]-gtemp[itemp]);
      /* FINDME: first index is the isotope ID.  Remove it                  */
      kiso[0][r][i] += mol->molec[gmol[m]].d[r] * ext;
    }
  }
  free(ktmp);
  ex->computed[r] = 1;
  return 0;
}


