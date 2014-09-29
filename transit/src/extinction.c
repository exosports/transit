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
  PREC_NREC **profsize=op->profsize;  /* Voigt-profile half-size            */
  double *aDop=op->aDop,          /* Doppler-width sample                   */
         *aLor=op->aLor;          /* Lorentz-width sample                   */
  int nDop=op->nDop,              /* Number of Doppler samples              */
      nLor=op->nLor;              /* Number of Lorentz samples              */

  PREC_NREC ilinewn, subw,
            nlines=tr->ds.li->n_l; /* Number of line transitions            */
  PREC_RES wavn;
  double fdoppler, florentz, /* Doppler and Lorentz-broadening factors      */
         csdiameter;         /* Collision diameter                          */
  double propto_k;

  PREC_VOIGTP *alphal, *alphad;

  int voigtfine = tr->voigtfine,  /* Voigt profile oversampling             */
      niso = iso->n_i,            /* Number of isotopes                     */
      nmol = mol->nmol;

  /* Wavenumber array variables:                                            */
  PREC_RES  *wn = tr->wns.v;
  PREC_NREC nwn = tr->wns.n;
  PREC_RES  dwn = tr->wns.d/tr->wns.o;

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
      if (i==0)
        transitprint(20, verblevel, "AlphaL[%li] = %.3e\n", j,
                    mol->molec[j].d[r]/mol->mass[j] * csdiameter * csdiameter *
                    sqrt(1/iso->isof[i].m + 1/mol->mass[j]));
    }
    alphal[i] *= florentz;

    /* Doppler profile width (divided by central wavenumber):               */
    alphad[i] = fdoppler/sqrt(iso->isof[i].m);

    /* Print Lorentz and Doppler broadening widths:                         */
    if(i <= 0)
      transitprint(1, verblevel, "Lorentz: %.9f, Doppler: %.9f broadening "
              "(T=%d, r=%li).\n", alphal[i], alphad[i]*wn[0], (int)temp, r);

    /* Search for aDop and aLor indices for alphal[i] and alphad[i]:        */
    idop[i] = binsearchapprox(aDop, alphad[i]*wn[0], 0, nDop);
    ilor[i] = binsearchapprox(aLor, alphal[i],       0, nLor);
  }

  /* Compute the spectra, proceed for every line:                           */
  for(ln=0; ln<nlines; ln++){
    /* Skip lines with strength lower than minelow:                         */
    if(lt->elow[ln] < tr->minelow)
      continue;

    /* Wavenumber of line transition:                                       */
    wavn = 1.0/(lt->wl[ln]*lt->wfct);

    /* If it is beyond the lower limit, skip to next line transition:       */
    if(wavn < tr->wns.i)
      continue;
    else
      /* Closest wavenumber index to transition, but no larger than:        */
      ilinewn = (wavn - tr->wns.i)/dwn;
    transitDEBUG(25, verblevel, "wavn: %g,  lgf: %g.\n", wavn, lt->gf[ln]);

    /* If it is beyond the last index, skip to next line transition:        */
    if(ilinewn >= nwn)
      continue;

    i = lt->isoid[ln]; /* Isotope ID of line                                */

    /* FINDME: de-hard code this threshold                                  */
    if (alphad[i]*wn[ilinewn]/alphal[i] >= 1e-1){
      /* Recalculate index for Doppler width:                               */
      idop[i] = binsearchapprox(aDop, alphad[i]*wn[ilinewn], 0, nDop);
    }

    /* Calculate opacity coefficient except the broadening factor
       nor the molecular abundance:                                         */
    propto_k = mol->molec[iso->imol[i]].d[r] * iso->isoratio[i] * /* Density */
               SIGCTE     * lt->gf[ln]           *    /* Constant * gf      */
               exp(-EXPCTE*lt->efct*lt->elow[ln]/temp) * /* Level population */
               (1-exp(-EXPCTE*wavn/temp))        /    /* Induced emission   */
               iso->isof[i].m                    /    /* Isotope mass       */
               iso->isov[i].z[r];                     /* Partition function */

    if (r== 100 && ln >= 1 && ln <= 19){
      transitprint(100, verblevel, "k=%.10e, d=%.4e, rat=%.4e, gf=%.4e, "
                   "elow=%.4e, T=%.4e, w=%.4e, m=%.4e, z=%.4e\n", propto_k,
                   mol->molec[iso->imol[i]].d[r], iso->isoratio[i], lt->gf[ln],
                   lt->elow[ln], temp, wavn, iso->isof[i].m, iso->isov[i].z[r]);
    }

    /* Set the lower and upper indices of the profile to be used:           */
    minj = ilinewn - profsize[idop[i]][ilor[i]] + 1;
    offset = minj;
    if(minj < 0)
      minj = 0;
    maxj = ilinewn + profsize[idop[i]][ilor[i]] + 1;
    if(maxj > nwn)
      maxj = nwn;

    /* Distance from center of line to sampled wavenumber by 'ilinewn',
       in number of fine bins:                                              */
    subw = voigtfine*(wavn - ilinewn*dwn - tr->wns.i)/dwn;

    /* Add the contribution from this line to the opacity spectrum:         */
    for(j=minj; j<maxj; j++){
      kiso[i][r][j] += propto_k * profile[idop[i]][ilor[i]][subw][j-offset];
    }

  }
  ex->computed[r] = 1;
  return 0;
}
