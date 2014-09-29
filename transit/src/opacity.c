/*
 * opacity.c - Calculate the opacity for a grid of temperature and pressures
 *             (or radius?) and store to file.
 *
 * Copyright (C) 2014 Patricio Cubillos
 *

List of functions defined:

int opacity(struct transit *tr)
   Driver routine to calculate or read the opacity.

int calcopacity(struct transit *tr, FILE *fp)
   Calculate a grid of opacities and Voigt profiles.

int readopacity(struct transit *tr, FILE *fp)
   Read a opacity grid.

int extinction(struct transit *tr, int r, int t)
   Calculate the opacity spectrum at an specific layer.

int freemem_opacity(struct opacity *op, long *pi)
   Free index of refraction array.

*/


#include <transit.h>

/* \fcnfh
   Calculate the opacity due to molecular transitions.

   Return: 0 on success                                                     */
int
opacity(struct transit *tr){
  struct transithint *th = tr->ds.th; /* transithint struct                 */
  static struct opacity op;           /* The opacity struct                 */

  /* Set the opacity struct's mem to 0:                                     */
  memset(&op, 0, sizeof(struct opacity));
  tr->ds.op = &op;

  /* Check that the radius array has been sampled:                          */
  transitcheckcalled(tr->pi, "opacity", 1, "makeradsample", TRPI_MAKERAD);

  /* Check if the opacity file was specified:                               */
  if(th->f_opa == NULL){
    transitprint(1, verblevel, "No opacity file specified.\n");
  }

  /* Check if the opacity file exists:                                      */
  int fe = fileexistopen(th->f_opa, &tr->fp_opa);
  tr->f_opa = th->f_opa;  /* Set file name in transit struct                */
  transitprint(10, verblevel, "File-exist status = %d\n", fe);

  /* If the file exists, read the opacities from it:                        */
  if(fe == 1){
    transitprint(1, verblevel, "Reading opacity file: '%s'.\n", tr->f_opa);
    /* Read the grid of opacities from file:                                */
    readopacity(tr, tr->fp_opa);
  }
  /* Else, calculate Voigt profiles and the opacity grid:                   */
  else{
    /* Open file for writing:                                               */
    tr->fp_opa = fopen(tr->f_opa, "wb");
    /* File given but cannot be opened:                                     */
    if(fe == -1 && tr->fp_opa == NULL){
        transiterror(TERR_WARNING, "Opacity filename '%s' cannot be opened for "
                                   "writing.\n", tr->f_opa);
        return -1;
    }
    /* Calculate the Voigt profiles (and grid of opacities if requested):   */
    if (tr->fp_opa != NULL)
      transitprint(1, verblevel, "Calculating new grid of opacities: '%s'.\n",
                               tr->f_opa);
    else
      transitprint(1, verblevel, "Calculating grid of Voigt profiles.\n");

    calcopacity(tr, tr->fp_opa);
  }

  /* Set progress indicator and return success:                             */
  tr->pi |= TRPI_OPACITY;
  return 0;
}


/* Calculate opacities for the grid of wavenumber, radius, and temperature
   arrays for each molecule.                                                */
int
calcopacity(struct transit *tr,
            FILE *fp){
  struct opacity *op=tr->ds.op;     /* Opacity struct                       */
  struct isotopes  *iso=tr->ds.iso; /* Isotopes struct                      */
  struct lineinfo *li=tr->ds.li;    /* Lineinfo struct                      */
  long Nmol, Ntemp, Nlayer, Nwave,  /* Opacity-grid  dimension sizes        */
       flag;                        /* Interpolation flag                   */
  int i, j, k, t, r,                /* for-loop indices                     */
      rn, iso1db;
  int nDop, nLor;            /* Number of Doppler and lorentz-width samples */
  double Lmin, Lmax, Dmin, Dmax;   /* Minimum and maximum widths            */
  PREC_VOIGT ****profile;          /* Grid of Voigt profiles                */
  int   voigtfine =tr->voigtfine;  /* Voigt profile wn oversampling factor  */
  float timesalpha=tr->timesalpha; /* Voigt wings width                     */

  struct timeval tv;  /* Time-keeping variables                             */
  double t0=0.0;

  /* Run transitcheckcalled: TBD                                            */

  /* Make logscale grid for the profile widths:                             */
  /* FINDME: Hardcoded values                                               */
  nDop = op->nDop = 20;  //100;
  nLor = op->nLor = 20;  //100;
  Lmin = 1e-4;
  Lmax = 1.0; //10.0;
  Dmin = 1e-3;
  Dmax = 0.25;
  op->aDop = logspace(Dmin, Dmax, nDop);
  op->aLor = logspace(Lmin, Lmax, nLor);

  /* Allocate array for the profile half-size:                              */
  op->profsize    = (PREC_NREC **)calloc(nDop,      sizeof(PREC_NREC *));
  op->profsize[0] = (PREC_NREC  *)calloc(nDop*nLor, sizeof(PREC_NREC));
  for (i=1; i<nDop; i++)
    op->profsize[i] = op->profsize[0] + i*nLor;

  /* Allocate grid of Voigt profiles:                                       */
  op->profile       = (PREC_VOIGT ****)calloc(nDop,     sizeof(PREC_VOIGT ***));
  op->profile[0]    = (PREC_VOIGT  ***)calloc(nDop*nLor, sizeof(PREC_VOIGT **));
  op->profile[0][0] = (PREC_VOIGT   **)calloc(nDop*nLor*voigtfine,
                                                         sizeof(PREC_VOIGT *));
  for (i=0; i<nDop; i++){
    op->profile[i] = op->profile[0] + i*nLor;
    for (j=0; j<nLor; j++){
      op->profile[i][j] = op->profile[0][0] + voigtfine*(j + i*nLor);
    }
  }
  profile = op->profile;
  transitprint(10, verblevel, "Number of Voigt profiles: %d.\n",
                              nDop*nLor*voigtfine);

  t0 = timestart(tv, "Begin Voigt profiles calculation.");
  /* Evaluate the profiles for the array of widths:                         */
  for   (i=0; i<nDop; i++){
    for (j=0; j<nLor; j++){
      /* Skip calculation if Doppler width << Lorentz width:                */
      /* Set size and pointer to previous profile:                          */
      if (op->aDop[i]*10.0 < op->aLor[j]  &&  i != 0){
        op->profsize[i][j] = op->profsize[i-1][j];
        for (k=0; k<voigtfine; k++)
          profile[i][j][k] = profile[i-1][j][k];
      }
      else{ /* Calculate a new profile for given widths:                    */
        op->profsize[i][j] = newprofile(profile[i][j], voigtfine, tr->wns.d,
                               op->aDop[i], op->aLor[j], timesalpha, tr->wns.n);
      }
      transitprint(5, verblevel, "Profile[%2d][%2d] size = %4li  (D=%.3g, "
                                  "L=%.3g).\n", i, j, 2*op->profsize[i][j]+1,
                                   op->aDop[i], op->aLor[j]);
    }
  }
  t0 = timecheck(verblevel, 0, 0, "End Voigt-profile calculation", tv, t0);

  /* Make temperature array from hinted values:                             */
  maketempsample(tr);
  Ntemp = op->Ntemp = tr->temp.n;
  op->temp = (PREC_RES *)calloc(Ntemp, sizeof(PREC_RES));
  for (i=0; i<Ntemp; i++)
    op->temp[i] = tr->temp.v[i];
  transitprint(1, verblevel, "There are %li temperature samples.\n", Ntemp);

  /* Evaluate the partition at these temperatures:                          */
  op->ziso    = (PREC_ATM **)calloc(iso->n_i,       sizeof(PREC_ATM *));
  op->ziso[0] = (PREC_ATM  *)calloc(iso->n_i*Ntemp, sizeof(PREC_ATM));
  for(i=1; i<iso->n_i; i++)
    op->ziso[i] = op->ziso[0] + i*Ntemp;

  /* Set interpolation function flag:                                       */
  flag = tr->interpflag;
  flag = 1;  /* FINDME: Temporary hack                                      */
  /* FINDME: resample if throwing a segfault for flag=1 (spline interp), I
             can't see what's going on.  Must be on the internals of the
             routine.  Side note: Pato also found some fishy behavior around
             resample (in makeradsample).                                   */
  /* Interpolate isotopic partition function:                               */
  for(i=0; i<iso->n_db; i++){  /* For each database separately:             */
    iso1db = iso->db[i].s;     /* Index of first isotope in current DB      */

    resamplex(flag, li->db[i].t, li->db[i].T, Ntemp, op->temp);
    for(j=0; j < iso->db[i].i; j++){
      transitASSERT(iso1db + j > iso->n_i-1, "Trying to reference an isotope "
             "(%i) outside the extended limit (%i).\n", iso1db+j, iso->n_i-1);
      resampley(flag, 1, li->isov[iso1db+j].z, op->ziso[iso1db+j]);
    }
  }
  resample_free();

  /* Get pressure array from transit (save in CGS units):                   */
  Nlayer = op->Nlayer = tr->rads.n;
  op->press = (PREC_RES *)calloc(Nlayer, sizeof(PREC_RES));
  for (i=0; i<Nlayer; i++)
    op->press[i] = tr->atm.p[i]*tr->atm.pfct;
  transitprint(1, verblevel, "There are %li radius samples.\n", Nlayer);

  /* Make molecules array from transit:                                     */
  Nmol = op->Nmol = tr->ds.iso->nmol;
  op->molID = (int *)calloc(Nmol, sizeof(int));
  transitprint(1, verblevel, "There are %li molecules with line "
                             "transitions.\n", Nmol);
  for (i=0, j=0; i<iso->n_i; i++){
    /* If this molecule is not yet in molID array, add it:                  */
    if (valueinarray(op->molID, iso->imol[i], j) < 0){
      op->molID[j++] = iso->imol[i];
      transitprint(10, verblevel, "Isotope's (%d) molecule ID: %d added at "
                   "position %d.\n", i, iso->imol[i], j-1);
    }
  }

  /* Get wavenumber array from transit:                                     */
  Nwave = op->Nwave = tr->wns.n;
  op->wns = (PREC_RES *)calloc(Nwave, sizeof(PREC_RES));
  for (i=0; i<Nwave; i++)
    op->wns[i] = tr->wns.v[i];
  transitprint(1, verblevel, "There are %li wavenumber samples.\n", Nwave);

  /* Allocate opacity array:                                                */
  if (fp != NULL){
    op->o      = (PREC_RES ****)       calloc(Nmol,   sizeof(PREC_RES ***));
    for (i=0; i<Nmol; i++){
      op->o[i] = (PREC_RES  ***)       calloc(Ntemp,  sizeof(PREC_RES **));
      for (t=0; t<Ntemp; t++){
        op->o[i][t] = (PREC_RES **)    calloc(Nlayer, sizeof(PREC_RES *));
        for (r=0; r<Nlayer; r++){
          op->o[i][t][r] = (PREC_RES *)calloc(Nwave,  sizeof(PREC_RES));
        }
      }
    }

    if (!op->o[0][0][0])
      transitprint(1, verblevel, "Allocation fail\n");

    /* Compute extinction:                                                  */
    for (r=0;   r<Nlayer; r++){  /* For each layer:                         */
      for (t=0; t<Ntemp;  t++){  /* For each temperature:                   */
      if((rn=extinction(tr, r, t)) != 0)
        transiterror(TERR_CRITICAL, "extinction() returned error code %i.\n",
                                    rn);
      }
    }

    /* Save dimension sizes:                                                */
    fwrite(&Nmol,   sizeof(long), 1, fp);
    fwrite(&Ntemp,  sizeof(long), 1, fp);
    fwrite(&Nlayer, sizeof(long), 1, fp);
    fwrite(&Nwave,  sizeof(long), 1, fp);

    /* Save arrays:                                                         */
    fwrite(&op->molID[0], sizeof(int),      Nmol,   fp);
    fwrite(&op->temp[0],  sizeof(PREC_RES), Ntemp,  fp);
    fwrite(&op->press[0], sizeof(PREC_RES), Nlayer, fp);
    fwrite(&op->wns[0],   sizeof(PREC_RES), Nwave,  fp);
    
    /* Save opacity:                                                        */
    for (i=0; i<Nmol; i++)
      for (t=0; t<Ntemp; t++)
        for (r=0; r<Nlayer; r++)
          fwrite(op->o[i][t][r], sizeof(PREC_RES), Nwave, fp);

    fclose(fp);
  }
  transitprint(2, verblevel, "Done.\n"); 
  return 0;
}


/* Read the opacity file and store values in the transit structure.         */
int
readopacity(struct transit *tr,  /* transit struct                          */
            FILE *fp){           /* Pointer to file to read                 */
  struct opacity *op=tr->ds.op;  /* opacity struct                          */
  int i, t, r;  /* for-loop indices                                         */

  transitprint(1, verblevel, "Reading opacity file:\n");
  /* Read file dimension sizes:                                             */
  fread(&op->Nmol,   sizeof(long), 1, fp);
  fread(&op->Ntemp,  sizeof(long), 1, fp);
  fread(&op->Nlayer, sizeof(long), 1, fp);
  fread(&op->Nwave,  sizeof(long), 1, fp);
  transitprint(1, verblevel, "Opacity grid size: Nmol=%li   Ntemp=%li   "
               "Nlayer=%li   Nwave=%li.\n", op->Nmol, op->Ntemp, op->Nlayer,
                                            op->Nwave);
  transitprint(10, verblevel, "ftell=%li\n", ftell(fp));

  /* Allocate and read arrays:                                              */
  op->molID = (int      *)calloc(op->Nmol,   sizeof(int));
  op->temp  = (PREC_RES *)calloc(op->Ntemp,  sizeof(PREC_RES));
  op->press = (PREC_RES *)calloc(op->Nlayer, sizeof(PREC_RES));
  op->wns   = (PREC_RES *)calloc(op->Nwave,  sizeof(PREC_RES));
  fread(op->molID, sizeof(int),      op->Nmol,   fp);
  fread(op->temp,  sizeof(PREC_RES), op->Ntemp,  fp);
  fread(op->press, sizeof(PREC_RES), op->Nlayer, fp);
  fread(op->wns,   sizeof(PREC_RES), op->Nwave,  fp);

  /* DEBUGGING: Print temperature array                                     */
  for (i=0; i<op->Ntemp; i++)
    transitprint(10, verblevel, "T[%d]=%4d  ", i, (int)op->temp[i]);
  transitprint(10, verblevel, "\n");

  /* Allocate and read the opacity grid:                                    */
  //op->o      = (PREC_RES ****)calloc(op->Nmol,          sizeof(PREC_RES ***));
  //op->o[0]   = (PREC_RES  ***)calloc(op->Nmol*op->Ntemp, sizeof(PREC_RES **));
  //op->o[0][0] = (PREC_RES  **)calloc(op->Nmol*op->Ntemp*op->Nlayer,
  //                                                        sizeof(PREC_RES *));
  //op->o[0][0][0] = (PREC_RES *)calloc(op->Nmol*op->Ntemp*op->Nlayer*op->Nwave,
  //                                                        sizeof(PREC_RES));
  //for (i=0; t<op->Nmol; i++){
  //  op->o[i] = op->o[0] + op->Ntemp*i;
  //  for (t=0; t<op->Ntemp; t++){
  //    op->o[i][t] = op->o[0][0] + op->Nlayer*(t + op->Ntemp*i);
  //    for (r=0; r<op->Nlayer; r++)
  //      op->o[i][t][r] = op->o[0][0][0] + op->Nwave*(r + op->Nlayer*(t + op->Ntemp*i));
  //  }
  //}
  //fread(op->o[0][0][0], sizeof(PREC_RES), op->Nmol*op->Ntemp*op->Nlayer*op->Nwave, fp);

  op->o      = (PREC_RES ****)       calloc(op->Nmol,   sizeof(PREC_RES ***));
  for (i=0; i<op->Nmol; i++){
    op->o[i] = (PREC_RES  ***)       calloc(op->Ntemp,  sizeof(PREC_RES **));
    for (t=0; t<op->Ntemp; t++){
      op->o[i][t] = (PREC_RES **)    calloc(op->Nlayer, sizeof(PREC_RES *));
      for (r=0; r<op->Nlayer; r++){
        op->o[i][t][r] = (PREC_RES *)calloc(op->Nwave,  sizeof(PREC_RES));
      }
    }
  }

  /* Read the opacity grid:                                                 */
  for (i=0; i<op->Nmol; i++)
    for (t=0; t<op->Ntemp; t++)
      for (r=0; r<op->Nlayer; r++)
        fread(op->o[i][t][r], sizeof(PREC_RES), op->Nwave, fp);

  return 0;
}


/* Calculate the extinction coefficient/density for the given pressure and
   temperature.                                                             */
int
extinction(struct transit *tr,      /* transit struct                       */
           int r,                   /* Layer index                          */
           int t){                  /* Temperature index                    */

  struct opacity   *op =tr->ds.op;  /* Opacity struct                       */
  struct isotopes  *iso=tr->ds.iso; /* Isotopes struct                      */
  struct molecules *mol=tr->ds.mol; /* Molecules struct                     */
  struct atm_data  *at =tr->ds.at;  /* Atmosphere struct                    */
  struct line_transition *lt=&(tr->ds.li->lt); /* Line transition struct    */ 

  /* Voigt profile variables:                                               */
  PREC_VOIGT ****profile=op->profile; /* Voigt profile                      */
  PREC_NREC **profsize=op->profsize;  /* Voigt-profile half-size            */
  double *aDop=op->aDop,              /* Doppler-width sample               */
         *aLor=op->aLor;              /* Lorentz-width sample               */
  int nDop=op->nDop,                  /* Number of Doppler samples          */
      nLor=op->nLor;                  /* Number of Lorentz samples          */

  double fdoppler, florentz, /* Doppler and Lorentz-broadening factors      */
         csdiameter,         /* Collision cross-section diameter            */
         opac;               /* Extinction coefficient                      */

  int i, j, m,      /* for-loop indices                                     */
      *idop, *ilor, /* Doppler and Lorentz width indices per isotope        */
      niso, nmol,   /* Number of isotopes and molecules                     */
      nwn, neval;

  long minj, maxj,  /* Indices of min and max wn array of line profile      */
       offset,      /* Index offset between profile and wn arrays           */
       ncalc;
  PREC_ATM moldensity,          /* Molecular density                        */
           temp;                /* Layer temperature                        */
  PREC_VOIGTP *alphal, *alphad; /* Lorentz and Doppler widths               */
  PREC_RES *wn, dwn, wavn;      /* Wavenumber sampling info                 */
  PREC_NREC ilinewn,  /* Wavenumber index for line transition               */
            nlines,   /* Number of line transitions                         */
            ln,       /* line-transition for-loop index                     */
            subw;     /* Fine binning index for profile of line transition  */
  int voigtfine  = tr->voigtfine;  /* Voigt profile oversampling            */  

  struct timeval tv;
  double t0 = timestart(tv, "\nOpacity calculation:");
  //transitprint(1, verblevel, "FLAG 210.\n");

  /* Set up of variables:                                                   */

  /* Layer temperature:                                                     */
  temp = op->temp[t];   /* DIFF */

  niso = iso->n_i;     /* Number of isotopes                                */
  nmol = mol->nmol;    /* Number of molecules in atmosphere                 */

  wn  = op->wns;    /* Wavenumber array                                     */
  nwn = op->Nwave;  /* Number of wavenumber samples                         */
  dwn = tr->wns.d;  /* Wavenumber sampling stepsize                         */

  nlines = tr->ds.li->n_l; /* Number of line transitions                    */

  /* Allocate alpha Lorentz and Doppler arrays:                             */
  alphal = (PREC_VOIGTP *)calloc(niso, sizeof(PREC_VOIGTP));
  alphad = (PREC_VOIGTP *)calloc(niso, sizeof(PREC_VOIGTP));

  /* Allocate width indices array:                                          */
  idop = (int *)calloc(niso, sizeof(int));
  ilor = (int *)calloc(niso, sizeof(int));

  /* Constant factors for line widths:                                      */
  fdoppler = sqrt(2*KB*temp/AMU) * SQRTLN2 / LS;
  florentz = sqrt(2*KB*temp/PI/AMU) / (AMU*LS);
  //transitprint(1, verblevel, "FLAG 230.\n");

  /* Calculate the isotope's widths for this layer:                         */
  for(i=0; i<niso; i++){
    /* Lorentz profile width:                                               */
    alphal[i] = 0.0;
    for(j=0; j<nmol; j++){
      /* Isotope's collision diameter:                                      */
      csdiameter = (mol->radius[j] + mol->radius[iso->imol[i]]);
      /* Density of molecule colliding with current isotope:                */
      moldensity = stateeqnford(at->mass, mol->molec[j].q[r], tr->atm.mm[r],
                         mol->mass[j], op->press[r], temp);
      /* Line width:                                                        */
      alphal[i] += moldensity/mol->mass[j] * csdiameter * csdiameter *
                   sqrt(1/iso->isof[i].m + 1/mol->mass[j]);
      if (i==0)
        transitprint(20, verblevel, "AlphaL[%d] = %.3e\n", j,
                    moldensity/mol->mass[j] * csdiameter * csdiameter *
                    sqrt(1/iso->isof[i].m + 1/mol->mass[j]));
    }
    alphal[i] *= florentz;

    /* Doppler profile width (divided by central wavenumber):               */
    alphad[i] = fdoppler/sqrt(iso->isof[i].m);

    /* Print Lorentz and Doppler broadening widths:                         */
    if(i <= 0)
      transitprint(1, verblevel, "Lorentz: %.9f cm-1, Doppler: %.9f cm-1 "
        "broadening (T=%d, r=%d).\n", alphal[i], alphad[i]*wn[0], (int)temp, r);

    /* Search for aDop and aLor indices for alphal[i] and alphad[i]:        */
    idop[i] = binsearchapprox(aDop, alphad[i]*wn[0], 0, nDop);
    ilor[i] = binsearchapprox(aLor, alphal[i],       0, nLor);
  }

  t0 = timestart(tv, "Begin opacity calculation:");
  neval = 4;
  ncalc = 0;
  /* Compute the opacity, proceed for every line:                           */
  for(ln=0; ln<nlines; ln++){
    /* Skip lines with strength lower than transit's minelow threshold:     */
    if(lt->elow[ln] < tr->minelow)
      continue;

    /* Wavenumber of line transition:                                       */
    wavn = 1.0 / (lt->wl[ln] * lt->wfct);

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

    /* Isotope ID of line transition:                                       */
    i = lt->isoid[ln];
    /* Molecule index in opacity grid for this isotope:                     */
    m = valueinarray(op->molID, iso->imol[i], op->Nmol);

    /* FINDME: de-hard code this threshold                                  */
    if (alphad[i]*wn[ilinewn]/alphal[i] >= 1e-1){
    neval++;
      /* Recalculate index for Doppler width:                               */
      idop[i] = binsearchapprox(aDop, alphad[i]*wn[ilinewn], 0, nDop);
    }

    /* Calculate opacity coefficient except the broadening factor
       and the molecular abundance:                                         */
    opac = iso->isoratio[i]                        * /* Density             */
           SIGCTE     * lt->gf[ln]                 * /* Constant * gf       */
           exp(-EXPCTE*lt->efct*lt->elow[ln]/temp) * /* Level population    */
           (1-exp(-EXPCTE*wavn/temp))              / /* Induced emission    */
           iso->isof[i].m                          / /* Isotope mass        */
           op->ziso[i][t];                           /* Partition function  */
 
    transitprint(24, verblevel,
                 "i=%i   temp=%g   Elow=%g\n"
                 "aD=%.7g   aL=%.7g\n"
                 "wl=%.10g  wn=%.10g\n"
                 "k =%12.5g   // densiso[imol]\n"
                 "  *%12.5g   // isoratio\n"
                 "  *%12.5g   // SIGCTE\n"
                 "  *%12.5g   // lt->gf[ln]\n"
                 "  *%12.5g   // exp(-EXPCTE*lt->elow[ln]/temp)\n"
                 "  *%12.5g   // (1-exp(-EXPCTE*wavn/temp))\n"
                 "  /%12.5g   // mass[i]\n"
                 "  /%12.5g   // ziso[i][t]\n"
                 " = %12.5g   // extinction\n\n",
                 i, temp, lt->elow[ln],
                 alphad[i]*wavn, alphal[i],
                 lt->wl[ln], 1.0 / (lt->wfct * lt->wl[ln] * tr->wavs.fct),
                 mol->molec[iso->imol[i]].d[r],
                 iso->isoratio[i],
                 SIGCTE,
                 lt->gf[ln],
                 exp(-EXPCTE*lt->elow[ln]/temp),
                 (1-exp(-EXPCTE*wavn/temp)),
                 iso->isof[i].m,
                 op->ziso[i][t],
                 opac);

    /* profsize is (number of points in the profile)/2.                     */

    /* profsize-ilinewn-1 is the starting point of the wavenumber array with
       respect to the starting index of the profile.                        */
    /* Set profwn such that the index mimic wavenumber's array:             */

    //transitprint(1, 2, "FLAG 279: idop=%d,  ilor=%d\n", idop[i], ilor[i]);
    //transitprint(1, 2, "FLAG 280: profsize=%li\n", profsize[idop[i]][ilor[i]]);

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
      //transitprint(1, 2, "%li  ", j-offset);
      //transitprint(1, 2, "j=%d, p[%li]=%.2g   ", j, j-offset,
      //                    profile[idop[i]][ilor[i]][subw][j-offset]);
      op->o[m][t][r][j] += opac * profile[idop[i]][ilor[i]][subw][j-offset];
    }
    ncalc++;
  }

  t0 = timecheck(verblevel, t, r, "End opacity calculation", tv, t0);
  transitprint(1, verblevel, "There were %d profile evaluations,\n"
               "       and %li line profiles calculated.\n", neval, ncalc);

  transitprint(2, verblevel, "Done.\n");
  return 0;
}

/* \fcnfh
   Free opacity array

   Return: 0 on success                                                     */
int
freemem_opacity(struct opacity *op, /* Opacity structure                    */
                long *pi){          /* transit progress flag                */
  /* Free arrays:                                                           */
  free(op->o[0][0][0]); /* The opacity                                      */
  free(op->o[0][0]);
  free(op->o[0]);
  free(op->o);

  free(op->profile[0][0][0]); /* The Voigt profiles                         */
  free(op->profile[0][0]);
  free(op->profile[0]);
  free(op->profile);

  free(op->profsize[0]);  /* The Voigt-profile half-size                    */
  free(op->profsize);

  //free();

  /* Update progress indicator and return:                                  */
  *pi &= ~(TRPI_OPACITY | TRPI_TAU);
  return 0;
}

