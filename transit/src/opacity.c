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
  transitprint(10, verblevel, "Opacity-file exist status = %d\n", fe);

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
    /* Create file:                                                         */
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
  //PREC_VOIGT ***vprofile;          /* Grid of Voigt profiles                */
  int   voigtfine =tr->voigtfine;  /* Voigt profile wn oversampling factor  */
  float timesalpha=tr->timesalpha; /* Voigt wings width                     */

  struct timeval tv;  /* Time-keeping variables                             */
  double t0=0.0;

  /* Run transitcheckcalled: TBD                                            */

  /* Make logscale grid for the profile widths:                             */
  /* FINDME: Hardcoded values                                               */
  nDop = op->nDop = 40;  //100;
  nLor = op->nLor = 40;  //100;
  Lmin = 1e-4;
  Lmax = 10.0; //10.0;
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

  /* Re-work of Voigt profile:                                              */
  /* FINDME: This is not really working.  Check why the hell it isn't       */
  //op->vprofile    = (PREC_VOIGT ***)calloc(nDop,      sizeof(PREC_VOIGT **));
  //op->vprofile[0] = (PREC_VOIGT  **)calloc(nDop*nLor, sizeof(PREC_VOIGT  *));
  //for (i=0; i<nDop; i++){
  //  op->vprofile[i] = op->vprofile[0] + i*nLor;
  //}
  //vprofile = op->vprofile;

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
        //vprofile[i][j] = vprofile[i-1][j];
      }
      else{ /* Calculate a new profile for given widths:                    */
        op->profsize[i][j] = newprofile(profile[i][j], voigtfine,
                             tr->wns.d/tr->owns.o, op->aDop[i], op->aLor[j],
                             timesalpha, tr->owns.n);
        //op->profsize[i][j] = getprofile(vprofile[i][j], tr->wns.d/tr->owns.o,
        //                    op->aDop[i], op->aLor[j], timesalpha, tr->owns.n);
      }
      transitprint(5, verblevel, "Profile[%2d][%2d] size = %4li  (D=%.3g, "
                                  "L=%.3g).\n", i, j, 2*op->profsize[i][j]+1,
                                   op->aDop[i], op->aLor[j]);
      /* Calculate the integrated area of the lines:                        */
      /* FINDME: I should probably normalize to have the integral = 1.0     */
      // double area;
      //transitprint(1,2, "%.5e, %.5e\n", op->vprofile[i][j][0], op->vprofile[i][j][1]);

      /* Print profile to screen:                                           */
      //if (i == 0 && j <1 ){
      //  for (k=0; k<2*op->profsize[i][j]+1; k++){
      //    transitprint(1,2, "%.5e, ", profile[i][j][0][k]);
      //  }
      //  transitprint(1, 2, "\n");
      //}

      // //transitprint(1, 2, "%.5e  %.5e  %.5e  %.5e  %.5e\n", profile[0][0]);
      // transitprint(1, 2, "Profile Areas: [");
      // for (k=0; k<voigtfine; k++){
      //   area = 0.0;
      //   for (int l=0; l<2*op->profsize[i][j]+1; l++){
      //     area += profile[i][j][k][l];
      //    //transitprint(1,2, "%.5e  ", profile[i][j][k][l]);
      //   }
      //   transitprint(1, 2, "%.8e  ", area*tr->wns.d);
      // }
      // transitprint(1, 2, "]  (%.3e)\n", tr->wns.d);
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
    for (r=0;   r<Nlayer; r++){  /* For each layer:                       */
    //for (r=0;   r<Nlayer; r++){  /* For each layer:                         */
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
  //for     (i=0; i < op->Nmol;   i++){
  //  op->o[i] = op->o[0] + op->Ntemp*i;
  //  for   (t=0; t < op->Ntemp;  t++){
  //    op->o[i][t] = op->o[0][0] + op->Nlayer*(t + op->Ntemp*i);
  //    for (r=0; r < op->Nlayer; r++)
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
  //struct extinction *ex=tr->ds.ex;

  /* Voigt profile variables:                                               */
  PREC_VOIGT ****profile =op->profile;   /* Voigt profile                   */
  PREC_NREC    **profsize=op->profsize;  /* Voigt-profile half-size         */
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
      iown, idwn;

  PREC_NREC onwn, dnwn;

  long minj, maxj, /* Indices of min and max wn array of line profile       */
       offset,      /* Index offset between profile and wn arrays           */
       ncalc;
  PREC_ATM moldensity,          /* Molecular density                        */
           temp;                /* Layer temperature                        */
  double *kmax, *kmin;
  PREC_VOIGTP *alphal, *alphad; /* Lorentz and Doppler widths               */

  PREC_RES dwn, odwn, ddwn, /* Wavenumber sampling interval */
           wavn, next_wn;  
  PREC_NREC nlines,   /* Number of line transitions                         */
            ln,       /* line-transition for-loop index                     */
            subw;     /* Fine binning index for profile of line transition  */

  struct timeval tv;
  double t0 = timestart(tv, "\nOpacity calculation:");
  //transitprint(1, verblevel, "FLAG 210.\n");

  /* per-spectrum minimum and maximum line width:                           */
  double minwidth=1e5, maxwidth=0.0;
 
  /* Temporal extinction array (per molecule):                              */
  double **ktmp;
  int ofactor;  /* Dynamic oversampling factor                              */

  long nadd =0, /* Number of co-added lines                                 */
       nskip=0, /* Number of skipped lines                                  */
       neval=0; /* Number of evaluated profiles                             */


  /* Set up of variables:                                                   */
  nlines = tr->ds.li->n_l; /* Number of line transitions                    */

  niso = iso->n_i;     /* Number of isotopes                                */
  nmol = mol->nmol;    /* Number of molecules in atmosphere                 */

  /* Number of wavenumber samples:                                          */
  onwn = tr->owns.n; /* Oversampled array                                   */

  /* Wavenumber sampling intervals:                                         */
  dwn  = tr->wns.d /tr->wns.o;  /* Wavenumber array                         */
  odwn = tr->owns.d/tr->owns.o; /* Oversampling array                       */

  temp = op->temp[t];  /* Layer temperature                                 */

  /* Constant factors for line widths:                                      */
  fdoppler = sqrt(2*KB*temp/AMU) * SQRTLN2 / LS;
  florentz = sqrt(2*KB*temp/PI/AMU) / (AMU*LS);

  /* Allocate alpha Lorentz and Doppler arrays:                             */
  alphal = (PREC_VOIGTP *)calloc(niso, sizeof(PREC_VOIGTP));
  alphad = (PREC_VOIGTP *)calloc(niso, sizeof(PREC_VOIGTP));

  /* Allocate width indices array:                                          */
  idop = (int *)calloc(niso, sizeof(int));
  ilor = (int *)calloc(niso, sizeof(int));

  kmax = (double *)calloc(op->Nmol, sizeof(double));
  kmin = (double *)calloc(op->Nmol, sizeof(double));

  ktmp = (double **)calloc(op->Nmol, sizeof(double *));
  ktmp[0] = (double *)calloc(op->Nmol*tr->owns.n, sizeof(double));
  for (i=1; i<op->Nmol; i++)
    ktmp[i] = ktmp[0] + tr->owns.n * i;

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
        transitprint(200, verblevel, "AlphaL[%d] = %.3e\n", j,
                     moldensity/mol->mass[j] * csdiameter * csdiameter *
                     sqrt(1/iso->isof[i].m + 1/mol->mass[j]));
    }
    alphal[i] *= florentz;

    /* Doppler profile width (divided by central wavenumber):               */
    alphad[i] = fdoppler / sqrt(iso->isof[i].m);

    /* Print Lorentz and Doppler broadening widths:                         */
    if(i <= 0)
      transitprint(1, verblevel, "Lorentz: %.9f cm-1, Doppler: %.9f cm-1 "
        "broadening (T=%d, r=%d).\n", alphal[i], alphad[i]*tr->wns.v[0],
                                      (int)temp, r);
    /* Max between Doppler and Lorentz widths:                              */
    maxwidth = fmax(alphal[i], alphad[i]*tr->wns.v[0]);
    minwidth = fmin(minwidth, maxwidth);

    /* Search for aDop and aLor indices for alphal[i] and alphad[i]:        */
    idop[i] = binsearchapprox(aDop, alphad[i]*tr->wns.v[0], 0, nDop);
    ilor[i] = binsearchapprox(aLor, alphal[i],              0, nLor);
  }

  /* Set oversampling resolution:                                           */
  for (i=1; i < tr->ndivs; i++)
    if (tr->odivs[i]*(dwn/tr->owns.o) >= 0.5 * minwidth){
      break;
    }
  ofactor = tr->odivs[i-1];         /* Dynamic-sampling oversampling factor */
  ddwn    = odwn * ofactor;         /* Dynamic-sampling grid interval       */
  dnwn    = 1 + (onwn-1) / ofactor; /* Number of dynamic-sampling values    */
  transitprint(1, verblevel, "Dynamic-sampling grid interval: %.9f  "
               "(scale factor:%i)\n", ddwn, ofactor);
  transitprint(1, verblevel, "Number of dynamic-sampling values: %li\n",
                                dnwn);

  t0 = timestart(tv, "Begin opacity calculation:");

  /* Compute the opacity, proceed for every line:                           */
  for(ln=0; ln<nlines; ln++){

    /* Wavenumber of line transition:                                       */
    wavn = 1.0 / (lt->wl[ln] * lt->wfct);
    /* Isotope ID of line transition:                                       */
    i = lt->isoid[ln];
    /* Molecule index in opacity grid for this isotope:                     */
    m = valueinarray(op->molID, iso->imol[i], op->Nmol);

    /* If this line falls beyond limits, skip to next line transition:      */
    if ((wavn < tr->wns.i) || (wavn > tr->owns.v[onwn-1]))
      continue;

    /* Calculate the extinction coefficient divided by the molecular
       abundance:                                                           */
    opac = iso->isoratio[i]                        * /* Density             */
           SIGCTE     * lt->gf[ln]                 * /* Constant * gf       */
           exp(-EXPCTE*lt->efct*lt->elow[ln]/temp) * /* Level population    */
           (1-exp(-EXPCTE*wavn/temp))              / /* Induced emission    */
           iso->isof[i].m                          / /* Isotope mass        */
           op->ziso[i][t];

    if (kmax[m] == 0)
      kmax[m] = kmin[m] = opac;
    else{
      kmax[m] = fmax(kmax[m], opac);
      kmin[m] = fmin(kmin[m], opac);
    }
  }

  transitprint(1, verblevel, "Half way through\n");
  /* Compute the opacity for each molecule:                                 */
  for(ln=0; ln<nlines; ln++){
    wavn = 1.0/(lt->wl[ln]*lt->wfct); /* Wavenumber  */
    i    = lt->isoid[ln];             /* Isotope ID  */
    /* Molecule index in opacity grid for this isotope:                     */
    m = valueinarray(op->molID, iso->imol[i], op->Nmol);

    if((wavn < tr->wns.i) || (wavn > tr->owns.v[onwn-1]))
      continue;

    opac = iso->isoratio[i] * SIGCTE * lt->gf[ln]        *
           exp(-EXPCTE * lt->efct * lt->elow[ln] / temp) *
           (1-exp(-EXPCTE*wavn/temp)) / iso->isof[i].m / op->ziso[i][t];

    /* Index of closest oversampled wavenumber:                             */
    iown = (wavn - tr->wns.i)/odwn;
    if (fabs(wavn - tr->owns.v[iown+1]) < fabs(wavn - tr->owns.v[iown]))
      iown++;

    /* Check if the next line falls on the same sampling index:             */
    while (ln != nlines-1 && lt->isoid[ln+1] == i){
      next_wn = 1.0 / (lt->wl[ln+1] * lt->wfct);
      if (fabs(next_wn - tr->owns.v[iown]) < odwn){
        nadd++;
        ln++;
        /* Add the contribution from this line into the opacity:            */
        opac += iso->isoratio[i] * SIGCTE * lt->gf[ln] *
                exp(-EXPCTE * lt->efct * lt->elow[ln] / temp)  *
                (1-exp(-EXPCTE*next_wn/temp)) / iso->isof[i].m /
                iso->isov[i].z[r];
      }
      else
        break;
    }

    //if (opac < ex->ethresh * kmax[imol]){
    if (opac < tr->ds.th->ethresh * kmax[m]){
      nskip++;
      continue;
    }

    /* Index of closest (but not larger than) dynamic-sampling wavenumber:  */
    idwn = (wavn - tr->wns.i)/ddwn;

    /* FINDME: de-hard code this threshold                                  */
    if (alphad[i]*wavn/alphal[i] >= 1e-1){
      /* Recalculate index for Doppler width:                               */
      idop[i] = binsearchapprox(aDop, alphad[i]*wavn, 0, nDop);
      ncalc++;
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


    /* Add the contribution from this line to the opacity spectrum:         */
    for(j=minj; j<maxj; j++){
      //transitprint(1, 2, "%li  ", j-offset);
      //transitprint(1, 2, "j=%d, p[%li]=%.2g   ", j, j-offset,
      //                    profile[idop[i]][ilor[i]][subw][j-offset]);
      ktmp[m][j] += opac * profile[idop[i]][ilor[i]][0][ofactor*j - offset];
      //op->o[m][t][r][j] += opac * profile[idop[i]][ilor[i]][subw][j-offset];
    }
    neval++;
  }

  transitprint(1, verblevel, "Almost there!\n");

  /* Downsample ktmp to the final sampling size:                            */
  for (m=0; m < op->Nmol; m++)
    downsample(ktmp[m], op->o[m][t][r], dnwn, tr->owns.o/ofactor);

  t0 = timecheck(verblevel, t, r, "End opacity calculation", tv, t0);

  transitprint(9, verblevel, "Number of co-added lines:     %8li  (%5.2f%%)\n",
                             nadd,  nadd *100.0/nlines);
  transitprint(9, verblevel, "Number of skipped profiles:   %8li  (%5.2f%%)\n",
                             nskip, nskip*100.0/nlines);
  transitprint(9, verblevel, "Number of evaluated profiles: %8li  (%5.2f%%)\n",
                             neval, neval*100.0/nlines);

  /* Free allocated memory:                                                 */
  free(alphal);
  free(alphad);
  free(idop);
  free(ilor);
  free(kmax);
  free(kmin);
  free(ktmp[0]);
  free(ktmp);

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

