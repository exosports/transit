// Copyright (C) 2015-2016 University of Central Florida. All rights reserved.
// Transit is under an open-source, reproducible-research license (see LICENSE).

#include <readatm.h>
#define ROUNDTHRESH 1e-5
static double zerorad=0;

char *atmfilename;

/* \fcnfh
   Initialize ds.at (atm_data).  Set abundance mass and allowrq parameters.
   Check existence, open, and set pointer to atmosphere file.
   Get keyword variables from atm file (list of isotopes among others).
   Get temperature and isotopes abundances per radius from atm file.

   Return:
     0 on succes, elses
    -1 if no atmospheric info file was specified and no defaults are allowed
    -2 if default handling mode does not exist
    -3 if something really bad happened
    -4 if sum of abundances add up to more than 1
    -5 if requested isotope is an ignored one                               */
int
getatm(struct transit *tr){
  int nmol,                           /* Number of molecules                */
      nrad,                           /* Number of radius samples           */
      i;                              /* for-loop index                     */
  struct transithint *th = tr->ds.th; /* Get transithint                    */
  static struct atm_data at;          /* transit's atm_data structure       */
  static struct molecules mol;        /* transit's molecules struct         */
  prop_samp *rads = &at.rads;         /* Radius sample                      */

  /* Set atm_data and molecules struct's mem to 0:                          */
  memset(&at,  0, sizeof(struct atm_data));
  memset(&mol, 0, sizeof(struct molecules));
  tr->ds.at  = &at;
  tr->ds.mol = &mol;

  FILE *fp = NULL;        /* Pointer to atmospheric file                    */

  at.mass     = th->mass;  /* bool: abundance by mass (1) or by number (0)  */
  tr->allowrq = th->allowrq;

  /* Check if atmospheric file was specified:                               */
  if(th->f_atm == NULL  ||  strcmp(th->f_atm, "-") == 0){
    tr_output(TOUT_ERROR,
      "getatm() :: No atmospheric file specified.\n");
    return -1;
  }
  else{
    /* Check that the file exists and can be opened.  Set name and pointer
       in transit structure:                                                */
    if((tr->fp_atm=verbfileopen(th->f_atm, "Atmospheric info ")) == NULL)
      exit(EXIT_FAILURE);
    atmfilename = tr->f_atm = th->f_atm;   /* Set file name                 */
    fp  = tr->fp_atm;        /* Pointer to file                             */
    tr_output(TOUT_INFO, "Reading atmosphere file: '%s'.\n", tr->f_atm);

    /* nrad will be the amount of allocated radii so far:                   */
    nrad      = 8;
    rads->n   = 0;
    rads->fct = 1.0;
  }

  /* Initialize atmosphere temperature-pressure arrays:                     */
  at.atm.tfct = 1;  /* Default temperature units are cgs                    */
  at.atm.pfct = 1;  /* Default pressure    units are cgs                    */
  rads->v  = (PREC_ATM *)calloc(nrad, sizeof(PREC_ATM));
  at.atm.t = (PREC_ATM *)calloc(nrad, sizeof(PREC_ATM));
  at.atm.p = (PREC_ATM *)calloc(nrad, sizeof(PREC_ATM));
  rads->v[0] = 1.0;

  /* Read keyword-variables from file:                                      */
  if((i=getmnfromfile(fp, &at, tr))<1){
    tr_output(TOUT_ERROR, "getmnfromfile() returned error code %i\n", i);
    exit(EXIT_FAILURE);
  }
  nmol = at.n_aiso; /* Number of molecules in atmospheric file              */

  /* Allocate mass and radius array of the molecules:                       */
  mol.nmol   = nmol;
  mol.ID     = (int       *)calloc(nmol, sizeof(int));
  mol.mass   = (PREC_ZREC *)calloc(nmol, sizeof(PREC_ZREC));
  mol.radius = (PREC_ZREC *)calloc(nmol, sizeof(PREC_ZREC));
  mol.pol    = (PREC_ZREC *)calloc(nmol, sizeof(PREC_ZREC));
  mol.molec  = (prop_mol  *)calloc(nmol, sizeof(prop_mol));

  /* Get (pseudo-fixed) molecular data values from 'molecules.dat' file:    */
  getmoldata(&at, &mol, tr->f_molfile);

  /* Allocate arrays for the mean molecular mass, density, and abundance:   */
  at.molec        = (prop_mol *)calloc(nmol, sizeof(prop_mol));
  at.mm           = (double   *)calloc(nrad, sizeof(double));
  for(i=0; i<nmol; i++){
    at.molec[i].d = (PREC_ATM *)calloc(nrad, sizeof(PREC_ATM));
    at.molec[i].q = (PREC_ATM *)calloc(nrad, sizeof(PREC_ATM));
    at.molec[i].n = nrad;
  }

  /* Read isotopic abundances:                                              */
  nrad = readatmfile(fp, tr, &at, rads, nrad);
  tr_output(TOUT_INFO, "Done.\n\n");
  fclose(fp);

  /* Set required values in 'rads' structure:                               */
  rads->i = rads->v[0];
  rads->f = rads->v[rads->n-1];
  rads->o = 1;
  rads->d = 0;

  /* Return succes and set progress indicator:                              */
  tr->pi |= TRPI_GETATM;
  return 0;
}


/* \fcnfh
   Compute the mean molecular mass, check that sum of abundances
   is no larger than 1.

   Return: Sum of abundances.                                               */
double
checkaddmm(double *mm,            /* Mean molecular mass stored             */
           PREC_NREC r,           /* Radius position                        */
           prop_mol *molec,       /* Molecule info                          */
           struct molecules *mol, /* Molecules                              */
           int n,                 /* Number of molecules                    */
           _Bool mass){           /* Mass abundances?                       */
  double sumq;  /* Fractional abundance sum                                 */
  int i;        /* for-loop index                                           */

  if(r >= molec[0].n) {
    tr_output(TOUT_ERROR,
      "In file %s (line %li) a radius beyond the allocated "
      "has been requested.", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  }

  /* Compute the mean molecular mass:                                       */
  sumq = *mm = 0;
  for(i=0; i<n; i++){
    if(mass)
      *mm += (molec[i].q[r])/(mol->mass[i]);
    else
      *mm += (molec[i].q[r])*(mol->mass[i]);
    sumq += molec[i].q[r];
  }
  if(mass)
    *mm = 1.0/(*mm);

  /* Check that sum of proportional abundances make sense:                  */
  if(sumq>1.001){
    tr_output(TOUT_ERROR, "Sum of abundances of isotopes adds up to "
                               "more than 1: %g\n", sumq);
    exit(EXIT_FAILURE);
  }

  return sumq;
}


/* \fcnfh
   Print out default values when only one radius is being selected.         */
void
telldefaults(struct isotopes *iso,
             struct atm_data *at){
  int i;
  tr_output(TOUT_INFO,
    "You are using one point atmospheric conditions:\n"
    " Temperature:         %g K\n"
    " Pressure:            %g dyne/cm2\n"
    " Mean molecular mass: %g AMU\n",
    at->atm.t[0]*at->atm.tfct, at->atm.p[0]*at->atm.pfct, at->mm[0]);
  /* Densities for all isotopes: */
  for(i=0; i<iso->n_i; i++)
      tr_output(TOUT_INFO, " %-8s: density %8g g/cm3\n",
        iso->isof[i].n, at->molec[iso->imol[i]].d[0]);
}


/* \fcnfh
   Free memory from atmosphere structure
   Return: 0 on success                                                     */
int
freemem_atmosphere(struct atm_data *at,
                   long *pi){
  /* Free structures:                                                       */
  free_samp(&at->rads);
  for(int i=0; i<at->n_aiso; i++)
    free_mol(at->molec+i);
  free_atm(&at->atm);

  /* Free arrays:                                                           */
  free(at->mm);
  free(at->info);

  /* Clear progress indicator and return success:                           */
  *pi &= ~(TRPI_GETATM);
  return 0;
}


/* \fcnfh
   Store info about the atmopshere file                                     */
void
storename(struct atm_data *at,
          char *line){
  /* Get the length of the line char array:                                 */
  int len;
  while(*line==' ' || *line=='\t')
    line++;
  len = strlen(line);

  /* Allocate and store info in atm_data.info:                              */
  if(!at->info){  /* Only if it hasn't been stored before:                  */
    at->info = calloc(len+1, sizeof(char));
    strcpy(at->info, line);
  }
}


/* \fcnfh
   Print error message when a line of file is longer than max characters.   */
static void
atmerr(int max,     /* Maximum length of an accepted line                   */
       char *file,  /* File from which we were reading                      */
       int line){   /* Line being read                                      */
  tr_output(TOUT_ERROR,
    "Line %i of file '%s' has more than %i characters, "
    "that is not allowed\n", file, max);
  exit(EXIT_FAILURE);
}


/* \fcnfh
   Print error message when a field with transition info is invalid.        */
static void
invalidfield(char *line,   /* Contents of the line                          */
             int nmb,      /* File number                                   */
             int fld,      /* Field with the error                          */
             char *fldn){  /* Name of the field                             */
  tr_output(TOUT_ERROR,
    "Line %i of file '%s': Field %i (%s) does not have a valid "
    "value:\n%s.\n", nmb, atmfilename, fld, fldn, line);
  exit(EXIT_FAILURE);
}


/* \fcnfh
   Check that val is positive. Throw error message if not.                  */
static void
checkposvalue(PREC_RES val, /* Value to check                               */
              int field,    /* Field where it was read                      */
              long line){   /* Line from which it was read                  */

  if(val<0) {
    tr_output(TOUT_ERROR,
      "While reading the %i-th field in line %li of atmosphere "
      "file %s, a negative value was found (%g).\n", field,
      line-1, atmfilename, val);
    exit(EXIT_FAILURE);
  }
}


/* \fcnfh
    Get keyword variables from atmosphere file (mass/number abundance bool;
    zero-radius offset; radius, temperature, and pressure units factor;
    atmfile name/info; list isotopes; list of proportional-abundance isotopes).
    Store molecules and proportional isotopes in atm_data struct.
    Determine which linedb isotope corresponds to such atm_data isotope.
    Solve non-matched linedb isotope cases.
    Put all non-ignore isotopes in transit.ds.iso structure.

    Return: Number of lines read                                            */
int
getmnfromfile(FILE *fp,                /* Pointer to atmospheric file       */
              struct atm_data *at,     /* atmosphere structure              */
              struct transit *tr){     /* Remainder molecules' factor       */
  struct molecules *mol=tr->ds.mol;
  char line[maxline], *lp, keyword[maxline];
  int nmol=0,            /* Number of species in atmospheric file           */
      i, j;              /* Auxiliary for-loop index                        */

  at->begline = 0; /* Line where the info begins                            */

  /* Read and store the keyword atmospheric variables:                      */
  while(1){
    switch(fgetupto_err(line, maxline, fp, &atmerr, tr->f_atm,
                        at->begline++)){
    /* Ignore comments and blank lines:                                     */
    case '\n':
      continue;
    case '#':
      getname(line+1, keyword);
      /* Molecule names with an abundance profile:                          */
      if (strcmp(keyword, "SPECIES") == 0){
        /* Go to next line:                                                 */
        fgetupto_err(line, maxline, fp, &atmerr, tr->f_atm, at->begline++);
        /* Count the number of words (molecules) in line:                   */
        nmol = countfields(line, ' ');
        tr_output(TOUT_RESULT, "The number of molecules is %d.\n", nmol);

        /* Allocate Molecules names:                                        */
        mol->name    = (char **)calloc(nmol,             sizeof(char *));
        mol->name[0] = (char  *)calloc(nmol*maxeisoname, sizeof(char));
        for(i=1; i<nmol; i++)
          mol->name[i] = mol->name[0] + i*maxeisoname;
        /* Read and store the molecule names:                               */
        tr_output(TOUT_DEBUG, "Molecules with abundance profile:\n  ");
        lp = line;
        /* Read and store names:                                            */
        for (i=0; i<nmol; i++){
          getname(lp, mol->name[i]);
          lp = nextfield(lp);
          if (i < nmol-1)
            tr_output(TOUT_DEBUG, "%s, ", mol->name[i]);
          else
            tr_output(TOUT_DEBUG, "%s.\n", mol->name[i]);
        }
      }
      continue;
    case 0:     /* Throw error if EOF:                                      */
      tr_output(TOUT_ERROR,
        "readatm :: EOF unexpectedly found at line %i "
        "of file %s while no t,p data points have been read.\n",
        at->begline, tr->f_atm);
      exit(EXIT_FAILURE);
      continue;

    /* Determine whether abundance is by mass or number:                    */
    case 'q':
      lp = line + 1;
      while(*lp++ == ' '); /* Skip blank spaces                             */
      lp--;
      switch(*lp|0x20){
      case 'n':
        at->mass = 0;      /* Number abundance (mixing ratio)               */
        break;
      case 'm':
        at->mass = 1;      /* Mass abundance (mass mixing ratio)            */
        break;
      default:
        tr_output(TOUT_ERROR,
          "'q' option in the atmosphere file can only be followed "
          "by 'm' (for abundances by mass) or 'n' (for abundances "
          "by number). '%s' is invalid.\n", line);
        break;
      }
      continue;

    /* Zero radius value:                                                   */
    case 'z':
      zerorad = atof(line+1);
      continue;

    /* Radius, temperature, or pressure units factor:                       */
    case 'u':
      switch(line[1]){
      case 'r':
        at->rads.fct = atof(line+2);
        break;
      case 'p':
        at->atm.pfct = atof(line+2);
        break;
      case 't':
        at->atm.tfct = atof(line+2);
        break;
      default:
        tr_output(TOUT_ERROR,
          "Invalid unit factor indication in atmosphere file.\n");
        exit(EXIT_FAILURE);
      }
      continue;

    case 'n':  /* Name or identifier for file data:                         */
      storename(at, line+1);
      continue;

    /* End of keyword variables:                                            */
    default:
      break;
    }
    break;
  }

  tr_output(TOUT_RESULT, "Read all keywords in atmosphere file without "
    "problems.\n");

  /* Set total number of molecules in atmosphere:                           */
  mol->nmol = at->n_aiso = nmol;

  /* Apply abundance scale factor if any:                                   */
  if (tr->nqmol > 0){
    tr->qmol = (int *)calloc(tr->nqmol, sizeof(int));
    lp = tr->ds.th->qmol;
    for (i=0; i < tr->nqmol; i++){
      tr->qmol[i] = -1;
      getname(lp, keyword);
      for (j=0; j<nmol; j++)
        if (strcmp(keyword, mol->name[j])== 0){
          tr->qmol[i] = j;
          break;
        }
      lp = nextfield(lp);
    }
  }
  /* Check that there was at least one isotope defined and re-allocate
     array sizes to their final size:                                       */
  if(!nmol) {
    tr_output(TOUT_ERROR,
      "No species were found in the atmospheric file, "
      "make sure to specify them with the comment/"
      "header in the previous line '#SPECIES'.\n");
    exit(EXIT_FAILURE);
  }

  /* Set position of beginning of data:                                     */
  at->begpos = ftell(fp) - strlen(line) - 1;

  /* Resolve what to do with those isotopes that appear in the
     line transition database, but not in the atmosphere file. Get
     the number of non-ignored isotopes in atm_data without linelist:       */
  //at->n_niso = checknonmatch(tr, at, isodo);
  /* FINDME: This will be a task in readline (if actually needed).          */

  return at->begline;
}


/* \fcnfh
    Read radius, pressure, temperature, and abundances and store it into
    at_data of transit.  Calculate mean molecular mass and densities.

    Detailed:
    Read and store radius, pressure, and temperature from file.
    Read abundances for each (non other-factor) isotope.
    Sum fractional abundances. Calculate ramaining (other-factor) abundances.
    Calculate mean molecular mass per radius.
    Calculate densities per isotope at each radius.

    Returns: number of sample radius                                        */
int
readatmfile(FILE *fp,                /* Atmospheric file                    */
            struct transit *tr,      /* transit struct                      */
            struct atm_data *at,     /* Atmosphere struct                   */
            prop_samp *rads,         /* Radius sampling                     */
            int nrad){               /* Remainder molecules' factor         */

  tr_output(TOUT_INFO, "Start reading abundances.\n");
  /* Find abundance related quantities for each radius                      */
  int lines = at->begline;
  PREC_NREC r = 0;  /* Radius index (number of radii being read)            */
  char rc;          /* File reading output                                  */
  float allowq = tr->allowrq;
  int nabundances;  /* Number of abundances in list                         */
  double sumq, /* Sum of abundances per line                         */
         ratio, sumq2=0.0;
  char line[maxline], *lp, *lp2;
  prop_mol *molec = at->molec;
  struct molecules *mol = tr->ds.mol;
  int sorted=1,    /* Flag to check if the layers are bottom-top sorted     */
      reversed=1;  /* Flag to check if the layers are top-bottom sorted     */
  int i, j,        /* Auxiliary for-loop indices                            */
      iH2=valueinarray(mol->ID, 105, mol->nmol),
      iHe=valueinarray(mol->ID,   2, mol->nmol);

  /* Count the number of abundances in each line:                           */
  fseek(fp, at->begpos, SEEK_SET); /* Go to position where data begins      */
  /* Skip comments:                                                         */
  while((rc=fgetupto_err(lp=line, maxline, fp, &atmerr, tr->f_atm, lines++))
        =='#' || rc=='\n');
  /* Count values per line:                                                 */
  nabundances = countfields(lp, ' ') - 3; /* Subtract rad, p, and T columns */
  fseek(fp, at->begpos, SEEK_SET); /* Go to position where data begins      */
  while(1){
    /* Reallocate if necessary:                                             */
    if(r==nrad){
      nrad <<= 1;
      rads->v     = (PREC_ATM *)realloc(rads->v,   nrad*sizeof(PREC_ATM));
      at->atm.t   = (PREC_ATM *)realloc(at->atm.t, nrad*sizeof(PREC_ATM));
      at->atm.p   = (PREC_ATM *)realloc(at->atm.p, nrad*sizeof(PREC_ATM));
      at->mm      = (double   *)realloc(at->mm,    nrad*sizeof(double));
      for(i=0; i<at->n_aiso; i++){
        molec[i].d = (PREC_ATM *)realloc(molec[i].d, nrad*sizeof(PREC_ATM));
        molec[i].q = (PREC_ATM *)realloc(molec[i].q, nrad*sizeof(PREC_ATM));
        molec[i].n = nrad;
      }
    }

    /* Skip comments and read next line:                                    */
    while((rc=fgetupto_err(lp=line, maxline, fp, &atmerr, tr->f_atm, lines++))
          =='#' || rc=='\n');
    /* If it is end of file, stop loop:                                     */
    if(!rc)
      break;

    /* Read and store radius, pressure, and temperature from file:          */
    rads->v[r] = strtod(lp, &lp2) + zerorad;  /* Radius                     */
    checkposvalue(rads->v[r], 1, lines);      /* Check value is positive    */
    if(lp==lp2)
      invalidfield(line, lines, 1, "radius");
    at->atm.p[r] = strtod(lp2, &lp);          /* Pressure                   */
    checkposvalue(at->atm.p[r], 2, lines);
    if(lp==lp2)
      invalidfield(line, lines, 2, "pressure");
    at->atm.t[r] = strtod(lp, &lp2);          /* Temperature                */
    checkposvalue(at->atm.t[r], 3, lines);
    if(lp==lp2)
      invalidfield(line, lines, 3, "temperature");

    /* Read abundances for each isotope.  Keep reading-in values
       while there are numbers in line:                                     */
    for(i=0, sumq=0; i<nabundances; i++){
      lp = lp2;
      /* Read the abundance of the isotope:                                 */
      molec[i].q[r] = strtod(lp, &lp2);

      if (tr->nqmol > 0){
        if ( (j=valueinarray(tr->qmol, i, tr->nqmol)) >= 0)
          molec[i].q[r] *= pow(10.0, tr->qscale[j]);
      }
      if (r==0)
        tr_output(TOUT_DEBUG, "density[%d, %li]: %.9f.\n",
                                    i, r, molec[i].q[r]);

      sumq += molec[i].q[r];      /* Add the abundances       */
      if (i != iH2 && i != iHe)
        sumq2 += molec[i].q[r];   /* Sum of metals abundance  */

      checkposvalue(molec[i].q[r], i+4, lines); /* Check that tmp is > 0    */
      if(lp==lp2)
        invalidfield(line, lines, 4+i, "isotope abundance");
    }

    if (tr->nqmol >0){
      ratio = molec[iH2].q[r] / molec[iHe].q[r];
      molec[iHe].q[r] =         (1-sumq2) / (1.0 + ratio);  /* He */
      molec[iH2].q[r] = ratio * (1-sumq2) / (1.0 + ratio);  /* H2 */
      sumq2 = 0.0;
    }
    transitASSERT(i!=at->n_aiso, "The line %s of file %s contains %d abundance "
                                 "values, when there were %d expected.\n",
                                 __LINE__, __FILE__, i, at->n_aiso);

    /* Calculate mean molecular mass and check whether abundances add up
       to one (within allowq threshold):                                    */
    sumq = checkaddmm(at->mm+r, r, molec, mol, at->n_aiso, at->mass);
    if(fabs(sumq-1.0) > allowq)
      tr_output(TOUT_WARN,
        "In radius %i (%g km), abundances don't add "
        "up to 1.0: %.9g\n", r, at->rads.v[r], sumq);

    /* Calculate densities using ideal gas law:                             */
    if (r>=0){
      tr_output(TOUT_DEBUG, "Abund: %.9f, mmm: %.3f, mass: %.3f, "
        "p: %.3f, T: %.3f.\n", molec[2].q[r], at->mm[r], mol->mass[2],
        at->atm.p[r]*at->atm.pfct, at->atm.t[r]*at->atm.tfct);
    }
    for(i=0; i<at->n_aiso; i++){
      molec[i].d[r] = stateeqnford(at->mass, molec[i].q[r], at->mm[r],
                                   mol->mass[i], at->atm.p[r]*at->atm.pfct,
                                   at->atm.t[r]*at->atm.tfct);
    }
    r++;
  }

  /* Re-allocate arrays to final size (nrad):                               */
  rads->n = nrad = r;
  rads->v   = (PREC_ATM *)realloc(rads->v,   nrad*sizeof(PREC_ATM));
  at->atm.t = (PREC_ATM *)realloc(at->atm.t, nrad*sizeof(PREC_ATM));
  at->atm.p = (PREC_ATM *)realloc(at->atm.p, nrad*sizeof(PREC_ATM));
  at->mm    = (double   *)realloc(at->mm,    nrad*sizeof(double));
  for(i=0; i<at->n_aiso; i++){
    molec[i].d = (PREC_ATM *)realloc(molec[i].d, nrad*sizeof(PREC_ATM));
    molec[i].q = (PREC_ATM *)realloc(molec[i].q, nrad*sizeof(PREC_ATM));
    molec[i].n = nrad;
  }
  tr_output(TOUT_RESULT, "Read %li layers between pressures %.3e and "
    "%.3e barye.\n", r, at->atm.p[0  ]*at->atm.pfct,
    at->atm.p[r-1]*at->atm.pfct);

  /* Check that the atmospheric layers are sorted from the bottom to
     the top layer:                                                         */
  for (i=0; i<nrad-1; i++){
    /* Sorted atmospheric layers:                                           */
    if ((rads->v[i] >= rads->v[i+1]) || (at->atm.p[i] <= at->atm.p[i+1])){
      sorted = 0;
    }
    /* Reversed atmospheric layers:                                         */
    if ((rads->v[i] <= rads->v[i+1]) || (at->atm.p[i] >= at->atm.p[i+1])){
      reversed = 0;
    }
  }

  if (sorted == 0 && reversed == 0) {
    tr_output(TOUT_ERROR,
      "The atmospheric layers are neither sorted "
      "from the bottom up, not from the top down.\n");
    exit(EXIT_FAILURE);
  }
  else if (reversed == 1){
    tr_output(TOUT_WARN,
      "The atmospheric layers are in reversed order (top-bottom). "
      "Resorting to be from the bottom-up.\n");
    /* Swap atmospheric-layer values:                                       */
    for (i=0; i<nrad/2; i++){
      swap(rads->v,   i, nrad-1-i);
      swap(at->atm.p, i, nrad-1-i);
      swap(at->atm.t, i, nrad-1-i);
      swap(at->mm,    i, nrad-1-i);
      for (j=0; j<at->n_aiso; j++){
        swap(molec[j].d, i, nrad-1-i);
        swap(molec[j].q, i, nrad-1-i);
      }
    }
  }

  return nrad;
}


/* Read and store non-layer-dependent molecular data (mass, radius, ID, 
   polarizability) and store in mol struct.                                 */
void
getmoldata(struct atm_data *at, struct molecules *mol, char *filename){
  int nmol = at->n_aiso;
  FILE *elist;

  double *mmass,       /* Molecular mass from list                          */
         *radius,      /* Molecular radii from list                         */
         *pol;         /* Molecular polarizability from list                */
  int *molID;          /* Molecular universal ID                            */
  long dinit;          /* Data initial position in File pointer             */
  int ndatamol=0,      /* Number of species in file                         */
      maxlinelen=501,  /* Line maximum length                               */
      i, j;            /* Auxiliary for-loop index                          */
  char **rname,        /* Molecules names for listed radii                  */
       line[maxlinelen], *lp;


  /* Open Molecular data file:                                              */
  if((elist=verbfileopen(filename, "Molecular info ")) == NULL)
    exit(EXIT_FAILURE);

  /* Read lines, skipping comments and blank lines:                         */
  do{
    dinit = ftell(elist);
    lp = fgets(line, maxlinelen, elist);
  }while (lp[0] == '\0' || lp[0] == '\n' || lp[0] == '#');

  /* Count the number of species:                                           */
  do{
    lp = fgets(line, maxlinelen, elist);
    ndatamol += 1;
    if (lp == NULL)
      break;
  }while (lp[0] != '\0' && lp[0] != '\n' && lp[0] != '#');

  /* Re-locate pointer to the data initial position and read first line:    */
  fseek(elist, dinit, SEEK_SET);
  lp = fgets(line, maxlinelen, elist);

  /* Allocate molecule ID, mass, names, polarizability, and radii:          */
  molID    = (int    *)calloc(ndatamol,            sizeof(int));
  mmass    = (double *)calloc(ndatamol,            sizeof(double));
  radius   = (double *)calloc(ndatamol,            sizeof(double));
  pol      = (double *)calloc(ndatamol,            sizeof(double));
  rname    = (char  **)calloc(ndatamol,            sizeof(char *));
  rname[0] = (char   *)calloc(ndatamol*MAXNAMELEN, sizeof(char));
  for (i=1; i<ndatamol; i++)
    rname[i] = rname[0] + i*MAXNAMELEN;

  /* Read molecular info from file:                                         */
  for (i=0; i<ndatamol; i++){
    /* Get universal molecule's ID:                                         */
    molID[i] = (int)strtol(lp, NULL, 10);
    /* Get molecule's name:                                                 */
    lp = nextfield(lp);
    getname(lp, rname[i]);
    /* Get molecule's mass:                                                 */
    lp = nextfield(lp);
    mmass[i] = strtod(lp, NULL);
    /* Get molecule's radius:                                               */
    lp = nextfield(lp);
    radius[i] = strtod(lp, NULL)/2.0;
    /* Get molecule's polarizability:                                       */
    lp = nextfield(lp); // Skip over the radius source
    lp = nextfield(lp);
    pol[i] = strtod(lp, NULL);
    lp = fgets(line, maxlinelen, elist);
    tr_output(TOUT_DEBUG, "Read species:'%6s',  ID:%3d,  mass:%7.4f,  "
      "radius:%5.2f.\n", rname[i], molID[i], mmass[i], radius[i]);
  }

  /* Assign info for each molecule in atosphere:                            */
  for (i=0; i<nmol; i++){
    /* Set the radius:                                                      */
    j = findstring(mol->name[i], rname, ndatamol);
    if (j == -1) {
      tr_output(TOUT_ERROR,
        "The atmospheric species '%s' is not present in the list "
        "of known species:\n '%s'.\n", mol->name[i], filename);
      exit(EXIT_FAILURE);
    }
    mol->radius[i] = radius[j] * ANGSTROM;
    /* Set the universal molecular ID:                                      */
    mol->ID[i]     = molID[j];
    /* Set the mass:                                                        */
    mol->mass[i]   = mmass[j];
    /* Set the polarizability:                                              */
    mol->pol[i]    = pol[j];
    tr_output(TOUT_RESULT, "Molecule %9s (ID: %3i) has radius %4.2f A "
      "and mass %7.4f u.\n", mol->name[i], mol->ID[i],
      mol->radius[i]/ANGSTROM, mol->mass[i]);
  }
}


/*
Re-load data from array into transit's atm structure.                       */
int
reloadatm(struct transit *tr,
          double *input){ /* Array with updated temp. and abund. profiles   */

  struct atm_data   *at=tr->ds.at;   /* Atmosphere data struct              */
  struct molecules *mol=tr->ds.mol;  /* Molecules stucture                  */
  int nlayers = at->rads.n,          /* Number of layers                    */
      nmol    = mol->nmol;           /* Number of molecules                 */
  int i, j;                          /* Auxilliary for-loop indices         */
  double sumq,
	 allowq = tr->allowrq;

  /* Update atmfile temperature array:                                      */
  for (i=0; i<nlayers; i++){
    at->atm.t[i] = input[i];
    tr_output(TOUT_DEBUG, "%8.3f  ", at->atm.t[i]);
  }
  tr_output(TOUT_DEBUG, "\n");

  /* Update atmfile abundance arrays:                                       */
  for (j=0; j<nmol; j++)
    for (i=0; i<nlayers; i++)
      at->molec[j].q[i] = input[nlayers*(j+1) + i];

  /* Re-calculate the mean molecular mass and check whether abundances add up
     to one (within allowq threshold):                                    */
  for (i=0; i<nlayers; i++){
    sumq = checkaddmm(at->mm+i, i, at->molec, mol, nmol, at->mass);
    if(fabs(sumq-1.0) > allowq)
      tr_output(TOUT_WARN,
        "In radius %i (%g km), abundances don't add "
        "up to 1.0: %.9g\n", i, at->rads.v[i], sumq);

    /* Re-calculate densities:                                              */
    for(j=0; j<nmol; j++)
      at->molec[j].d[i] = stateeqnford(at->mass, at->molec[j].q[i], at->mm[i],
                                       mol->mass[j], at->atm.p[i]*at->atm.pfct,
                                                     at->atm.t[i]*at->atm.tfct);
  }

  /* Re-calculate radius:                                                   */
  tr_output(TOUT_DEBUG, "Old radius boundaries: [%.1f, %.1f]\n",
    at->rads.v[0], at->rads.v[nlayers-1]);
  /* Check that r0, p0, and gsurf were defined:                             */
  if (tr->p0 == 0 || tr->r0 == 0 || tr->gsurf == 0) {
    tr_output(TOUT_ERROR,
      "Surface gravity (%.1f) or reference pressure "
      "(%.3e) or radius (%.1f) were not defined.\n",
      tr->gsurf, tr->p0, tr->r0);
    exit(EXIT_FAILURE);
  }
  radpress(tr->gsurf, tr->p0, tr->r0, at->atm.t,
           at->mm, at->atm.p, at->rads.v, nlayers, at->rads.fct);
  /* Actualize other variables of at->rads:                                 */
  at->rads.i = at->rads.v[0];
  at->rads.f = at->rads.v[nlayers-1];
  tr_output(TOUT_DEBUG, "New radius boundaries: [%.1f, %.1f]\n",
    tr->ds.at->rads.v[0], tr->ds.at->rads.v[nlayers-1]);

  /* Re-resample transit arrays:                                            */
  makeradsample(tr);
  return 0;
}


int radpress(double g0,        /* Surface gravity (cm/s^2)                  */
             double p0,        /* Reference pressure                        */
             double r0,        /* Reference height                          */
             double *temp,     /* Temperature (K)                           */
             double *mu,       /* Mean molecular mass (g/mol)               */
             double *pressure, /* Pressure array                            */
             double *radius,   /* Radius array                              */
             int nlayers,      /* Number of layers                          */
             double rfct){     /* Radius-array units                        */
  double temp0, mu0; /* Temperature, mean moecular mass at p0               */
  double g;       /* Gravity at a given layer                               */
  int i, i0=-1;   /* Indices                                                */
  double minimum=1E+37;   /* Running minimum, to determine i0               */

  /* Find pressure index closest to p0                                      */
  for (i=0; i<nlayers; i++){
    if (fabs(pressure[i]-p0) < minimum){
      i0 = i;
      minimum = fabs(pressure[i]-p0);
    }
  }
  /* Temp, mean mol. mass, radius & gravity at this layer                   */
  if (pressure[i0] > p0){
    temp0 = temp[i0] + ((temp[i0+1]-temp[i0])/log(pressure[i0+1]/pressure[i0])) * 
                       log(p0/pressure[i0]);
    mu0   = mu[i0]   + ((mu[i0+1]  -mu[i0]  )/log(pressure[i0+1]/pressure[i0])) * 
                       log(p0/pressure[i0]);
    radius[i0] = r0 + 0.5 * (temp[i0] / mu[i0] + temp0 / mu0) * 
                      (KB/AMU * log(p0/pressure[i0]) / g0) / rfct;
  }
  else{
    temp0 = temp[i0] + ((temp[i0-1]-temp[i0])/log(pressure[i0-1]/pressure[i0])) * 
                       log(p0/pressure[i0]);
    mu0   = mu[i0]   + ((mu[i0-1]  -mu[i0]  )/log(pressure[i0-1]/pressure[i0])) * 
                       log(p0/pressure[i0]);
    radius[i0] = r0 - 0.5 * (temp[i0] / mu[i0] + temp0 / mu0) * 
                      (KB/AMU * log(pressure[i0]/p0) / g0) / rfct;
  }
  g = g0 * pow(r0 / radius[i0], 2);

  /* Reference pressure is not in given range:                              */
  if (i0 == -1){
    tr_output(TOUT_ERROR,
      "Reference pressure level (%.3e) not found "
      "in range [%.3e, %.3e].\n", p0, pressure[0], pressure[nlayers-1]);
    exit(EXIT_FAILURE);
    return 0;
  }

  /* Calculate radii below p0                                               */
  for (i=i0-1; i>=0; i--){
    radius[i] = radius[i+1] - 0.5 * (temp[i] / mu[i] + temp[i+1] / mu[i+1]) * 
                              (KB/AMU * log(pressure[i]/pressure[i+1]) / g) / 
                              rfct;
    g = g * pow(radius[i+1] / radius[i], 2);
  }

  g = g0 * pow(r0 / radius[i0], 2);

  /* Calculate radii above p0                                               */
  for (i=i0+1; i<nlayers; i++){
    radius[i] = radius[i-1] + 0.5 * (temp[i] / mu[i] + temp[i-1] / mu[i-1]) * 
                              (KB/AMU * log(pressure[i-1]/pressure[i]) / g) / 
                              rfct;
    g = g * pow(radius[i-1] / radius[i], 2);
  }

  tr_output(TOUT_DEBUG, "PRESSURE:\n");
  for (i=0; i<nlayers; i++)
    tr_output(TOUT_DEBUG, "%.3e, ", pressure[i]);
  tr_output(TOUT_DEBUG, "\n");

  tr_output(TOUT_DEBUG, "RADIUS:\n");
  for (i=0; i<nlayers; i++){
    tr_output(TOUT_DEBUG, "%.1f, ", radius[i]);
  }
  tr_output(TOUT_DEBUG, "\n");
  return 1;
}
