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
    transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
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
    transitprint(1, verblevel, "Reading atmosphere file: '%s'.\n", tr->f_atm);

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
    transiterror(TERR_SERIOUS, "getmnfromfile() returned error code %i\n", i);
    exit(EXIT_FAILURE);
  }
  nmol = at.n_aiso; /* Number of molecules in atmospheric file              */

  /* Allocate mass and radius array of the molecules:                       */
  mol.nmol   = nmol;
  mol.ID     = (int       *)calloc(nmol, sizeof(int));
  mol.mass   = (PREC_ZREC *)calloc(nmol, sizeof(PREC_ZREC));
  mol.radius = (PREC_ZREC *)calloc(nmol, sizeof(PREC_ZREC));
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
  transitprint(1, verblevel, "Done.\n\n");
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

  if(r >= molec[0].n)
    transiterror(TERR_CRITICAL,
                 "In file %s (line %li) a radius beyond the allocated "
                 "has been requested.", __FILE__, __LINE__);

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
    transiterror(TERR_SERIOUS, "Sum of abundances of isotopes adds up to "
                               "more than 1: %g\n", sumq);
  }

  return sumq;
}


/* \fcnfh
   Print out default values when only one radius is being selected.         */
void
telldefaults(struct isotopes *iso,
             struct atm_data *at){
  int i;
  transitprint(1, verblevel,
               "You are using one point atmospheric conditions:\n"
               " Temperature:         %g K\n"
               " Pressure:            %g dyne/cm2\n"
               " Mean molecular mass: %g AMU\n",
               at->atm.t[0]*at->atm.tfct, at->atm.p[0]*at->atm.pfct, at->mm[0]);
  /* Densities for all isotopes: */
  for(i=0; i<iso->n_i; i++)
      transitprint(1, verblevel, " %-8s: density %8g g/cm3\n",
                                   iso->isof[i].n,
                                   at->molec[iso->imol[i]].d[0]);
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
  transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
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
  transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
               "Line %i of file '%s': Field %i (%s) does not have a valid "
               "value:\n%s.\n", nmb, atmfilename, fld, fldn, line);
  exit(EXIT_FAILURE);
}


/* \fcnfh
   Check that val is positive. Throw error message if not.                  */
static inline void
checkposvalue(PREC_RES val, /* Value to check                               */
              int field,    /* Field where it was read                      */
              long line){   /* Line from which it was read                  */

  if(val<0)
    transiterror(TERR_SERIOUS,
                 "While reading the %i-th field in line %li of atmosphere "
                 "file %s, a negative value was found (%g).\n", field,
                 line-1, atmfilename, val);
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
        transitprint(15, verblevel, "The number of molecules is %d.\n", nmol);

        /* Allocate Molecules names:                                        */
        mol->name    = (char **)calloc(nmol,             sizeof(char *));
        mol->name[0] = (char  *)calloc(nmol*maxeisoname, sizeof(char));
        for(i=1; i<nmol; i++)
          mol->name[i] = mol->name[0] + i*maxeisoname;
        /* Read and store the molecule names:                               */
        transitprint(1, verblevel, "Molecules with abundance profile:\n  ");
        lp = line;
        /* Read and store names:                                            */
        for (i=0; i<nmol; i++){
          getname(lp, mol->name[i]);
          lp = nextfield(lp);
          if (i < nmol-1)
            transitprint(1, verblevel, "%s, ", mol->name[i]);
          else
            transitprint(1, verblevel, "%s.\n", mol->name[i]);
        }
      }
      continue;
    case 0:     /* Throw error if EOF:                                      */
      transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
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
        transiterror(TERR_SERIOUS,
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
        transiterror(TERR_SERIOUS, "Invalid unit factor indication in "
                                   "atmosphere file.\n");
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

  transitprint(3, verblevel, "Read all keywords in atmosphere file without "
                             "problems.\n");

  /* Set total number of molecules in atmosphere:                           */
  mol->nmol = at->n_aiso = nmol;

  /**/
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
  if(!nmol)
    transiterror(TERR_SERIOUS, "No species were found in the atmospheric file, "
                               "make sure to specify them with the comment/"
                               "header in the previous line '#SPECIES'.\n");

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

  transitprint(1, verblevel, "Start reading abundances.\n");
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
  int i, j,         /* Auxiliary for-loop indices                           */
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
        transitprint(30, verblevel, "density[%d, %li]: %.9f.\n",
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
      transiterror(TERR_WARNING, "In radius %i (%g km), abundances don't add "
                                 "up to 1.0: %.9g\n", r, at->rads.v[r], sumq);

    /* Calculate densities using ideal gas law:                             */
    if (r>=0){
      transitprint(30, verblevel, "Abund: %.9f, mmm: %.3f, mass: %.3f, "
                                "p: %.3f, T: %.3f.\n", molec[2].q[r], at->mm[r],
                                   mol->mass[2], at->atm.p[r]*at->atm.pfct,
                                   at->atm.t[r]*at->atm.tfct);
    }
    for(i=0; i<at->n_aiso; i++)
      molec[i].d[r] = stateeqnford(at->mass, molec[i].q[r], at->mm[r],
                                   mol->mass[i], at->atm.p[r]*at->atm.pfct,
                                   at->atm.t[r]*at->atm.tfct);
    transitprint(30, verblevel, "dens[%2li]: %.14f,   ", r, molec[2].d[r]);
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
  transitprint(1, verblevel, "Read %li layers between pressures %.3e and "
                             "%.3e barye.\n", r, at->atm.p[0  ]*at->atm.pfct,
                                                 at->atm.p[r-1]*at->atm.pfct);

  return nrad;
}


/* Read and store non-layer-dependent molecular data (mass, radius, ID)
   and store in mol struct.                                                 */
void
getmoldata(struct atm_data *at, struct molecules *mol, char *filename){
  int nmol = at->n_aiso;
  FILE *elist;

  double *mmass,       /* Molecular mass from list                          */
         *radius;      /* Molecular radii from list                         */
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

  /* Allocate molecule ID, mass, names, and radii:                          */
  molID    = (int    *)calloc(ndatamol,            sizeof(int));
  mmass    = (double *)calloc(ndatamol,            sizeof(double));
  radius   = (double *)calloc(ndatamol,            sizeof(double));
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
    lp = fgets(line, maxlinelen, elist);
    transitprint(30, verblevel, "Read species:'%6s',  ID:%3d,  mass:%7.4f,  "
      "radius:%5.2f.\n", rname[i], molID[i], mmass[i], radius[i]);
  }

  /* Assign info for each molecule in atosphere:                            */
  for (i=0; i<nmol; i++){
    /* Set the radius:                                                      */
    j = findstring(mol->name[i], rname, ndatamol);
    if (j == -1){
      transiterror(TERR_SERIOUS, "The atmospheric species '%s' is not "
                   "present in the list of known species:\n '%s'.\n",
                    mol->name[i], filename);
    }
    mol->radius[i] = radius[j] * ANGSTROM;
    /* Set the universal molecular ID:                                      */
    mol->ID[i]     = molID[j];
    /* Set the mass:                                                        */
    mol->mass[i]   = mmass[j];
    transitprint(1, verblevel, "Molecule %9s (ID: %3i) has radius %4.2f A "
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
    transitprint(20, verblevel, "%8.3f  ", at->atm.t[i]);
  }
  transitprint(20, verblevel, "\n");

  /* Update atmfile abundance arrays:                                       */
  for (j=0; j<nmol; j++)
    for (i=0; i<nlayers; i++)
      at->molec[j].q[i] = input[nlayers*(j+1) + i];

  /* Re-calculate the mean molecular mass and check whether abundances add up
     to one (within allowq threshold):                                    */
  for (i=0; i<nlayers; i++){
    sumq = checkaddmm(at->mm+i, i, at->molec, mol, nmol, at->mass);
    if(fabs(sumq-1.0) > allowq)
      transiterror(TERR_WARNING, "In radius %i (%g km), abundances don't add "
                                 "up to 1.0: %.9g\n", i, at->rads.v[i], sumq);

    /* Re-calculate densities:                                              */
    for(j=0; j<nmol; j++)
      at->molec[j].d[i] = stateeqnford(at->mass, at->molec[j].q[i], at->mm[i],
                                       mol->mass[j], at->atm.p[i]*at->atm.pfct,
                                                     at->atm.t[i]*at->atm.tfct);
  }

  /* Re-calculate radius:                                                   */
  transitprint(30, verblevel, "Old radius boundaries: [%.1f, %.1f]\n",
                   at->rads.v[0], at->rads.v[nlayers-1]);
  /* Check that r0, p0, and gsurf were defined:                             */
  if (tr->p0 == 0 || tr->r0 == 0 || tr->gsurf == 0)
    transiterror(TERR_SERIOUS, "Surface gravity (%.1f) or reference pressure "
      "(%.3e) or radius (%.1f) were not defined.\n", tr->gsurf, tr->p0, tr->r0);
  radpress(tr->gsurf, tr->p0, tr->r0, at->atm.t,
           at->mm, at->atm.p, at->rads.v, nlayers, at->rads.fct);
  /* Actualize other variables of at->rads:                                 */
  at->rads.i = at->rads.v[0];
  at->rads.f = at->rads.v[nlayers-1];
  transitprint(30, verblevel, "New radius boundaries: [%.1f, %.1f]\n",
                   tr->ds.at->rads.v[0], tr->ds.at->rads.v[nlayers-1]);

  /* Re-resample transit arrays:                                            */
  makeradsample(tr);
  return 0;
}


int radpress(double g,         /* Surface gravity (m/s^2)                   */
             double p0,        /* Reference pressure                        */
             double r0,        /* Reference height                          */
             double *temp,     /* Temperature (K)                           */
             double *mu,       /* Mean molecular mass (g/mol)               */
             double *pressure, /* Pressure array                            */
             double *radius,   /* Radius array                              */
             int nlayer,       /* Number of layers                          */
             double rfct){     /* Radius-array units                        */
  double rad0;   /* Radius at reference pressure                            */
  int i, i0=-1;  /* Indices                                                 */

  /* Start from zero altitude, adjust later:                                */
  radius[0] = 0.0;
  for (i=1; i<nlayer; i++){
    /* Cumulative trapezoidal integration to solve: dz = -H*dlog(p):        */
    radius[i] = radius[i-1] - 0.5*log(pressure[i]/pressure[i-1]) *
                              (KB/g/AMU) * (temp[i]/mu[i] + temp[i-1]/mu[i-1]);
    /* Find index of layers around p0 such:  press[i-1] <= p0 < press[i]:   */
    if (fabs(pressure[i-1]-p0) <  fabs(pressure[i]-pressure[i-1]))
      if (fabs(pressure[i]-p0) <= fabs(pressure[i]-pressure[i-1]))
        i0 = i-1;
  }

  /* Reference pressure is no in given range:                               */
  if (i0 == -1){
    transiterror(TERR_CRITICAL, "Reference pressure level (%.3e) not found "
               "in range [%.3e, %.3e].\n", p0, pressure[0], pressure[nlayer-1]);
    return 0;
  }

  /* Log-linearly interpolate to get radius at p0:                          */
  rad0 = radius[i0] + (radius[i0+1]-radius[i0]) *
         log(p0/pressure[i0]) / log(pressure[i0+1]/pressure[i0]);

  transitprint(20, verblevel, "PRESSURE:\n");
  for (i=0; i<nlayer; i++)
    transitprint(20,verblevel, "%.3e, ", pressure[i]);
  transitprint(20, verblevel, "\n");

  /* Adjust radius array to force radius(p0) = r0:                          */
  transitprint(20, verblevel, "RADIUS:\n");
  for (i=0; i<nlayer; i++){
    radius[i] = r0 + (radius[i] - rad0)/rfct;
    transitprint(20, verblevel, "%.1f, ", radius[i]);
  }
  transitprint(20, verblevel, "\n");
  return 1;
}
