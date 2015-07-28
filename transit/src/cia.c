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

/* \fcnfh
   Read CIA info from tabulated files.
   Return: 0 on success                                                     */
int
readcia(struct transit *tr){
  FILE *fp;       /* Pointer to CIA file                                    */
  char *file,     /* CIA file name                                          */
       *colname;  /* CIA isotope names                                      */
  PREC_CIA **a,   /* CIA cross sections sample                              */
           *wn;   /* CIA sampled wavenumber array                           */

  static struct cia st_cia;  /* CIA structure                               */
  tr->ds.cia = &st_cia;
  int npairs = tr->ds.cia->nfiles = tr->ds.th->ncia; /* Number of CIA files */
  int p;                 /* Auxiliary wavenumber index                      */
  long nt = 0, wa;       /* Number of temperature, wn samples in CIA file   */
  char rc;
  char *lp, *lpa;        /* Pointers in file                                */
  int maxline=1000, n=0;  /* Max length of line. Counter                     */
  long lines;            /* Lines read counter                              */
  long i;                /* Auxiliary for indices                           */
  char line[maxline+1];  /* Array to hold line being read                   */
  struct molecules *mol=tr->ds.mol;

  /* Make sure that radius and wavenumber samples exist:                    */
  transitcheckcalled(tr->pi, "interpolatecia", 2, "makewnsample", TRPI_MAKEWN,
                                               "makeradsample", TRPI_MAKERAD);

  /* Allocate (output) transit extinction array (in cm-1):                  */
  st_cia.e    = (PREC_CIA **)calloc(tr->wns.n,            sizeof(PREC_CIA *));
  st_cia.e[0] = (PREC_CIA  *)calloc(tr->wns.n*tr->rads.n, sizeof(PREC_CIA));
  for(p=1; p < tr->wns.n; p++)
    st_cia.e[p] = st_cia.e[0] + p*tr->rads.n;
  memset(st_cia.e[0], 0, tr->wns.n*tr->rads.n*sizeof(double));

  /* If there are no files, allocate tr.ds.cia.e (extinction) and return:   */
  if(!npairs){
    return 0;
  }
  transitprint(1, verblevel, "Computing CIA opacities for %i database%s:\n",
               npairs, npairs>1 ? "s":"");

  /* Allocate string for molecule names:                                    */
  colname = (char *)calloc(maxline, sizeof(char));

  /* Allocate molecules' ID:                                                */
  st_cia.mol1 = (int   *)calloc(npairs, sizeof(int));
  st_cia.mol2 = (int   *)calloc(npairs, sizeof(int));
  /* Number of temperature and wavenumber samples per file:                 */
  st_cia.ntemp = (int  *)calloc(npairs, sizeof(int));
  st_cia.nwave = (int  *)calloc(npairs, sizeof(int));
  /* CIA, temperature and wavenumber samples:                               */
  st_cia.cia  = (PREC_CIA ***)calloc(npairs, sizeof(PREC_CIA **));
  st_cia.temp = (PREC_CIA  **)calloc(npairs, sizeof(PREC_CIA  *));
  st_cia.wn   = (PREC_CIA  **)calloc(npairs, sizeof(PREC_CIA  *));
  /* Min and max allowed temperatures in CIA files:                         */
  st_cia.tmin = -1.0;      /* (will update later as we read CIA files)      */
  st_cia.tmax = 100000.0;

  for (p=0; p < npairs; p++){
    /* Copy file names from hint:                                           */
    file = xstrdup(tr->ds.th->ciafile[p]);

    /* Attempt to open the files:                                           */
    if((fp=fopen(file, "r")) == NULL)
      transiterror(TERR_SERIOUS, "Cannot read CIA file '%s'.\n", file);
    transitprint(10, verblevel, "  CIA file (%d/%d): '%s'\n",
                                p+1, npairs, file);
    lines = 0; /* lines read counter                                        */
    lpa   = 0;
    /* Read the file headers:                                               */
    while(1){
      /* Skip comments, blanks and read next line:                          */
      while((rc=fgetupto_err(lp=line, maxline, fp, &ciaerr, file, lines++))
             =='#' || rc=='\n');
      /* If it is end of file, stop loop:                                   */
      if(!rc)
        transiterror(TERR_SERIOUS, "File '%s' finished before opacity info.\n",
                     file);

      switch(rc){
      case 'i': /* Read the name of the isotopes:                           */
        while(isblank(*++lp));
        /* Check that there are exactly two isotopes:                       */
        if(countfields(lp, ' ') != 2)
          transiterror(TERR_SERIOUS,
                       "Wrong line %i in CIA file '%s', if it begins with a "
                       "'i', it should have the species separated by blank "
                       "spaces.  Rest of line:\n'%s'\n", lines, file, lp);

        st_cia.mol1[p] = st_cia.mol2[p] = -1;
        /* Allocate and copy the name of the first moleculee:               */
        getname(lp, colname);
        /* Find the ID of the first molecule:                               */
        for(i=0; i<mol->nmol; i++)
          if(strcmp(mol->name[i], colname)==0)
            st_cia.mol1[p] = i;
        /* If the molecule is not in the atmosphere file:                   */
        if(st_cia.mol1[p] == -1)
          transiterror(TERR_SERIOUS, "CIA molecule '%s' from file '%s' does "
                    "not match any in the atmsopheric file.\n", colname, file);

        /* Allocate and store the name of the second isotope:               */
        lp = nextfield(lp);
        getname(lp, colname);
        for(i=0; i < mol->nmol; i++)
          if(strcmp(mol->name[i], colname)==0)
            st_cia.mol2[p] = i;
        if(st_cia.mol2[p] == -1)
          transiterror(TERR_SERIOUS, "CIA molecule '%s' from file '%s' does "
                    "not match any in the atmsopheric file.\n", colname, file);
        transitprint(10, verblevel, "  CIA molecules: [%s, %s]\n",
                         mol->name[st_cia.mol1[p]], mol->name[st_cia.mol2[p]]);
        continue;

      case 't': /* Read the sampling temperatures array:                    */
        while(isblank(*++lp));
        nt = st_cia.ntemp[p] = countfields(lp, ' '); /* Number of temps.    */
        transitprint(5, verblevel, "  Number of temperature samples: %ld\n",
                                   nt);
        if(!nt)
          transiterror(TERR_SERIOUS, "Wrong line %i in CIA file '%s', if it "
                       "begins with a 't' then it should have the "
                       "blank-separated fields with the temperatures. Rest "
                       "of line: %s.\n", lines, file, lp);
        /* Allocate and store the temperatures array:                       */
        st_cia.temp[p] = (PREC_CIA *)calloc(nt, sizeof(PREC_CIA));
        n = 0;    /* Count temperatures per line                            */
        lpa = lp; /* Pointer in line                                        */
        transitprint(20, verblevel, "  Temperatures (K) = [");
        while(n < nt){
          while(isblank(*lpa++));
          st_cia.temp[p][n] = strtod(--lpa, &lp); /* Get value */
          transitprint(20, verblevel, "%d, ", (int)st_cia.temp[p][n]);
          if(lp==lpa)
            transiterror(TERR_CRITICAL, "Less fields (%i) than expected (%i) "
                         "were read for temperature in the CIA file '%s'.\n",
                         n, nt, file);
          if((lp[0]|0x20) == 'k') lp++; /* Remove trailing K if exists      */
          lpa = lp;
          n++;
        }
        transitprint(20, verblevel, "\b\b]\n");
        continue;
      default:
        break;
      }
      break;
    }
    /* Set tmin and tmax:                                                   */
    st_cia.tmin = fmax(st_cia.tmin, st_cia.temp[p][  0]);
    st_cia.tmax = fmin(st_cia.tmax, st_cia.temp[p][n-1]);

    /* Set an initial value for allocated wavenumber fields:                */
    wa = 32;

    /* Allocate wavenumber array:                                           */
    wn   = (PREC_CIA  *)calloc(wa,    sizeof(PREC_CIA));
    /* Allocate input extinction array (in cm-1 amagat-2):                  */
    a    = (PREC_CIA **)calloc(wa,    sizeof(PREC_CIA *));
    a[0] = (PREC_CIA  *)calloc(wa*nt, sizeof(PREC_CIA));
    for(i=1; i<wa; i++)
      a[i] = a[0] + i*nt;

    n=0;
    /* Read information for each wavenumber sample:                         */
    while(1){
      /* Skip comments and blanks; read next line:                          */
      if (n)
        while((rc=fgetupto_err(lp=line, maxline, fp, &ciaerr, file, lines++))
              =='#'||rc=='\n');
      /* Stop, if it is end of file:                                        */
      if(!rc)
        break;

      /* Re-allocate (double the size) if necessary:                        */
      if(n==wa){
        wn   = (PREC_CIA  *)realloc(wn,  (wa<<=1) * sizeof(PREC_CIA));
        a    = (PREC_CIA **)realloc(a,    wa *      sizeof(PREC_CIA *));
        a[0] = (PREC_CIA  *)realloc(a[0], wa * nt * sizeof(PREC_CIA));
        for(i=1; i<wa; i++)
          a[i] = a[0] + i*nt;
      }

      /* Store new line: wavenumber first, then loop over cross sections:   */
      while(isblank(*lp++));
      wn[n] = strtod(lp-1, &lpa);  /* Store wavenumber                      */
      if(lp==lpa+1)
        transiterror(TERR_CRITICAL, "Invalid fields for the %ith wavenumber "
                                    "in the CIA file '%s'.\n", n+1, file);
      i = 0;
      while(i<nt){
        a[n][i] = strtod(lpa, &lp); /* Store cross section                  */
        if(lp==lpa)
          transiterror(TERR_CRITICAL, "Less fields (%i) than expected (%i) "
                       "were read for the %ith wavenumber in the CIA "
                       "file '%s'.\n", i, nt, n+1, file);
        lpa = lp;
        i++;
      }
      n++;
    }

    /* Re-allocate arrays to their final sizes:                             */
    if(n<wa){
      st_cia.wn[p] = (PREC_CIA  *)realloc(wn,   n*   sizeof(PREC_CIA));
      a            = (PREC_CIA **)realloc(a,    n*   sizeof(PREC_CIA *));
      a[0]         = (PREC_CIA  *)realloc(a[0], n*nt*sizeof(PREC_CIA));
      for(i=1; i<n; i++)
        a[i] = a[0] + i*nt;
    }
    transitprint(5, verblevel, "  Number of wavenumber samples: %d\n", n);
    transitprint(5, verblevel, "  Wavenumber array (cm-1) = [%.1f, %.1f, "
         "%.1f, ..., %.1f, %.1f, %.1f]\n", st_cia.wn[p][0],   st_cia.wn[p][1],
      st_cia.wn[p][2], st_cia.wn[p][n-3], st_cia.wn[p][n-2], st_cia.wn[p][n-1]);
    /* Wavenumber boundaries check:                                         */
    if ((st_cia.wn[p][  0] > tr->wns.v[          0]) ||
        (st_cia.wn[p][n-1] < tr->wns.v[tr->wns.n-1]) ){
      transiterror(TERR_SERIOUS, "The wavelength range [%.2f, %.2f] cm-1 of "
                   "the CIA file:\n  '%s',\ndoes not cover Transit's "
                   "wavelength range [%.2f, %.2f] cm-1.\n", file,
                   st_cia.wn[p][0], st_cia.wn[p][n-1],
                   tr->wns.v[0],    tr->wns.v[tr->wns.n-1]);
    }
    st_cia.cia[p] = a;
    st_cia.nwave[p] = n;
    fclose(fp);
  }
  /* FINDME: The program breaks when I free colname, it makes no sense      */
  free(colname);
  transitprint(1, verblevel, "Done.\n");
  tr->pi |= TRPI_CIA;
  return 0;
}


int
interpolatecia(struct transit *tr){
  struct molecules *mol=tr->ds.mol;
  struct cia       *cia=tr->ds.cia;
  prop_atm *atm = &tr->atm;
  PREC_CIA **e; /* Temporary interpolated cia                               */
  double *tmpw = malloc(tr->wns.n  * sizeof(double)); /* temperatures array */
  double *tmpt = malloc(tr->rads.n * sizeof(double)); /* wavenumber array   */
  double *densiso1, *densiso2, amagat2; /* Density arrays and square amagat */
  int i, j, n;

  /* Reset CIA opacity to zero:                                             */
  memset(cia->e[0], 0, tr->wns.n*tr->rads.n*sizeof(double));

  /* Allocate temporary array for opacity:                                  */
  e    = (PREC_CIA **)calloc(tr->wns.n,            sizeof(PREC_CIA *));
  e[0] = (PREC_CIA  *)calloc(tr->wns.n*tr->rads.n, sizeof(PREC_CIA));
  for(i=1; i < tr->wns.n; i++)
    e[i] = e[0] + i*tr->rads.n;

  /* Get transit temperatures and wavenumber arrays:                        */
  for(i=0; i<tr->rads.n; i++){
    tmpt[i] = atm->tfct * atm->t[i];
    /* Check for temperature boundaries:                                    */
    if (tmpt[i] < cia->tmin)
      transiterror(TERR_SERIOUS, "The layer %d in the atmospheric model has "
                   "a lower temperature (%.1f K) than the lowest allowed "
                   "CIA temperature (%.1f K).\n", i, tmpt[i], cia->tmin);
    if (tmpt[i] > cia->tmax)
      transiterror(TERR_SERIOUS, "The layer %d in the atmospheric model has "
                   "a higher temperature (%.1f K) than the highest allowed "
                   "CIA temperature (%.1f K).\n", i, tmpt[i], cia->tmax);
  }
  for(i=0; i<tr->wns.n; i++)
    tmpw[i] = tr->wns.fct * tr->wns.v[i];

  for (n=0; n < cia->nfiles; n++){
    /* Interpolate data to transit's wavenumber and temperature sampling:   */
    bicubicinterpolate(e, cia->cia[n], cia->wn[n],   cia->nwave[n],
                                       cia->temp[n], cia->ntemp[n],
                                       tmpw, tr->wns.n, tmpt, tr->rads.n);

    /* Get density profile of isotopes from atmosphere-file isotopes:       */
    densiso1 = mol->molec[cia->mol1[n]].d;
    densiso2 = mol->molec[cia->mol2[n]].d;

    /* Calculate absorption coefficients in cm-1 units:                     */
    for(i=0; i < tr->rads.n; i++){
      amagat2 = densiso1[i]*densiso2[i]/(AMU*mol->mass[cia->mol1[n]] * AMAGAT *
                                         AMU*mol->mass[cia->mol2[n]] * AMAGAT);
      for(j=0; j < tr->wns.n; j++)
        cia->e[j][i] += e[j][i]*amagat2;
    }
  }
  free(e[0]);
  free(e);
  free(tmpw);
  free(tmpt);

  return 0;
}


/* \fcnfh
   Interpolates 'src' into 'res' according to the new dimensions, first
   interpolates the second dimension and then the first. The result is
   added to whatever it is already existent in 'res'

   Return: 0 on success                                                     */
int
bicubicinterpolate(double **res,  /* target array [t1][t2]                  */
                   double **src,  /* Source array [x1][x2]                  */
                   double *x1,    /* Source first array                     */
                   long nx1,      /* Size of x1                             */
                   double *x2,    /* Source second array                    */
                   long nx2,      /* Size of x2                             */
                   double *t1,    /* Requested first array                  */
                   long nt1,      /* Size of t1                             */
                   double *t2,    /* Requested second array                 */
                   long nt2){     /* Size of t2                             */

  long i, j;           /* Auxiliary for-loop indices                        */
  double fx1 = x1[0],  /* First and last values of source arrays:           */
         fx2 = x2[0],
         lx1 = x1[nx1-1],
         lx2 = x2[nx2-1];
  long lj=nt2, fj=0;  /* Indices of edges of 2nd dimension                  */
  long li=nt1, fi=0;  /* Indices of edges of 1st dimension                  */
  double *z1, *z2;

  memset(res[0], 0, nt1*nt2*sizeof(double));
  /* Return if sampling regions don't match:                                */
  if(t1[0]>lx1 || t1[nt1-1]<fx1 || t2[0]>lx2 || t2[nt2-1]<fx2)
    return 0;

  /* Find indices where requested array values are within source boundaries.
     (i.e., do not extrapolate):                                            */
  while(t1[fi++]<fx1);
  fi--;

  for(i=0; i<li; i++)
    if(t1[i]>lx1)
      li = i;

  while(t2[fj++]<fx2);
  fj--;

  for(j=0; j<lj; j++)
    if(t2[j]>lx2)
      lj = j;

  /* Arrays created by spline_init to be used in interpolation calculation: */
  z1 = calloc(nx2, sizeof(double));
  z2 = calloc(nx1, sizeof(double));

  /* Temporary middle array to hold data that has been interpolated in one
     direction:                                                             */
  double **f2 = (double **)malloc(nt2    *sizeof(double *));
  f2[0]       = (double  *)malloc(nt2*nx1*sizeof(double  ));
  for(i=1; i<nt2; i++)
    f2[i] = f2[0] + i*nx1;

  /* Interpolate the 2nd dimension:                                         */
  for(i=0; i<nx1; i++){
    spline_init(z1, x2, src[i], nx2);
    for(j=fj; j<lj; j++){
      f2[j][i] = splinterp_pt(z1, nx2, x2, src[i], t2[j]);
    }
  }

  /* Interpolate the 1st dimension:                                         */
  for(j=fj; j<lj; j++){
    spline_init(z2, x1, f2[j], nx1);
    for(i=fi; i<li; i++){
      res[i][j] += splinterp_pt(z2, nx1, x1, f2[j], t1[i]);
    }
  }

  free(f2[0]);
  free(f2);
  free(z1);
  free(z2);

  return 0;
}


/* \fcnfh
   Error printing function for lines longer than maxline in the CIA file    */
void
ciaerr(int max,    /* Max line length                                       */
       char *name, /* CIA file name                                         */
       int line){  /* Line number                                           */
  transiterror(TERR_SERIOUS,
               "Line %i of CIA file '%s' is longer than %i characters ...\n"
               " hard coded values in file '%s' need to be changed.\n",
               line, name, max, __FILE__);
}


/* \fcnfh
   Free cia structure
   Return: 0 on success                                                     */
int
freemem_cia(struct cia *cia,
            long *pi){
  int i=0;
  /* Free arrays                                                            */
  free(cia->e[0]);
  free(cia->e);
  for (i=0; i<cia->nfiles; i++){
    free(cia->cia[i][0]);
    free(cia->cia[i]);
    free(cia->wn[i]);
    free(cia->temp[i]);
  }
  if(cia->nfiles){
    free(cia->cia);
    free(cia->wn);
    free(cia->temp);

    free(cia->mol1);
    free(cia->mol2);
    free(cia->nwave);
    free(cia->ntemp);
  }

  /* Unset appropiate flags:                                                */
  *pi &= ~(TRPI_CIA);
  return 0;
}
