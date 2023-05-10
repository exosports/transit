// Copyright (C) 2015-2016 University of Central Florida. All rights reserved.
// Transit is under an open-source, reproducible-research license (see LICENSE).

#include <transit.h>

/* \fcnfh
   Read Cross-section data from tabulated files.
   Return: 0 on success                                                     */
int
readcs(struct transit *tr){
  FILE *fp;       /* Pointer to cross-section file                          */
  char *file,     /* Cross-section file name                                */
       *colname;  /* Cross-section isotope names                            */
  PREC_CS **a,    /* Cross-section cross sections sample                    */
           *wn;   /* Cross-section sampled wavenumber array                 */

  static struct cross st_cross;  /* Cross-section structure                 */
  tr->ds.cross = &st_cross;
  /* Number of Cross-section files:                                         */
  int nfiles = tr->ds.cross->nfiles = tr->ds.th->ncross;
  long nt = 0, wa;        /* Number of temperature & wn samples in CSfile   */
  char rc;
  char *lp, *lpa;         /* Pointers in file                               */
  int maxline=3000,       /* Max length of line                             */
      j, k, n=0;          /* Counters                                       */
  long lines;             /* Lines read counter                             */
  long i,                 /* Auxiliary for indices                          */
       nspec;             /* Number of species in cross-section file        */
  char line[maxline+1];   /* Array to hold line being read                  */
  struct molecules *mol=tr->ds.mol;

  /* Make sure that radius and wavenumber samples exist:                    */
  transitcheckcalled(tr->pi, "interpcs", 2, "makewnsample", TRPI_MAKEWN,
                                            "makeradsample", TRPI_MAKERAD);

  /* Allocate Transit extinction array (in cm-1):                           */
  st_cross.e    = (PREC_CS **)calloc(tr->wns.n,            sizeof(PREC_CS *));
  st_cross.e[0] = (PREC_CS  *)calloc(tr->wns.n*tr->rads.n, sizeof(PREC_CS));
  for(j=1; j < tr->wns.n; j++)
    st_cross.e[j] = st_cross.e[0] + j*tr->rads.n;
  memset(st_cross.e[0], 0, tr->wns.n*tr->rads.n*sizeof(double));

  /* Min and max allowed temperatures in CS files:                          */
  st_cross.tmin =     0.0;
  st_cross.tmax = 70000.0;

  /* If there are no files, allocate tr.ds.cross.e (extinction) and return: */
  if(!nfiles){
    return 0;
  }
  tr_output(TOUT_RESULT, "Computing cross-section opacities for %i "
    "database%s:\n", nfiles, (nfiles > 1 ? "s" : ""));

  /* Allocate string for molecule names:                                    */
  colname = (char *)calloc(maxline, sizeof(char));

  /* Allocate species' ID array:                                            */
  st_cross.mol    = (int **)calloc(nfiles, sizeof(int *));
  st_cross.mol[0] = (int  *)calloc(nfiles*2, sizeof(int));
  for(i=1; i < nfiles; i++)
    st_cross.mol[i] = st_cross.mol[0] + i*2;

  /* Number of temperature and wavenumber samples per file:                 */
  st_cross.ntemp = (int  *)calloc(nfiles, sizeof(int));
  st_cross.nwave = (int  *)calloc(nfiles, sizeof(int));
  st_cross.nspec = (int  *)calloc(nfiles, sizeof(int));
  /* Temperature and wavenumber samples:                                    */
  st_cross.cs   = (PREC_CS ***)calloc(nfiles, sizeof(PREC_CS **));
  st_cross.temp = (PREC_CS  **)calloc(nfiles, sizeof(PREC_CS  *));
  st_cross.wn   = (PREC_CS  **)calloc(nfiles, sizeof(PREC_CS  *));

  for (j=0; j < nfiles; j++){
    /* Copy file names from hint:                                           */
    file = xstrdup(tr->ds.th->csfile[j]);

    /* Attempt to open the files:                                           */
    if((fp=fopen(file, "r")) == NULL) {
      tr_output(TOUT_ERROR, "Cannot read cross-section file '%s'.\n",file);
      exit(EXIT_FAILURE);
    }
    tr_output(TOUT_DEBUG,
      "  Cross-section file (%d/%d): '%s'\n", (j + 1), nfiles, file);

    lines = 0; /* lines read counter                                        */
    lpa   = 0;
    /* Read the file headers:                                               */
    while(1){
      /* Skip comments, blanks and read next line:                          */
      while((rc=fgetupto_err(lp=line, maxline, fp, &cserr, file, lines++))
             =='#' || rc=='\n');
      /* If it is end of file, stop loop:                                   */
      if(!rc) {
        tr_output(TOUT_ERROR,
          "File '%s' finished before opacity info.\n", file);
        exit(EXIT_FAILURE);
      }

      switch(rc){
      case 'i': /* Read the name of the isotopes:                           */
        while(isblank(*++lp));
        /* Count the number of species:                                     */
        nspec = st_cross.nspec[j] = countfields(lp, ' ');
        if (nspec != 1 && nspec != 2) {
          tr_output(TOUT_ERROR,
            "Wrong header in cross section file '%s', The 'i'-line "
            "should contain either one or two species, separated by "
            "blank spaces. The line reads:\n  '%s'\n", file, lp);
          exit(EXIT_FAILURE);
        }

        for (k=0; k<nspec; k++){
          /* Read the name of the species:                                  */
          getname(lp, colname);
          /* Find the ID of the species:                                    */
          for(i=0; i<mol->nmol; i++)
            if(strcmp(mol->name[i], colname)==0)
              st_cross.mol[j][k] = i;
          /* If the species is not in the atmosphere file:                  */
          if(st_cross.mol[j][k] == -1) {
            tr_output(TOUT_ERROR,
              "Cross-section species '%s' from file '%s' does not match "
              "any in the atmsopheric file.\n", colname, file);
            exit(EXIT_FAILURE);
          }
          lp = nextfield(lp);
        }

        tr_output(TOUT_DEBUG, "  Cross-section species: ");
        for (k=0; k<nspec; k++)
          tr_output(TOUT_DEBUG, "%s, ", mol->name[st_cross.mol[j][k]]);
        tr_output(TOUT_DEBUG, "\n");
        continue;

      case 't': /* Read the sampling temperatures array:                    */
        while(isblank(*++lp));
        nt = st_cross.ntemp[j] = countfields(lp, ' ');  /* Number of temps. */
        tr_output(TOUT_DEBUG, "  Number of temperature samples: %ld\n",
                                   nt);
        if(!nt) {
          tr_output(TOUT_ERROR,
            "Wrong line %i in cross-section file '%s', if it begins with "
            "a 't' then it should have the blank-separated fields with "
            "the temperatures. Rest of line: '%s'.\n", lines, file, lp);
          exit(EXIT_FAILURE);
        }

        /* Allocate and store the temperatures array:                       */
        st_cross.temp[j] = (PREC_CS *)calloc(nt, sizeof(PREC_CS));
        n = 0;    /* Count temperatures per line                            */
        lpa = lp; /* Pointer in line                                        */
        tr_output(TOUT_DEBUG, "  Temperatures (K) = [");
        while(n < nt){
          while(isblank(*lpa++));
          st_cross.temp[j][n] = strtod(--lpa, &lp);  /* Get value           */
          tr_output(TOUT_DEBUG, "%d, ", (int)st_cross.temp[j][n]);
          if(lp==lpa) {
            tr_output(TOUT_ERROR,
              "Less fields (%i) than expected (%i) were read for "
              "temperature in the cross-section file '%s'.\n", n, nt, file);
            exit(EXIT_FAILURE);
          }

          if((lp[0]|0x20) == 'k') lp++; /* Remove trailing K if exists      */
          lpa = lp;
          n++;
        }
        tr_output(TOUT_DEBUG, "\b\b]\n");
        continue;
      default:
        break;
      }
      break;
    }
    /* Set tmin and tmax:                                                   */
    st_cross.tmin = fmax(st_cross.tmin, st_cross.temp[j][  0]);
    st_cross.tmax = fmin(st_cross.tmax, st_cross.temp[j][n-1]);

    /* Set an initial value for allocated wavenumber fields:                */
    wa = 32;

    /* Allocate wavenumber array:                                           */
    wn   = (PREC_CS  *)calloc(wa,    sizeof(PREC_CS));
    /* Allocate input extinction array (in cm-1 amagat-2):                  */
    a    = (PREC_CS **)calloc(wa,    sizeof(PREC_CS *));
    a[0] = (PREC_CS  *)calloc(wa*nt, sizeof(PREC_CS));
    for(i=1; i<wa; i++)
      a[i] = a[0] + i*nt;

    n=0;
    /* Read information for each wavenumber sample:                         */
    while(1){
      /* Skip comments and blanks; read next line:                          */
      if (n)
        while((rc=fgetupto_err(lp=line, maxline, fp, &cserr, file, lines++))
              =='#'||rc=='\n');
      /* Stop, if it is end of file:                                        */
      if(!rc)
        break;

      /* Re-allocate (double the size) if necessary:                        */
      if(n==wa){
        wn   = (PREC_CS  *)realloc(wn,  (wa<<=1) * sizeof(PREC_CS));
        a    = (PREC_CS **)realloc(a,    wa *      sizeof(PREC_CS *));
        a[0] = (PREC_CS  *)realloc(a[0], wa * nt * sizeof(PREC_CS));
        for(i=1; i<wa; i++)
          a[i] = a[0] + i*nt;
      }

      /* Store new line: wavenumber first, then loop over cross sections:   */
      while(isblank(*lp++));
      wn[n] = strtod(lp-1, &lpa);  /* Store wavenumber                      */
      if(lp==lpa+1) {
        tr_output(TOUT_ERROR,
          "Invalid fields for the %ith wavenumber in the cross-section "
          "file '%s'.\n", n+1, file);
        exit(EXIT_FAILURE);
      }

      i = 0;
      while(i<nt){
        a[n][i] = strtod(lpa, &lp); /* Store cross-section extinction       */
        if(lp==lpa) {
          tr_output(TOUT_ERROR,
            "Less fields (%i) than expected (%i) were read for the %ith "
            "wavenumber in the cross-section file '%s'.\n", i, nt, n+1, file);
          exit(EXIT_FAILURE);
        }

        lpa = lp;
        i++;
      }
      n++;
    }

    /* Re-allocate arrays to their final sizes:                             */
    if(n<wa){
      st_cross.wn[j] = (PREC_CS  *)realloc(wn,   n*   sizeof(PREC_CS));
      a              = (PREC_CS **)realloc(a,    n*   sizeof(PREC_CS *));
      a[0]           = (PREC_CS  *)realloc(a[0], n*nt*sizeof(PREC_CS));
      for(i=1; i<n; i++)
        a[i] = a[0] + i*nt;
    }
    tr_output(TOUT_DEBUG, "  Number of wavenumber samples: %d\n", n);
    tr_output(TOUT_DEBUG, "  Wavenumber array (cm-1) = [%.1f, %.1f, "
      "%.1f, ..., %.1f, %.1f, %.1f]\n",
      st_cross.wn[j][  0], st_cross.wn[j][  1],
      st_cross.wn[j][  2], st_cross.wn[j][n-3],
      st_cross.wn[j][n-2], st_cross.wn[j][n-1]);

    /* Wavenumber boundaries check:                                         */
    if ((st_cross.wn[j][  0] > tr->wns.v[          0]) ||
        (st_cross.wn[j][n-1] < tr->wns.v[tr->wns.n-1]) ){
      tr_output(TOUT_ERROR,
        "The wavelength range [%.2f, %.2f] cm-1 of the cross-section "
        "file:\n  '%s',\ndoes not cover Transit's wavelength range "
        "[%.2f, %.2f] cm-1.\n", file, st_cross.wn[j][0], st_cross.wn[j][n-1],
        tr->wns.v[0], tr->wns.v[tr->wns.n-1]);
      exit(EXIT_FAILURE);
    }
    st_cross.cs[j] = a;
    st_cross.nwave[j] = n;
    fclose(fp);
  }
  free(colname);
  tr_output(TOUT_RESULT, "Done.\n");
  tr->pi |= TRPI_CS;
  return 0;
}


int
interpcs(struct transit *tr){
  struct molecules *mol=tr->ds.mol;
  struct cross     *cross=tr->ds.cross;
  prop_atm *atm = &tr->atm;
  PREC_CS **e;  /* Temporary interpolated cia                               */
  double *tmpw = malloc(tr->wns.n  * sizeof(double)); /* Temperatures array */
  double *tmpt = malloc(tr->rads.n * sizeof(double)); /* Wavenumber array   */
  double dens;  /* Density scaling factor                                   */
  int i, j, k, n,
      icsmol;  /* Cross-section species index                               */

  /* Reset the cross-section opacity to zero:                               */
  memset(cross->e[0], 0, tr->wns.n*tr->rads.n*sizeof(double));

  /* Allocate temporary array for opacity:                                  */
  e    = (PREC_CS **)calloc(tr->wns.n,            sizeof(PREC_CS *));
  e[0] = (PREC_CS  *)calloc(tr->wns.n*tr->rads.n, sizeof(PREC_CS));
  for(i=1; i < tr->wns.n; i++)
    e[i] = e[0] + i*tr->rads.n;

  /* Get transit temperatures and wavenumber arrays:                        */
  for(i=0; i<tr->rads.n; i++){
    tmpt[i] = atm->tfct * atm->t[i];
    /* Check for temperature boundaries:                                    */
    if (tmpt[i] < cross->tmin) {
      tr_output(TOUT_ERROR,
        "The layer %d in the atmospheric model has a lower temperature "
        "(%.1f K) than the lowest allowed cross-section temperature "
        "(%.1f K).\n", i, tmpt[i], cross->tmin);
      exit(EXIT_FAILURE);
    }
    if (tmpt[i] > cross->tmax) {
      tr_output(TOUT_ERROR,
        "The layer %d in the atmospheric model has a higher temperature "
        "(%.1f K) than the highest allowed  cross-section temperature "
        "(%.1f K).\n", i, tmpt[i], cross->tmax);
      exit(EXIT_FAILURE);
    }
  }
  for(i=0; i<tr->wns.n; i++)
    tmpw[i] = tr->wns.fct * tr->wns.v[i];

  for (n=0; n < cross->nfiles; n++){
    /* Interpolate data to transit's wavenumber and temperature sampling:   */
    bicubicinterpolate(e, cross->cs[n], cross->wn[n],   cross->nwave[n],
                                        cross->temp[n], cross->ntemp[n],
                                        tmpw, tr->wns.n, tmpt, tr->rads.n);

    /* Calculate absorption coefficients in cm-1 units:                     */
    for(i=0; i < tr->rads.n; i++){
      dens = 1.0;
      for(k=0; k < cross->nspec[n]; k++){
        icsmol = cross->mol[n][k];
        dens *= mol->molec[icsmol].d[i] / (AMU * mol->mass[icsmol] * AMAGAT);
      }
      for(j=0; j < tr->wns.n; j++){
        /* FINDME: Bicubic interplation on CIA files with zeros in the 
                   temperature/wavenumber grid can lead to bogus opacities.  
                   This bandage handles the case when the bogus opacities 
                   are negative.  bicubicinterpolate() should be changed to 
                   better handle this scenario.                             */
        if(e[j][i] > 0)
          cross->e[j][i] += e[j][i] * dens;
      }
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
  while(t1[fi++] < fx1);
  fi--;

  for(i=0; i<li; i++)
    if(t1[i] > lx1)
      li = i;

  while(t2[fj++] < fx2);
  fj--;

  for(j=0; j<lj; j++)
    if(t2[j] > lx2)
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
   Error printing function for lines longer than maxline in the CS file     */
void
cserr(int max,     /* Max line length                                       */
      char *name,  /* Cross-section file name                               */
      int line){   /* Line number                                           */
  tr_output(TOUT_ERROR,
    "Line %i of cross-section file '%s' is longer than %i characters ...\n "
    "hard coded values in file '%s' need to be changed.\n", line, name,
    max, __FILE__);
  exit(EXIT_FAILURE);
}


/* \fcnfh
   Free cross-section structure
   Return: 0 on success                                                     */
int
freemem_cs(struct cross *cross,
           long *pi){
  int i=0;
  free(cross->e[0]);
  free(cross->e);
  for (i=0; i<cross->nfiles; i++){
    free(cross->cs[i][0]);
    free(cross->cs[i]);
    free(cross->wn[i]);
    free(cross->temp[i]);
  }
  if(cross->nfiles){
    free(cross->cs);
    free(cross->wn);
    free(cross->temp);

    free(cross->mol[0]);
    free(cross->mol);
    free(cross->nwave);
    free(cross->ntemp);
    free(cross->nspec);
  }

  /* Unset appropiate flags:                                                */
  *pi &= ~(TRPI_CS);
  return 0;
}
