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

#define CIA_DOFLOAT  2
#define CIA_RADFIRST 1

/* FUNCTION
   Initialize the optical depth structure for transit or eclipse geometry   */
int
init_optdepth(struct transit *tr){
  struct transithint *th=tr->ds.th;  /* Transit hint structure              */
  long int nwn  = tr->wns.n;         /* Number of wavenumbers               */
  long int i;                        /* For counting angles                 */
  int an = tr->ann;                  /* Number of angles                    */
  static struct optdepth tau;        /* Optical depth                       */
  static struct grid intens;         /* Intensity grid                      */

  long int nrad = tr->rads.n;  /* Number of layers                          */
  if (strcmp(tr->sol->name, "transit") == 0)
    nrad = tr->ips.n;          /* Number of impact parameter samples        */

  /* Initialize tau (optical depth) structure:                              */
  tr->ds.tau    = &tau;
  /* Set maximum optical depth:                                             */
  if(th->toomuch > 0)
    tau.toomuch = th->toomuch;

  /* Allocate array with layer index where tau reaches toomuch:             */
  tau.last = (long      *)calloc(nwn,      sizeof(long));
  /* Allocate optical-depth array [rad][wn]:                                */
  tau.t    = (PREC_RES **)calloc(nwn,      sizeof(PREC_RES *));
  tau.t[0] = (PREC_RES  *)calloc(nwn*nrad, sizeof(PREC_RES  ));
  for(i=1; i<nwn; i++)
    tau.t[i] = tau.t[0] + i*nrad;

  /* Eclipse-only structures:                                               */
  if (strcmp(tr->sol->name, "eclipse") == 0){
    /* Initialize intensity grid structure:                                 */
    tr->ds.intens = &intens;
    memset(&intens, 0, sizeof(struct grid));

    /* Allocate 2D array of intensities [angle][wn]:                        */
    intens.a    = (PREC_RES **)calloc(an,     sizeof(PREC_RES *));
    intens.a[0] = (PREC_RES  *)calloc(an*nwn, sizeof(PREC_RES  ));
    for(i=1; i<an; i++)
      intens.a[i] = intens.a[0] + i*nwn;
  }

  return 0;
}


/* FUNCTION
   Calculate the optical depth as a function of radii for a spherically
   symmetric planet.
   Return: 0 on success                                                     */
int
tau(struct transit *tr){
  struct transithint *th = tr->ds.th;    /* transithint struct              */
  struct optdepth *tau=tr->ds.tau;       /* Def optical depth structure     */

  struct extinction *ex = tr->ds.ex;     /* Extinction struct               */
  PREC_RES **e = ex->e;                  /* Extinction coefficient          */
  PREC_RES (*fcn)() = tr->sol->optdepth; /* eclipsetau or transittau func.  */

  long wi, ri = 0; /* Indices for wavenumber, and radius                    */
  int rn;          /* Functions output code                                 */
  int i;

  FILE *totEx   = NULL,
       *cloudEx = NULL,
       *scattEx = NULL;

  PREC_ATM *density = (PREC_ATM *)calloc(tr->ds.mol->nmol, sizeof(PREC_ATM));
  double   *Z       = (double   *)calloc(tr->ds.iso->n_i,  sizeof(double));

  prop_samp *rad = &tr->rads;  /* Radius sampling                           */
  PREC_RES *r  = rad->v;       /* Radius array                              */
  long int rnn = rad->n;       /* Number of layers                          */
  double  rfct = rad->fct;     /* Radius array units factor                 */

  /* Store the height of each layer (eclipse) or impact parameter (transit)
     starting from the outermost layer:                                     */
  PREC_RES *h;
  long int nh; /* Number of layers / impact-parameter samples               */
  double hfct;
  if (strcmp(tr->sol->name, "eclipse") == 0){
    h = (PREC_RES *)calloc(rnn, sizeof(PREC_RES)); /* Reversed radius array */
    for (ri=0; ri <rnn; ri++)
      h[ri] = r[rnn-ri-1];
    nh = rnn;     /* Number of layers                        */
    hfct = rfct;
  }
  else{
    prop_samp *ip = &tr->ips;
    h  = ip->v;  /* Impact parameter array                   */
    nh = ip->n;  /* Number of impact parameter samples       */
    hfct = ip->fct;
  }
  /* Request at least four layers to calculate a spline interpolation:     */
  if(nh < 4)
    transiterror(TERR_SERIOUS, "At least four layers (%d given) are required "
                 "(three for spline, one for the analitical part).\n", nh);

  PREC_RES *tau_wn;              /* Optical depth array                     */
  PREC_ATM *temp = tr->atm.t,    /* Temperature array                       */
           tfct  = tr->atm.tfct; /* Temperature units                       */

  prop_samp *wn = &tr->wns; /* Wavenumber sampling                          */
  long int wnn  = wn->n;    /* Number of wavenumber samples                 */
  double wfct   = wn->fct;  /* Wavenumber units factor                      */

  PREC_RES er[rnn];   /* Array of extinction per radius                     */
  int lastr = rnn-1;  /* Radius index of last computed extinction           */

  int wnextout = (long)(wnn/10.0); /* (Wavenumber sample size)/10,
                                      used for progress printing            */

  double e_s[rnn],                  /* Extinction from scattering           */
         e_c[rnn];                  /* Extinction from clouds               */
  PREC_CIA **e_cia = tr->ds.cia->e; /* Extinction from CIA                  */
  struct extscat *sc = tr->ds.sc;   /* Scattering extinction struct         */

  /* FINDME: TRANSIT ONLY */
  tr->save.ext = th->save.ext;  /* Ext. coefficient save/restore filename   */

  /* Check idxrefrac and extwn have been called:                            */
  transitcheckcalled(tr->pi, "tau", 2, "idxrefrac", TRPI_IDXREFRAC,
                                       "extwn",     TRPI_EXTWN);

  /* TRU_*BITS mark the bits that are carrying flags for each aspect.
  i.e. TRU_EXTBITS is 1 for all the flag positions that are relevant
  for the extinction computation (only TRU_OUTTAU in that case);
  TRU_ATMBITS is 1 for all the flag positions that are relevant for
  the atmospheric computation; and so  on...therefore the following line
  passes all the TAU relevant flags from the user hint to
  the main structure.  In particular, it only passes whether TRU_OUTTAU
  was 1 or 0 as specified by the user. */
  transitacceptflag(tr->fl, th->fl, TRU_TAUBITS);

  /* Set cloud structure:                                                   */
  static struct extcloud cl;
  cl.cloudext = th->cl.cloudext;  /* Maximum cloud extinction               */
  cl.cloudtop = th->cl.cloudtop;  /* Top layer radius                       */
  cl.cloudbot = th->cl.cloudbot;  /* Radius of maxe                         */
  tr->ds.cl = &cl;

  /* Has the extinction coefficient been calculated boolean:                */
  _Bool *comp = ex->computed;

  /* Restore extinction savefile if exists:                                 */
  if(tr->save.ext)
    restfile_extinct(tr->save.ext, e, comp, rnn, wnn);

  /* Compute extinction at the outermost layer:                             */
  if(!comp[rnn-1]){
    transitprint(1, verblevel, "Computing extinction at outermost layer.\n");
    if (tr->fp_opa != NULL)
      rn = interpolmolext(tr, rnn-1, ex->e);
    else{
      for (i=0; i < tr->ds.mol->nmol; i++)
        density[i] = tr->ds.mol->molec[i].d[rnn-1];
      for (i=0; i < tr->ds.iso->n_i; i++)
        Z[i]       = tr->ds.iso->isov[i].z [rnn-1];

      if((rn=computemolext(tr, ex->e+(rnn-1), tr->atm.t[rnn-1]*tr->atm.tfct,
                           density, Z, 0)) != 0)
        transiterror(TERR_CRITICAL,  "computemolext() returned error "
                                     "code %i.\n", rn);
    }
    ex->computed[rnn-1] = 1;
  }

  /* Save total, cloud, and scattering extinction to file if requested:     */
  if (th->savefiles){
    totEx = openFile("total_extion.dat",
              "# 2D total extinction\n"
              "# er [wn][rad]; wn[0]=min(wn), row[0]=bottom (max(p))\n");
    cloudEx = openFile("cloud_extion.dat",
              "# 2D cloud extinction\n"
              "# e_c [wn][rad]; wn[0]=min(wn), row[0]=bottom (max(p))\n");
    scattEx = openFile("scatt_extion.dat",
              "# 2D scatt extinction\n"
              "# e_s [wn][rad]; wn[0]=min(wn), row[0]=bottom (max(p))\n");
  }

  transitprint(1, verblevel, "Calculating optical depth at various radii:\n");
  /* For each wavenumber:                                                   */
  for(wi=0; wi<wnn; wi++){
    tau_wn = tau->t[wi];

    /* Print output every 10% progress:                                     */
    if(wi > wnextout){
      transitprint(2, verblevel, "%i%%\n", (int)(100*(float)wi/wnn+0.5));
      wnextout += (long)(wnn/10.0);
    }

    /* Calculate extinction from scattering, clouds, and CIA at each level: */
    computeextscat(e_s,  rnn, sc, rad->v, rad->fct, temp, tfct, wn->v[wi]*wfct);
    computeextcloud(e_c, rnn, &cl, rad, temp, tfct, wn->v[wi]*wfct);

    /* Put the extinction values in a new array, the values may be
       temporarily overwritten by (fcn)(), but they should be restored:     */
    for(ri=0; ri < rnn; ri++)
      er[ri] = e[ri][wi] + e_s[ri] + e_c[ri] + e_cia[wi][ri];

    /* For each height:                                                     */
    for(ri=0; ri < nh; ri++){
      /* Compute extinction at new radius if the impact parameter is smaller
         than the radius of last calculated extinction:                     */
      if(h[ri]*hfct < r[lastr]*rfct){
        /* FINDME: What if the ray ends up going through a lower layer because
           of the refraction? */
        if(ri)
          transitprint(30, verblevel, "Last Tau (height=%9.4g, wn=%9.4g): "
                               "%10.4g.\n", h[ri-1], wn->v[wi], tau_wn[ri-1]);
        /* While the extinction at a radius bigger than the impact
           parameter is not computed, go for it: */
        do{
          if(!comp[--lastr]){
            /* Compute extinction at given radius:                          */
            transitprint(5, verblevel, "Radius %i: %.9g cm ... \n",
                                       lastr+1, r[lastr]*rfct);
            if (tr->fp_opa != NULL)
              rn = interpolmolext(tr, lastr, ex->e);
            else{
              for (i=0; i < tr->ds.mol->nmol; i++)
                density[i] = tr->ds.mol->molec[i].d[lastr];
              for (i=0; i < tr->ds.iso->n_i; i++)
                Z[i]       = tr->ds.iso->isov[i].z [lastr];
              if((rn=computemolext(tr, ex->e+lastr,
                         tr->atm.t[lastr]*tr->atm.tfct, density, Z, 0)) != 0)
                transiterror(TERR_CRITICAL,  "computemolext() returned error "
                                             "code %i.\n", rn);
            }
            ex->computed[lastr] = 1;
            /* Update the value of the extinction at the right place:       */
            er[lastr] = e[lastr][wi] + e_s[lastr] + e_c[lastr] +
                        e_cia[wi][lastr];
          }
        }while(h[ri]*hfct < r[lastr]*rfct);
      }
      /* :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: */
      /* Calculate the optical depth (call to transittau or eclipsetau):    */
      tau_wn[ri] = rfct * fcn(tr, h[ri]*hfct/rfct, er);

      /* Check if the optical depth reached toomuch:                        */
      if (tau_wn[ri] > tau->toomuch){
        tau->last[wi] = ri;   /* Set tau.last                               */
        if (ri < 3){
          transiterror(TERR_WARNING, "At wavenumber %g (cm-1), the optical "
                       "depth (%g) exceeded toomuch (%g) at the height "
                       "level %li (%g km), this should have happened in a "
                       "deeper layer.\n", wn->v[wi],
                       tau_wn[ri], tau->toomuch, ri, h[ri]*hfct/1e5);
        }
        break;  /* Exit height loop when the optical depth reached toomuch  */
      }
      transitDEBUG(22, verblevel, "Tau(lambda %li=%9.07g, r=%9.4g) : %g "
          "(toomuch: %g)\n", wi, wn->v[wi], r[ri], tau_wn[ri], tau->toomuch);
    }

    /* Write total, cloud, and scattering extinction to file if requested:  */
    if (th->savefiles){
      save1Darray(tr, totEx,    er, rnn, wi);
      save1Darray(tr, cloudEx, e_c, rnn, wi);
      save1Darray(tr, scattEx, e_s, rnn, wi);
    }

    if(ri==nh){
      transitprint(1, verblevel, "WARNING: At wavenumber %g cm-1, tau reached "
                   "the bottom of the atmosphere with tau: %g (tau max: %g).\n",
                   wn->v[wi], tau_wn[ri-1], tau->toomuch);
      tau->last[wi] = ri-1;
    }
  }
  transitprint(1, verblevel, "Done.\n");

  /* Save various files if requested in the config file:                 */

  /* Save files requested                                                */
  if (th->savefiles){
    /* 2D tau [wn][rad], wn[0] = min(wn), rad[0]=top of the atm (min(p)) */
    savetau(tr);
    /* 2D CIA  extinction [wn][rad],
                              wn[0] = min(wn), rad[0]=bottom (max(p))    */
    saveCIA(tr);
    /* 2D mol-line extinction [wn][rad],
                              wn[0] = min(wn), rad[0]=bottom (max(p))    */
    savemolExtion(tr, ri);
    /* 2D total extinction [wn][rad],
                               wn[0] = min(wn), rad[0]=bottom (max(p))   */
    closeFile(totEx);
    /* 2D cloud extinction [wn][rad],
                               wn[0] = min(wn), rad[0]=bottom (max(p))   */
    closeFile(cloudEx);
    /* 2D scatt extinction [wn][rad],
                              wn[0] = min(wn), rad[0]=bottom (max(p))   */
    closeFile(scattEx);
  }

  /* Print detailed output if requested:                                 */
  if(tr->ds.det->tau.n)
    detailout(&tr->wns, &tr->ips,  &tr->ds.det->tau, tau->t, 0);
  if(tr->ds.det->ext.n)
    detailout(&tr->wns, &tr->rads, &tr->ds.det->ext, e, CIA_RADFIRST);
  if(tr->ds.det->cia.n)
    detailout(&tr->wns, &tr->rads, &tr->ds.det->cia, (double **)e_cia,
              CIA_DOFLOAT);

  if(tr->save.ext)
    savefile_extinct(tr->save.ext, e, comp, rnn, wnn);

  /* Print lowest impact parameter before optical depth gets too big:       */
  if(tr->f_toomuch)
    printtoomuch(tr->f_toomuch, tr->ds.tau, &tr->wns, &tr->ips);

  /* Set progress indicator and output tau if requested:                    */
  tr->pi |= TRPI_TAU;

  /* Free allocated memory:                                                 */
  free(density);
  free(Z);
  if (strcmp(tr->sol->name, "eclipse") == 0)
    free(h);
  return 0;
}


/* \fcnfh
   Generalized function to print 1D array to a file                         */
void print1dArrayDouble(FILE *outf, double *array, int noColumns, char *format){
	if(outf == NULL)
		outf= stdout;
	for(int col=0; col<noColumns; col++)
		fprintf(outf, format, array[col]);
	fprintf(outf, "\n");
}

/* \fcnfh
   Generalized function to print 2D array to a file                         */
void print2dArrayDouble(FILE *outf, double **array,
     int noRows, int noColumns, char *format, prop_samp *wn){
	if(outf == NULL)
		outf= stdout;
	for(int row=0; row < noRows; row++){
		fprintf(outf, "wavenumber: ");
		fprintf(outf, format, wn->v[row]);
		fprintf(outf, "\n");
		print1dArrayDouble(outf, array[row], noColumns, format);
		fprintf(outf, "\n");
	}
}


/* \fcnfh
  Print to a file molecular line extinction                                */
void
savemolExtion(struct transit *tr, long ri){
  struct extinction *ex = tr->ds.ex;     /* Extinction struct              */
  PREC_RES **e = ex->e;                   /* Extinction coefficient         */

  prop_samp *wn = &tr->wns;   /* Wavenumber sampling                        */
  long int wnn = wn->n;      /* Number of wavenumber samples               */
  prop_samp *rad = &tr->rads; /* Radius sampling                            */
  long int rnn = rad->n;     /* Number of layers                           */

  /* open file to write                                                     */
  FILE *myFile = fopen("mol_extion.dat", "w");

  /* format of the characters written                                       */
  char *format= "%-20.10g";

  /* write header                                                           */
  fprintf(myFile, "\n");
  fprintf(myFile, "# mol-line extinction\n");
  fprintf(myFile, "# e [rad][wn]; rad[0]=bottom (max(p)); wn[0]=min(wn)\n");
  fprintf(myFile, "\n");

  /* write file, row --> [rnn], column --> [wnn]                            */
  for(int ri=0; ri < rnn; ri++){
	fprintf(myFile, "radius: %-20.10g\n", rad->v[ri]);
	print1dArrayDouble(myFile, e[ri], wnn, format);
	fprintf(myFile, "\n");
  }
  /* close the file                                                         */
  fflush(myFile);
  fclose(myFile);
}


/* \fcnfh
  Print to a file CIA                                                      */
void
saveCIA(struct transit *tr){
  PREC_CIA **e_cia = tr->ds.cia->e; /* Extinction from CIA                  */

  prop_samp *wn = &tr->wns;   /* Wavenumber sampling                        */
  long int wnn = wn->n;      /* Number of wavenumber samples               */
  prop_samp *rad = &tr->rads; /* Radius sampling                            */
  long int rnn = rad->n;     /* Number of layers                           */

  /* open file to write                                                     */
  FILE *myFile = fopen("CIA.dat", "w");

  /* format of the characters written                                       */
  char *format= "%-20.10g";

  /* write header                                                           */
  fprintf(myFile, "\n");
  fprintf(myFile, "# 2D CIA extinction\n");
  fprintf(myFile, "# e_cia [wn][rad]; wn[0]=min(wn); row[0]=bottom (max(p))\n");
  fprintf(myFile, "\n");

  /* call 2D array function, row --> [wnn], column --> [rnn]                */
  print2dArrayDouble(myFile, e_cia, wnn, rnn, format, wn);

  /* close the file                                                         */
  fflush(myFile);
  fclose(myFile);
}


/* \fcnfh
  Print to a file tau                                                      */
void
save1Darray(struct transit *tr, FILE *myFile, PREC_RES *array1d,
                                              int nrad, long wi){
  prop_samp *wn = &tr->wns;   /* Wavenumber sampling                        */

  /* format of the characters written                                       */
  char *format= "%-20.10g";

  /*                                                                        */
  fprintf(myFile, "\n");
  fprintf(myFile, "wavenumber: %-20.10g\n", wn->v[wi]);
  print1dArrayDouble(myFile, array1d, nrad, format);
}

/* \fcnfh
  Print to a file tau                                                      */
FILE *
openFile(char *filename, char *header){
  /* open file to write                                                     */
  FILE *myFile = fopen(filename, "w");
  fprintf(myFile, "\n");
  fputs(header, myFile);
  return myFile;
}

/* \fcnfh
  Print to a file tau                                                      */
void
closeFile(FILE *myFile){
  /* close the file                                                         */
  fflush(myFile);
  fclose(myFile);
}


/* \fcnfh
  Print 2D tau into a file                                                  */
void
savetau(struct transit *tr){
  struct optdepth *tau=tr->ds.tau;       /* Def optical depth structure    */
  prop_samp *wn = &tr->wns;   /* Wavenumber sampling                        */
  long int wnn = wn->n;      /* Number of wavenumber samples               */
  prop_samp *rad = &tr->rads; /* Radius sampling                            */
  long int rnn = rad->n;     /* Number of layers                           */

  /* open file to write                                                     */
  FILE *myFile = fopen("tau.dat", "w");

  /* format of the characters written                                       */
  char *format= "%-20.10g";

  /* write header                                                           */
  fprintf(myFile, "\n");
  fprintf(myFile, "# 2D optical depth\n");
  fprintf(myFile, "# tau [wn][rad]; wn[0]=min(wn); rad[0]=top (min(p))\n");
  fprintf(myFile, "\n");

  /* call 2D array function, row --> [wnn], column --> [rnn]                */
  print2dArrayDouble(myFile, tau->t, wnn, rnn, format, wn);

  /* close the file                                                         */
  fflush(myFile);
  fclose(myFile);
}


/* FUNCTION
   Print to file the optical depth, cia, or extinction at the requested
   wavenumbers
   Return: 0 on success                                                     */
int
detailout(prop_samp *wn,         /* transit's wavenumber array              */
          prop_samp *rad,        /* Radius array                            */
          struct detailfld *det, /* Detail field struct                     */
          PREC_RES **arr,        /* Array of values to store                */
          short flag){           /* Flags                                   */

  long i,        /* Auxiliary for-loop index                                */
       u, d, m;  /* Auxiliary binary search indices                         */

  /* The radius index is first in array:                                    */
  _Bool radfirst = (_Bool)(flag & CIA_RADFIRST);
  /* Print as float value:                                                  */
  _Bool dofloat  = (_Bool)(flag & CIA_DOFLOAT);

  long idx[det->n];  /* Wavenumber indices                                  */
  double val;
  float **arrf = (float **)arr;       /* Float-casted array                 */
  FILE *out = fopen(det->file, "w");  /* Pointer to file                    */
  if(!out)
    transiterror(TERR_SERIOUS, "Cannot open '%s' for writing fine detail.\n",
                               det->file);
  transitprint(1, verblevel, "\nPrinting in '%s'. Fine detail of %s at "
                             "selected wavenumbers.\n", det->file, det->name);

  fprintf(out, "#Radius-w=>    ");
  /* Binary search to find the index for the requested wavenumbers:         */
  for(i=0; i < det->n; i++){
    val = det->ref[i];
    u = wn->n-1;
    if(val == wn->v[u])
      d = u;
    else{
      d = 0;
      while(u-d > 1){
        m = (u+d)/2;
        if(wn->v[m] > val)
          u = m;
        else
          d = m;
      }
    }
    idx[i] = d; /* Wavenumber index in transit array                        */
    /* Print the wavenumber:                                                */
    fprintf(out, "%-15.8g", wn->v[idx[i]]);
  }
  fprintf(out, "\n");

  /* Print radii and value:                                                 */
  if(radfirst){
    for(m=0; m<rad->n; m++){
      fprintf(out, "%-15.7g", rad->v[m]);
      for(i=0; i<det->n; i++){
        if(dofloat)
          val = arrf[m][idx[i]];
        else
          val = arr[m][idx[i]];
        fprintf(out, "%-15.7g", val);
      }
      fprintf(out, "\n");
    }
  }
  else{
    for(m=0; m<rad->n; m++){
      fprintf(out, "%-15.7g", rad->v[m]);
      for(i=0; i<det->n; i++){
        if(dofloat)
          val = arrf[idx[i]][m];
        else
          val = arr[idx[i]][m];
        fprintf(out, "%-15.7g", val);
      }
      fprintf(out, "\n");
    }
  }

  fclose(out);
  return 0;
}


/* FUNCTION
   Print (to file or stdout) for each wavelength the maximum optical depth
   calculated, the radius of such optical depth, and the index of the layer */
void
printtoomuch(char *file,            /* Filename to save the info            */
             struct optdepth *tau,  /* Tau information                      */
             prop_samp *wn,         /* Wavenumber sampling                  */
             prop_samp *rad){       /* Radius sampling                      */

  long w;              /* Auxiliary for-loop index for wavenumber           */
  FILE *out = stdout;  /* File pointer                                      */

  /* Open file if it was specified, default to stdout:                      */
  if(file[0] != '-')
    out = fopen(file, "w");
  if(!out)
    transiterror(TERR_WARNING, "Cannot open '%s' for writing depth where the "
                               "optical depth reached toomuch.\n", file);

  transitprint(20, verblevel, "\nPrinting in '%s' the depth where the "
                   "optical depth got larger than %g.\n", file, tau->toomuch);

  /* Print header:                                                          */
  fprintf(out,
    "# Wavelength   Max Optical   Radius at the    Radius\n"
    "   (microns)         depth   max depth (km)    index\n");
  /* Print the wavenumber and radius:                                       */
  for(w=0; w < wn->n; w++)
    fprintf(out, "%12.7f   %.5e     %12.4f     %04ld\n",
                 1.0/wn->v[w]*wn->fct*1e4, tau->t[w][tau->last[w]],
                 rad->v[tau->last[w]]*rad->fct/1e5,  tau->last[w]);
  fclose(out);
}


/* FUNCTION
   Free tau structure
   Return: 0 on success                                                     */
int
freemem_tau(struct optdepth *tau,  /* Optical depth structure               */
            long *pi){             /* Progress flag                         */

  /* Free arrays:                                                           */
  free(tau->t[0]);
  free(tau->t);
  free(tau->last);

  /* Update progress indicator and return:                                  */
  *pi &= ~(TRPI_TAU);
  return 0;
}

/* FUNCTION                                                                 */
void
outdebtauex(char *name,
            PREC_RES **e,
            prop_samp *ip,
            PREC_RES **t,
            long rn,
            long w){
  FILE *fp = fopen(name, "w");

  long j;
  for(j=0; j<rn; j++)
    fprintf(fp, "%-15.10g%-15.10g\t%-15.10g\n", ip->v[j], t[w][rn-j-1],
            e[j][w]);
  fclose(fp);
}


/* FUNCTION
   Print to file (name) the extinction coefficient as a function of radius
   (up to layer index rn) for the specified wavenumber range (indices from
    wi to wf).                                                              */
void
outdebex(char *name,    /* File name to save values                         */
         PREC_RES **e,  /* Extinction-coefficient array [nwn][nlayers]      */
         PREC_RES *r,   /* Radius array [nlayers]                           */
         long rn,       /* FINDME */
         long wi,       /* Initial wavenumber index                         */
         long wf){      /* Final wavenumber index                           */

  FILE *fp = fopen(name, "w");  /* File pointer                             */
  long i, j;                    /* Auxiliary for-loop indices               */

  for(j=0; j < rn; j++){   /* */
    fprintf(fp, "%-15.10g\t", r[j]);
    for(i=wi; i<=wf; i++)  /* */
      fprintf(fp, "%-15.10g\t", e[j][i]);
    fprintf(fp, "\n");
  }
  fclose(fp);
}


/* FUNCTION
   Print to file (name) the optical depth as function of impact parameter
   for the specified wavenumber range (indices from wi to wf).              */
void
outdebtau(char *name,  /* File name to save values                          */
       prop_samp *ip,  /* Impact parameter sampling [nlayers]               */
       PREC_RES **t,   /* Optical depth array  [nwn][nlayers]               */
       long wi,        /* Initial wavenumber index                          */
       long wf){       /* Final wavenumber index                            */

  FILE *fp = fopen(name, "w");  /* File pointer                             */
  long i, j;                    /* Auxiliary for-loop indices               */

  for(j=0; j < ip->n; j++){  /* Loop over the layers                        */
    fprintf(fp, "%-15.10g\t", ip->v[j]);
    for(i=wi; i <= wf; i++)  /* Loop over the wavenumbers                   */
      fprintf(fp, "%-15.10g\t", t[i][j]);
    fprintf(fp, "\n");
  }
  fclose(fp);  /* Close pointer                                             */
}
