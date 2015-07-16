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

/* FUNCTION:  Calculate the opacity due to molecular transitions.
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

  /* Check if the opacity file exists:                                      */
  int file_exists = fileexistopen(th->f_opa, NULL);
  transitprint(10, verblevel, "Opacity-file exist status = %d\n", file_exists);

  /* Set the file name in the transit struct:                               */
  tr->f_opa = th->f_opa;

  /* Opacity file has some error that makes it unusable, or not specified:  */
  if (file_exists < -1 || file_exists == 0) {

    /* Calculate Voigt profiles:                                            */
    transitprint(1, verblevel, "Calculating grid of Voigt profiles.\n");
    calcprofiles(tr);

    /* Set progress indicator and return success:                           */
    tr->pi |= TRPI_OPACITY;
    return 0;
  }

  /* Opacity file specified, but it just doesn't exist yet:                 */
  if (file_exists == -1) {

    /* Open file for writing:                                               */
    tr->fp_opa = fopen(tr->f_opa, "wb");

    /* Immediately return if the file could not be opened:                  */
    if (tr->fp_opa == NULL){
      transiterror(TERR_WARNING, "Opacity filename '%s' cannot be opened "
                                 "for writing.\n", tr->f_opa);
      return -1;
    }

    /* Calculate Voigt profiles:                                            */
    transitprint(1, verblevel, "Calculating grid of Voigt profiles.\n");
    calcprofiles(tr);

    /* Calculate the grid of opacities:                                     */
    transitprint(1, verblevel, "Calculating new grid of opacities: '%s'.\n",
                               tr->f_opa);
    calcopacity(tr, tr->fp_opa);

    /* Free the line-transition memory:                                     */
    freemem_linetransition(&tr->ds.li->lt, &tr->pi);
    tr->pi |= TRPI_READDATA;
    tr->pi |= TRPI_READINFO;
    tr->pi |= TRPI_CHKRNG;

    /* Set progress indicator and return success:                           */
    tr->pi |= TRPI_OPACITY;
    return 0;
  }

  /* Implied: The opacity file exists:                                      */
  file_exists = fileexistopen(th->f_opa, &tr->fp_opa);
  tr->f_opa = th->f_opa;
  transitprint(10, verblevel, "Opacity-file exist status = %d\n", file_exists);

  if (file_exists != 1) {

    transitprint(1, verblevel, "Opening opacity file failed.\n");

    /* Calculate Voigt profiles:                                            */
    transitprint(1, verblevel, "Calculating grid of Voigt profiles.\n");
    calcprofiles(tr);

    /* Set progress indicator and return success:                           */
    tr->pi |= TRPI_OPACITY;
    return 0;
  }

  /* Should attempt to use shared memory:                                   */
  if (tr->opashare) {

    /* Get ID or create shared opacityhint struct:                          */
    key_t hintkey = ftok(tr->f_opa, 'a');
    op.hintID = shmget(hintkey, sizeof(struct opacityhint), 0644 | IPC_CREAT);

    /* If reserving the hint segment was unsuccessful, give up:             */
    if (op.hintID == -1) {

      transitprint(1, verblevel, "Could not get hint shared memory ID.\n");

      /* Read the grid of opacities from file:                              */
      transitprint(1, verblevel, "Reading opacity file: '%s'.\n", tr->f_opa);
      readopacity(tr, tr->fp_opa);

      /* Set progress indicator and return success:                         */
      tr->pi |= TRPI_OPACITY;
      return 0;
    }

    /* Attach to the hint shared memory segment:                            */
    op.hint = shmat(op.hintID, (void *) 0, 0);

    /* If attaching was unsuccessful, give up.                              */
    if (op.hint == (struct opacityhint *) -1) {

      transitprint(1, verblevel, "Could not attach to hint shared memory.\n");

      /* Read the grid of opacities from file:                              */
      transitprint(1, verblevel, "Reading opacity file: '%s'.\n", tr->f_opa);
      readopacity(tr, tr->fp_opa);

      /* Set progress indicator and return success:                         */
      tr->pi |= TRPI_OPACITY;
      return 0;
    }

    /* If no process has claimed to be the master, claim it:                */
    if (op.hint->master_PID == 0)
      op.hint->master_PID = getpid();

    /* Wait a few cycles in case multiple processes claimed it at once:     */
    int i;
    for (i = 0; i < 10; i++);

    /* Check again to see if this process claimed the memory:               */
    if (op.hint->master_PID == getpid()) {

      /* Share the grid of opacities from file:                             */
      transitprint(1, verblevel, "Sharing opacity file: '%s'.\n", tr->f_opa);

      if (shareopacity(tr, tr->fp_opa) || mountopacity(tr)) {

        transitprint(1, verblevel, "Shared memory sh/mt error. Aborting.\n");
        op.hint->status |= TSHM_ERROR;

        /* Read the grid of opacities from file:                            */
        transitprint(1, verblevel, "Reading opacity file: '%s'.\n",
                                    tr->f_opa);
        readopacity(tr, tr->fp_opa);
      }

      /* Assuming that all processes have attached to the two shared
         memory segments, mark them for destruction as soon as all of the
         processes have detached.                                           */
      struct shmid_ds buf1, buf2;
      shmctl(op.hintID, IPC_STAT, &buf1);

      do {
        shmctl(op.mainID, IPC_STAT, &buf2);
      }
      while (buf2.shm_nattch != buf1.shm_nattch);

      transitprint(1, verblevel, "Marking shared memory for removal.\n");
      shmctl(op.hintID, IPC_RMID, &buf1);
      shmctl(op.mainID, IPC_RMID, &buf2);

      /* Set progress indicator and return success:                         */
      tr->pi |= TRPI_OPACITY;
      return 0;
    }

    else {

      while ((op.hint->status & TSHM_WRITTEN) == 0) {

        /* If there's an error with the shared memory, abort.               */
        if (op.hint->status & TSHM_ERROR) {

          transitprint(1, verblevel, "Shared memory error. Aborting.\n");

          /* Read the grid of opacities from file:                          */
          transitprint(1, verblevel, "Reading opacity file: '%s'.\n",
                                      tr->f_opa);
          readopacity(tr, tr->fp_opa);

          /* Set progress indicator and return success:                     */
          tr->pi |= TRPI_OPACITY;
          return 0;
        }
      }
    }

    /* If there's an error attaching or mounting the memory, abort.     */
    if (attachopacity(tr) || mountopacity(tr)) {

      transitprint(1, verblevel, "Shared memory att/mt error. Aborting.\n");

      /* Read the grid of opacities from file:                          */
      transitprint(1, verblevel, "Reading opacity file: '%s'.\n",
                                  tr->f_opa);
      readopacity(tr, tr->fp_opa);
    }
  }

  /* Should not attempt to use shared memory:                               */
  else {

    /* Read the grid of opacities from file:                                */
    transitprint(1, verblevel, "Reading opacity file: '%s'.\n", tr->f_opa);
    readopacity(tr, tr->fp_opa);
  }

  /* Set progress indicator and return success:                           */
  tr->pi |= TRPI_OPACITY;
  return 0;
}


/*  FUNCTION:  Calculate a grid of Voigt profiles.                          */
int
calcprofiles(struct transit *tr){
  struct transithint *th = tr->ds.th; /* transithint struct                 */
  struct opacity *op=tr->ds.op;     /* Opacity struct                       */
  int i, j;                         /* for-loop indices                     */
  int nDop, nLor;                   /* Number of Doppler and Lorentz-widths */
  double Lmin, Lmax, Dmin, Dmax;    /* Minimum and maximum widths           */
  PREC_VOIGT ***profile;            /* Grid of Voigt profiles               */
  float timesalpha=tr->timesalpha;  /* Voigt wings width                    */
  struct timeval tv;  /* Time-keeping variables                             */
  double t0=0.0;

  /* Make logscale grid for the profile widths:                             */
  /* FINDME: Add check that these numbers make sense                        */
  nDop = op->nDop = th->nDop;
  nLor = op->nLor = th->nLor;
  Dmin = th->dmin;
  Dmax = th->dmax;
  Lmin = th->lmin;
  Lmax = th->lmax;
  op->aDop = logspace(Dmin, Dmax, nDop);
  op->aLor = logspace(Lmin, Lmax, nLor);

  /* Allocate array for the profile half-size:                              */
  op->profsize    = (PREC_NREC **)calloc(nDop,      sizeof(PREC_NREC *));
  op->profsize[0] = (PREC_NREC  *)calloc(nDop*nLor, sizeof(PREC_NREC));
  for (i=1; i<nDop; i++)
    op->profsize[i] = op->profsize[0] + i*nLor;

  /* Allocate grid of Voigt profiles:                                       */
  op->profile    = (PREC_VOIGT ***)calloc(nDop,      sizeof(PREC_VOIGT **));
  op->profile[0] = (PREC_VOIGT  **)calloc(nDop*nLor, sizeof(PREC_VOIGT *));
  for (i=1; i<nDop; i++){
    op->profile[i] = op->profile[0] + i*nLor;
  }
  profile = op->profile;
  transitprint(10, verblevel, "Number of Voigt profiles: %d.\n", nDop*nLor);

  t0 = timestart(tv, "Begin Voigt profiles calculation.");
  /* Evaluate the profiles for the array of widths:                         */
  for   (i=0; i<nDop; i++){
    for (j=0; j<nLor; j++){
      /* Skip calculation if Doppler width << Lorentz width:                */
      /* Set size and pointer to previous profile:                          */
      if (op->aDop[i]*10.0 < op->aLor[j]  &&  i != 0){
        op->profsize[i][j] = op->profsize[i-1][j];
        profile[i][j] = profile[i-1][j];
      }
      else{ /* Calculate a new profile for given widths:                    */
        op->profsize[i][j] = getprofile(&profile[i][j],
                             tr->wns.d/tr->owns.o, op->aDop[i], op->aLor[j],
                             timesalpha, tr->owns.n);
      }
      transitprint(25, verblevel, "Profile[%2d][%2d] size = %4li  (D=%.3g, "
                                  "L=%.3g).\n", i, j, 2*op->profsize[i][j]+1,
                                   op->aDop[i], op->aLor[j]);
    }
  }
  t0 = timecheck(verblevel, 0, 0, "End Voigt-profile calculation.", tv, t0);
  return 0;
}

/* FUNCTION:  Calculate opacities for the grid of wavenumber, radius,
   and temperature arrays for each molecule.                                */
int
calcopacity(struct transit *tr,
            FILE *fp){
  struct opacity *op=tr->ds.op;     /* Opacity struct                       */
  struct isotopes  *iso=tr->ds.iso; /* Isotopes struct                      */
  struct molecules *mol=tr->ds.mol; /* Molecules struct                     */
  struct lineinfo *li=tr->ds.li;    /* Lineinfo struct                      */
  long Nmol, Ntemp, Nlayer, Nwave;  /* Opacity-grid  dimension sizes        */
  int i, j, t, r,                   /* for-loop indices                     */
      rn, iso1db;
  double *z;
  int k;

  PREC_ATM *density = (PREC_ATM *)calloc(mol->nmol, sizeof(PREC_ATM));
  double   *Z       = (double   *)calloc(iso->n_i,  sizeof(double));

  /* Make temperature array from hinted values:                             */
  maketempsample(tr);
  Ntemp = op->Ntemp = tr->temp.n;
  op->temp = (PREC_RES *)calloc(Ntemp, sizeof(PREC_RES));
  for (i=0; i<Ntemp; i++)
    op->temp[i] = tr->temp.v[i];
  /* Temperature boundaries check:                                          */
  if (op->temp[0] < li->tmin)
    transiterror(TERR_SERIOUS, "The opacity file attempted to sample a "
                 "temperature (%.1f K) below the lowest allowed "
                 "TLI temperature (%.1f K).\n", op->temp[0], li->tmin);
  if (op->temp[Ntemp-1] > li->tmax)
    transiterror(TERR_SERIOUS, "The opacity file attempted to sample a "
                 "temperature (%.1f K) beyond the highest allowed "
                 "TLI temperature (%.1f K).\n", op->temp[Ntemp-1], li->tmax);
  transitprint(1, verblevel, "There are %li temperature samples.\n", Ntemp);

  /* Evaluate the partition at these temperatures:                          */
  op->ziso    = (PREC_ATM **)calloc(iso->n_i,       sizeof(PREC_ATM *));
  op->ziso[0] = (PREC_ATM  *)calloc(iso->n_i*Ntemp, sizeof(PREC_ATM));
  for(i=1; i<iso->n_i; i++)
    op->ziso[i] = op->ziso[0] + i*Ntemp;

  /* Interpolate the partition function:                                    */
  for(i=0; i<iso->n_db; i++){  /* For each database separately:             */
    iso1db = iso->db[i].s;     /* Index of first isotope in current DB      */

    for(j=0; j < iso->db[i].i; j++){
      transitASSERT(iso1db + j > iso->n_i-1, "Trying to reference an isotope "
             "(%i) outside the extended limit (%i).\n", iso1db+j, iso->n_i-1);

      z = calloc(li->db[i].t, sizeof(double));
      spline_init(z, li->db[i].T, li->isov[iso1db+j].z, li->db[i].t);
      for(k=0;k<Ntemp;k++)
        op->ziso[iso1db+j][k] = splinterp_pt(z, li->db[i].t, li->db[i].T,
                                       li->isov[iso1db+j].z, op->temp[k]);
      free(z);
    }
  }

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
    /* If this molecule is not yet in molID array, add it's universal ID:   */
    if (valueinarray(op->molID, mol->ID[iso->imol[i]], j) < 0){
      op->molID[j++] = mol->ID[iso->imol[i]];
      transitprint(10, verblevel, "Isotope's (%d) molecule ID: %d (%s) "
                   "added at position %d.\n", i, op->molID[j-1],
                   mol->name[iso->imol[i]], j-1);
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
    op->o      = (PREC_RES ****)       calloc(Nlayer, sizeof(PREC_RES ***));
    for (r=0; r<Nlayer; r++){
      op->o[r] = (PREC_RES  ***)       calloc(Ntemp,  sizeof(PREC_RES **));
      for (t=0; t<Ntemp; t++){
        op->o[r][t] = (PREC_RES **)    calloc(Nmol,   sizeof(PREC_RES *));
        for (i=0; i<Nmol; i++){
          op->o[r][t][i] = (PREC_RES *)calloc(Nwave,  sizeof(PREC_RES));
        }
      }
    }

    if (!op->o[0][0][0])
      transitprint(1, verblevel, "Allocation fail.\n");

    /* Compute extinction:                                                  */
    for (r=0;   r<Nlayer; r++){  /* For each layer:                         */
      transitprint(3, verblevel, "\nOpacity Grid at layer %03d/%03ld.\n",
                                 r+1, Nlayer);
      for (t=0; t<Ntemp;  t++){  /* For each temperature:                   */
        /* Get density and partition-function arrays:                       */
        for (j=0; j < mol->nmol; j++)
          density[j] = stateeqnford(tr->ds.at->mass, mol->molec[j].q[r],
                       tr->atm.mm[r], mol->mass[j], op->press[r], op->temp[t]);
        for (j=0; j < iso->n_i; j++)
          Z[j] = op->ziso[j][t];
        if((rn=computemolext(tr, op->o[r][t], op->temp[t], density, Z, 1)) != 0)
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
    for (r=0; r<Nlayer; r++)
      for (t=0; t<Ntemp; t++)
        for (i=0; i<Nmol; i++)
          fwrite(op->o[r][t][i], sizeof(PREC_RES), Nwave, fp);

    fclose(fp);
  }
  transitprint(2, verblevel, "Done.\n");
  return 0;
}


/* FUNCTION: Read the opacity file and store values in the transit
   structure.                                                               */
int
readopacity(struct transit *tr,  /* transit struct                          */
            FILE *fp){           /* Pointer to file to read                 */
  struct opacity *op=tr->ds.op;  /* opacity struct                          */
  int i, t, r;  /* for-loop indices                                         */

  /* Read file dimension sizes:                                             */
  fread(&op->Nmol,   sizeof(long), 1, fp);
  fread(&op->Ntemp,  sizeof(long), 1, fp);
  fread(&op->Nlayer, sizeof(long), 1, fp);
  fread(&op->Nwave,  sizeof(long), 1, fp);
  transitprint(1, verblevel, "Opacity grid size: Nmolecules    = %5li\n"
                             "                   Ntemperatures = %5li\n"
                             "                   Nlayers       = %5li\n"
                             "                   Nwavenumbers  = %5li\n",
                             op->Nmol, op->Ntemp, op->Nlayer, op->Nwave);
  transitprint(20, verblevel, "ftell = %li\n", ftell(fp));

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
  transitprint(10, verblevel, "Molecule IDs = [");
  for (i=0; i < op->Nmol; i++)
    transitprint(10, verblevel, "%3d, ", op->molID[i]);
  transitprint(10, verblevel, "\b\b]\n\n");

  transitprint(10, verblevel, "Temperatures (K)  =  [");
  for (i=0; i < op->Ntemp; i++)
    transitprint(10, verblevel, "%6d, ", (int)op->temp[i]);
  transitprint(10, verblevel, "\b\b] \n\n");

  transitprint(10, verblevel, "Pressures (bar) = [ ");
  for (i=0; i < op->Nlayer; i++)
    transitprint(10, verblevel, "%.2e, ", 1e-6*op->press[i]);
  transitprint(10, verblevel, "\b\b] \n\n");

  transitprint(10, verblevel, "Wavenumber (cm-1) = [");
  for (i=0; i < 4; i++)
    transitprint(10, verblevel, "%7.2f, ", op->wns[i]);
  transitprint(10, verblevel, "..., ");
  for (i=op->Nwave-4; i < op->Nwave; i++)
    transitprint(10, verblevel, "%7.2f, ", op->wns[i]);
  transitprint(10, verblevel, "\b\b]\n\n");

  /* Allocate and read the opacity grid:                                    */
  op->o      = (PREC_RES ****)       calloc(op->Nlayer, sizeof(PREC_RES ***));
  for     (r=0; r < op->Nlayer; r++){
    op->o[r] = (PREC_RES  ***)       calloc(op->Ntemp,  sizeof(PREC_RES **));
    for   (t=0; t < op->Ntemp; t++){
      op->o[r][t] = (PREC_RES **)    calloc(op->Nmol,   sizeof(PREC_RES *));
      for (i=0; i < op->Nmol; i++){
        op->o[r][t][i] = (PREC_RES *)calloc(op->Nwave,  sizeof(PREC_RES));
      }
    }
  }

  /* Read the opacity grid:                                                 */
  for     (r=0; r < op->Nlayer; r++)
    for   (t=0; t < op->Ntemp;  t++)
      for (i=0; i < op->Nmol;   i++)
        fread(op->o[r][t][i], sizeof(PREC_RES), op->Nwave, fp);

  return 0;
}


/* FUNCTION: Read the opacity file and store values in shared memory.       */
int
shareopacity(struct transit *tr, /* transit struct                          */
            FILE *fp){           /* Pointer to file to read                 */
  struct opacity *op=tr->ds.op;  /* opacity struct                          */
  struct opacityhint *oh=op->hint;  /* opacity hint struct                  */

  /* Read file dimension sizes:                                             */
  fread(&op->Nmol,   sizeof(long), 1, fp);
  fread(&op->Ntemp,  sizeof(long), 1, fp);
  fread(&op->Nlayer, sizeof(long), 1, fp);
  fread(&op->Nwave,  sizeof(long), 1, fp);
  transitprint(1, verblevel, "Opacity grid size: Nmolecules    = %5li\n"
                             "                   Ntemperatures = %5li\n"
                             "                   Nlayers       = %5li\n"
                             "                   Nwavenumbers  = %5li\n",
                             op->Nmol, op->Ntemp, op->Nlayer, op->Nwave);
  transitprint(20, verblevel, "ftell = %li\n", ftell(fp));

  /* Copy dimensional data into the shared hint struct:                     */
  oh->Nwave = op->Nwave;
  oh->Ntemp = op->Ntemp;
  oh->Nlayer = op->Nlayer;
  oh->Nmol = op->Nmol;

  /* If creating and attaching the main segment fails, return:              */
  if (attachopacity(tr))
    return 1;

  /* Read arrays:                                                           */
  char *p = op->mainaddr;
  fread(p, sizeof(int),      op->Nmol,   fp);
  p += sizeof(int) * op->Nmol;
  fread(p,  sizeof(PREC_RES), op->Ntemp,  fp);
  p += sizeof(PREC_RES) * op->Ntemp;
  fread(p, sizeof(PREC_RES), op->Nlayer, fp);
  p += sizeof(PREC_RES) * op->Nlayer;
  fread(p,   sizeof(PREC_RES), op->Nwave,  fp);
  p += sizeof(PREC_RES) * op->Nwave;

  /* Read opacity grid:                                                     */
  fread(p, sizeof(PREC_RES),
        op->Nmol * op->Ntemp * op->Nlayer * op->Nwave, fp);

  oh->status |= TSHM_WRITTEN;
  return 0;
}


/* FUNCTION: Locate (or create) and attach to the main shared memory.
   Return: 0 on success, 1 on failure                                       */
int
attachopacity(struct transit *tr){ /* transit struct                        */
  struct opacity *op=tr->ds.op;   /* opacity struct                         */
  struct opacityhint *oh=op->hint;  /* opacity hint struct                  */

  op->Nwave = oh->Nwave;
  op->Ntemp = oh->Ntemp;
  op->Nlayer = oh->Nlayer;
  op->Nmol = oh->Nmol;

  int main_shm_size  /* Size of the main shared memory, starting with: grid */
    = sizeof(PREC_RES) * op->Nmol * op->Ntemp * op->Nlayer * op->Nwave
    + sizeof(int) * op->Nmol          /* op->molID  */
    + sizeof(PREC_RES) * op->Ntemp    /* op->temp   */
    + sizeof(PREC_RES) * op->Nlayer   /* op->press  */
    + sizeof(PREC_RES) * op->Nwave;   /* op->wns    */

  /* Allocate or locate the main shared memory:                             */
  key_t mainkey = ftok(tr->f_opa, 'b');
  op->mainID = shmget(mainkey, main_shm_size, 0644 | IPC_CREAT);

  /* If allocation failed, abort:                                           */
  if (op->mainID == -1) {
    transiterror(TERR_WARNING, "Failed to allocate main shared memory.\n");
    transitprint(1, verblevel, "Attachment failed. shmget() returned -1.\n");
    return 1;
  }

  /* Attach to the main shared memory:                                      */
  op->mainaddr = shmat(op->mainID, (void *) 0, 0);

  if (op->mainaddr == (void *) -1) {
    transiterror(TERR_WARNING, "Failed to attach to main shared memory.\n");
    transitprint(1, verblevel, "Attachment failed. shmat() returned -1.\n");
    return 1;
  }

  oh->num_attached++;
  return 0;
}


/* FUNCTION: Create pointers to the opacity values in shared memory.
   Return: 0 on success                                                     */
int
mountopacity(struct transit *tr){ /* transit struct                         */
  struct opacity *op=tr->ds.op;   /* opacity struct                         */
  int i, t, r;  /* for-loop indices                                         */

  char *p = op->mainaddr;
  op->molID = (int *) p;
  p += sizeof(int) * op->Nmol;
  op->temp = (PREC_RES *) p;
  p += sizeof(PREC_RES) * op->Ntemp;
  op->press = (PREC_RES *) p;
  p += sizeof(PREC_RES) * op->Nlayer;
  op->wns = (PREC_RES *) p;
  p += sizeof(PREC_RES) * op->Nwave;

  op->o      = (PREC_RES ****)       calloc(op->Nlayer, sizeof(PREC_RES ***));
  for     (r=0; r < op->Nlayer; r++){
    op->o[r] = (PREC_RES  ***)       calloc(op->Ntemp,  sizeof(PREC_RES **));
    for   (t=0; t < op->Ntemp; t++){
      op->o[r][t] = (PREC_RES **)    calloc(op->Nmol,   sizeof(PREC_RES *));
      for (i=0; i < op->Nmol; i++){
        op->o[r][t][i] = (PREC_RES *) p /* Map the 4D structure to 1D.      */
                                      + r * op->Ntemp * op->Nmol * op->Nwave
                                      + t * op->Nmol * op->Nwave
                                      + i * op->Nwave;
      }
    }
  }

  return 0;
}


/* FUNCTION: Free opacity array
   Return: 0 on success                                                     */
int
freemem_opacity(struct opacity *op, /* Opacity structure                    */
                long *pi){          /* transit progress flag                */
  /* Free arrays:                                                           */
  free(op->o[0][0][0]); /* The opacity                                      */
  free(op->o[0][0]);
  free(op->o[0]);
  free(op->o);

  free(op->profile[0][0]); /* The Voigt profiles                            */
  free(op->profile[0]);
  free(op->profile);

  free(op->profsize[0]);  /* The Voigt-profile half-size                    */
  free(op->profsize);

  /* Update progress indicator and return:                                  */
  *pi &= ~(TRPI_OPACITY | TRPI_TAU);

  return 0;
}
