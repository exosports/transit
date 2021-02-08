// Copyright (C) 2015-2016 University of Central Florida. All rights reserved.
// Transit is under an open-source, reproducible-research license (see LICENSE).

#include <transit.h>

#define TLI_WAV_UNITS 1e-4 /* TLI wavelength (microns, as of v4) */
#define TLI_E_UNITS   1    /* TLI Elow units (cm-1, as of v4)    */

static double tli_to_microns = TLI_WAV_UNITS/1e-4;


/* FNUCTION:
  Do a binary search in file pointed by 'fp' between 'off' and 'off+nfields'
  looking for 'target' as the first item of a record of length 'reclength',
  result index (with respect to offs) is stored in 'resultp'.               */
void
datafileBS(FILE *fp,            /* File pointer                             */
           PREC_NREC offs,      /* Initial position of data in tli file     */
           PREC_NREC nfields,   /* Number of fields to search               */
           PREC_LNDATA target,  /* Target value                             */
           PREC_NREC *resultp,  /* Result index                             */
           int reclength,       /* Total length of record                   */
           int up){             /* Flag to search up, or down               */

  /* Variable to keep wavelength:                                           */
  PREC_LNDATA temp;
  /* Search boundaries:                                                     */
  PREC_NREC lo=0,          /* Starting point of binary search               */
            hi=nfields-1,  /* Starting point of binary search               */
            loc;           /* Current location of closest value             */

  tr_output(TOUT_DEBUG, "BS: Start looking from %li in %li fields "
    "for %f\n", offs, nfields, target);
  /* Binary search:                                                         */
  do{
    loc = (hi+lo)/2;                           /* Mid record's index        */
    fseek(fp, offs+reclength*loc, SEEK_SET);   /* Move pointer              */
    fread(&temp, sizeof(PREC_LNDATA), 1, fp);  /* Read value                */
    tr_output(TOUT_DEBUG, "BS: found wl %.8f microns at position %li\n",
                                temp*tli_to_microns, loc);
    /* Re-set lower or higher boundary:                                     */
    if(target > temp)
      lo = loc;
    else
      hi = loc;
  }while (hi-lo > 1);

  /* Linear search:                                                         */
  if (up){
    loc = lo;
    /* Linear search for the entries above loc:                             */
    while(loc < nfields-1){
      fseek(fp, offs + reclength*(loc+1), SEEK_SET);
      fread(&temp, sizeof(PREC_LNDATA), 1, fp);
      if (temp > target)
        break;
      loc++;
    }
  }else{
    loc = hi;
    /* Linear search for the entries below loc:                             */
    while(loc > 0){
      fseek(fp, offs + reclength*(loc-1), SEEK_SET);
      fread(&temp, sizeof(PREC_LNDATA), 1, fp);
      if (temp < target)
        break;
      loc--;
    }

  }
  /* Final remarks:                                                         */
  *resultp = loc;
  fseek(fp, offs + reclength*loc, SEEK_SET);
  fread(&temp, sizeof(PREC_LNDATA), 1, fp);
  tr_output(TOUT_RESULT, "Binary search found wavelength: %.8f at "
    "position %li.\n", temp*tli_to_microns, loc);
}


/* FUNCTION:
  Read initial and final wavelength limits and number of databases.
  Allocate pointers to databases, and isotopes arrays.  Get databases
  info: names, number of temperatures, temperatures, number of
  isotopes, isotope names and masses, partition function, and cross
  sections. Get cumulative number of isotopes.
  Returns 0 on success                                                      */
int
readtli_bin(FILE *fp,
            struct transit *tr,
            struct lineinfo *li){
  /* Declare varables:                                                      */
  double iniw, finw;     /* Initial and final wavelength of database        */
  unsigned short ndb;    /* Number of databases                             */
  unsigned short rs;     /* FINDME: auxiliary what??                        */
  unsigned short nT,     /* Number of temperatures per database             */
                 niso=0, /* Cumulative number of isotopes per database      */
                 nDBiso; /* Number of isotopes per database                 */
  PREC_ZREC *T, *Z;      /* Auxiliary temperature and part. func. pointers  */
  int correliso=0;       /* Isotopes correlative number                     */
  int i, j;
  struct isotopes *iso=tr->ds.iso;

  /* Read TLI version, lineread version, and lineread revision number:      */
  fread(&li->tli_ver, sizeof(unsigned short), 1, fp);
  fread(&li->lr_ver,  sizeof(unsigned short), 1, fp);
  fread(&li->lr_rev,  sizeof(unsigned short), 1, fp);
  /* Check compatibility of versions:                                       */
  if(li->tli_ver != compattliversion) {
    tr_output(TOUT_ERROR,
      "The version of the TLI file: %i (lineread v%i.%i) is not "
      "compatible with this version of transit, which can only "
      "read version %i.\n", li->tli_ver, li->lr_ver,
      li->lr_rev, compattliversion);
    exit(EXIT_FAILURE);
  }

  /* Read initial wavelength, final wavelength, and number of databases:    */
  fread(&iniw, sizeof(double), 1, fp);
  fread(&finw, sizeof(double), 1, fp);
  tr_output(TOUT_RESULT, "Initial wavelength: %.2f (um)\n"
    "Final   wavelength: %.2f (um)\n", iniw, finw);
  fread(&ndb, sizeof(unsigned short), 1, fp);
  tr_output(TOUT_RESULT, "Number of databases: %d.\n", ndb);

  /* Allocate pointers according to the number of databases:                */
  iso->db       = (prop_db      *)calloc(ndb, sizeof(prop_db     ));
  li->db        = (prop_dbnoext *)calloc(ndb, sizeof(prop_dbnoext));
  /* Allocate isotope info.  Start with size 1, then reallocate as needed   */
  iso->isof     = (prop_isof    *)calloc(1,   sizeof(prop_isof));
  li->isov      = (prop_isov    *)calloc(1,   sizeof(prop_isov));
  iso->isoratio = (double       *)calloc(1,   sizeof(double));

  /* Read info for each database:                                           */
  for(i=0; i<ndb; i++){
    /* Allocate and get DB's name:                                          */
    fread(&rs, sizeof(unsigned short), 1, fp);         /* Get DBname length */
    iso->db[i].n = (char *)calloc(rs+1, sizeof(char)); /* Allocate          */
    fread(iso->db[i].n, sizeof(char), rs, fp);         /* Read              */
    iso->db[i].n[rs] = '\0';
    tr_output(TOUT_DEBUG, "  DB name size: %d'\n", rs);
    /* Read and allocate the molecule's name:                               */
    fread(&rs, sizeof(unsigned short), 1, fp);
    iso->db[i].molname = (char *)calloc(rs+1, sizeof(char));
    fread(iso->db[i].molname, sizeof(char), rs, fp);
    iso->db[i].molname[rs] = '\0';

    tr_output(TOUT_RESULT, "Database (%d/%d) name: '%s' (%s molecule)\n",
                                i+1, ndb, iso->db[i].n, iso->db[i].molname);

    /* Get number of temperatures and isotopes:                             */
    fread(&nT,     sizeof(unsigned short), 1, fp);
    fread(&nDBiso, sizeof(unsigned short), 1, fp);
    tr_output(TOUT_RESULT, "  Number of temperatures: %d\n"
                               "  Number of isotopes:     %d\n", nT, nDBiso);
    li->db[i].t  = nT;
    iso->db[i].i = nDBiso;

    /* Allocate for the different temperature points and read:              */
    T = li->db[i].T = (double *)calloc(nT, sizeof(double));
    fread(T, sizeof(double), nT, fp);
    tr_output(TOUT_RESULT, "  Temperatures: [%6.1f, %6.1f, ..., %6.1f]\n",
                               T[0], T[1], T[nT-1]);
    li->tmin = fmax(li->tmin, T[   0]);
    li->tmax = fmin(li->tmax, T[nT-1]);

    /* Reallocate memory to account for new isotopes:                       */
    li->isov  = (prop_isov  *)realloc(li->isov,
                                      (correliso+nDBiso)*sizeof(prop_isov));
    iso->isof = (prop_isof  *)realloc(iso->isof,
                                      (correliso+nDBiso)*sizeof(prop_isof));
    iso->isoratio = (double *)realloc(iso->isoratio,
                                      (correliso+nDBiso)*sizeof(double));
    /* Allocate memory for this DB's partition function:                    */
    li->isov[correliso].z = (double *)calloc((correliso+nDBiso)*nT,
                                             sizeof(double));

    tr_output(TOUT_DEBUG, "So far, cumIsotopes: %i, at databases: %i, "
      "position %li.\n", correliso+nDBiso, i, ftell(fp));

    /* Display database general info:                                       */
    tr_output(TOUT_DEBUG, "DB %i: '%s' has %i (%i) temperatures, "
      "%i (%i) isotopes, and starts at cumulative isotope %i.\n",
      iso->isof[correliso].d, iso->db[i].n,
      li->db[i].t, nT, iso->db[i].i, nDBiso, iso->db[i].s);


    /* Read the isotopes from this data base:                               */
    for (j=0; j<nDBiso; j++){
      /* Store isotopes'  DB index number:                                  */
      iso->isof[correliso].d = i;
      tr_output(TOUT_DEBUG,"  Correlative isotope number: %d", correliso);

      /* Read isotopes' name:                                               */
      fread(&rs, sizeof(unsigned short), 1, fp);
      iso->isof[correliso].n = (char *)calloc(rs+1, sizeof(char));
      fread(iso->isof[correliso].n, sizeof(char), rs, fp);
      iso->isof[correliso].n[rs] = '\0';
      tr_output(TOUT_RESULT, "  Isotope (%i/%i): '%s'\n", j+1, nDBiso,
                                  iso->isof[correliso].n);
      tr_output(TOUT_DEBUG, "   Isotope name size: %d'\n", rs);

      /* Read mass and isotopic ratio:                                      */
      fread(&iso->isof[correliso].m,   sizeof(double), 1, fp);
      fread((iso->isoratio+correliso), sizeof(double), 1, fp);
      tr_output(TOUT_RESULT, "    Mass:  %g u (%g gr)\n",
                          iso->isof[correliso].m, iso->isof[correliso].m*AMU);
      tr_output(TOUT_RESULT, "    Isotopic ratio: %.4g\n",
                          iso->isoratio[correliso]);
      tr_output(TOUT_DEBUG, "    File position: %li.\n", ftell(fp));

      /* Read partition function:                                           */
      Z  = li->isov[correliso].z = li->isov[correliso-j].z + nT*j;
      fread(Z,  sizeof(double), nT, fp);
      li->isov[correliso].n = nT;

      tr_output(TOUT_DEBUG, "    Part Function:    [%.2e, %.2e, ..., "
                                  "%.2e]\n", Z[0],  Z[1],  Z[nT-1]);
      correliso++;
    }

    /* Update cumulative isotope count (index of first isitope in this DB): */
    iso->db[i].s = niso;
    niso += nDBiso;
  }

  tr_output(TOUT_RESULT, "Cumulative number of isotopes per DB: [");
  for (i=0; i<ndb; i++)
    tr_output(TOUT_RESULT,"%2d, ", iso->db[i].s);
  tr_output(TOUT_RESULT, "\b\b].\n");
  tr_output(TOUT_RESULT, "acum Iso: %2d.\n", niso);

  /* Update structure values:                                               */
  li->ni  = iso->n_i  = niso;  /* Number of isotopes                        */
  li->ndb = iso->n_db = ndb;   /* Number of databases                       */
  /* Position of first line data                                            */
  li->endinfo = ftell(fp);
  li->wi = iniw;            /* Initial wavelength                           */
  li->wf = finw;            /* Final wavelength                             */
  /* Allocate isotope's variable data                                       */
  iso->isov = (prop_isov *)calloc(iso->n_i, sizeof(prop_isov));

  tr->pi |= TRPI_READBIN;
  return 0;
}


/* FUCTION:
   Set the molecule's index of the isotopes:                                */
int
setimol(struct transit *tr){
  struct molecules *mol = tr->ds.mol;
  struct isotopes  *iso = tr->ds.iso;
  int i; /* Auxiliary for-loop index                                        */

  /* Return if there are no isotopes (i.e., no line transitions):           */
  if (iso->n_i == 0)
    return 0;

  iso->imol = (int *)calloc(iso->n_i, sizeof(int));
  iso->nmol = 0;
  for (i=0; i<iso->n_i; i++){
    //transitprint(1, verblevel, "Isotope %d is '%s', from DB %d: '%s'.\n",
    //                  i, iso->isof->n, iso->isof->d, iso->db[iso->isof->d].n);
    /* Search water molecule's index:                                       */
    iso->imol[i] = findstring(iso->db[iso->isof[i].d].molname, mol->name,
                              mol->nmol);
    tr_output(TOUT_DEBUG, "Isotope '%s', is mol %d: '%s'.\n",
                 iso->isof[i].n, iso->imol[i], mol->name[iso->imol[i]]);

    /* Find if current imol is already in imol array:                       */
    if (valueinarray(iso->imol, iso->imol[i], i) < 0){
      tr_output(TOUT_DEBUG, "Isotope %s (%d) is a new molecule (%s).\n",
                                  iso->isof[i].n, i, mol->name[iso->imol[i]]);
      iso->nmol++;
    }
  }
  return 0;
}


/* FUNCTION
   Initialize wavelength sample struct.
   Set initial and final wavelengths to use.

   Return:  0   if all hinted values were accepted, else
           (positive if something was changed):
           0x20 if hinted final wavelength was changed
           0x10 if hinted initial wavelength was changed
           (negative: If something bad happen):
           -1   if hinted initial wavelength is larger than sug. final
           -2   if hinted initial is larger than largest allowed final
           -3   if hinted final is shorter than smallest allowed initial    */
int
checkrange(struct transit *tr,   /* General parameters and  hints           */
           struct lineinfo *li){ /* Values returned by readinfo_tli         */

  int res=0;                                /* Return value                 */
  prop_samp *tsamp = &tr->wns;              /* Transit wavenumber sampling  */
  PREC_LNDATA dbini = li->wi*TLI_WAV_UNITS, /* Minimum DB wavelength        */
              dbfin = li->wf*TLI_WAV_UNITS; /* Maximum DB wavelength        */
  double cm_to_micron = 1e4;                /* Conversion factor to microns */
  double wlmin, wlmax;

  /* Transit wavelength limits in cgs units:                                */
  wlmin = 1.0/(tsamp->f*tsamp->fct);
  wlmax = 1.0/(tsamp->i*tsamp->fct);

  tr_output(TOUT_DEBUG,
    "Transit initial and final wavelengths are %6g and %6g cm.\n"
    "The database max and min wavelengths are  %6g and %6g cm.\n",
    wlmin, wlmax, dbini, dbfin);

  /* Check that it is not below the minimum value:                          */
  if(dbini > wlmax){
    tr_output(TOUT_ERROR,
      "Final wavelength (%g) is smaller than minimum wavelength "
      "in database (%g).\n", wlmax, dbini);
    return -3;
  }
  /* Warn if it is above the maximum TLI value:                             */
  if(wlmax > dbfin)
    tr_output(TOUT_WARN, "Final requested wavelength (%g microns) "
                 "is larger than the maximum informative value in database "
                 "(%g microns).\n", wlmax*cm_to_micron, dbfin*cm_to_micron);

  /* Check that it is not larger than the maximum db wavelength:            */
  if(dbfin < wlmin){
    tr_output(TOUT_ERROR, "Initial wavelength (%g cm) is larger than "
      "maximum wavelength in database (%g cm).\n", wlmin, dbfin);
    return -2;
  }
  /* Warn if it is below the maximum TLI value:                             */
  if(wlmin < dbini)
    tr_output(TOUT_WARN, "Initial requested wavelength (%g microns) "
                 "is smaller than the minimum informative value in database "
                 "(%g microns).\n", wlmin*cm_to_micron, dbini*cm_to_micron);

  /* Return status:                                                         */
  return res;
}


/* FUNCTION:
   Check TLI file exists.  Check that machine formating is compatible
   with lineread.  Read either TLI file.  Declare line_transition.

  Return: 1 on success
         -1 unavailable file
         -2 Filename not hinted                                             */
int
readinfo_tli(struct transit *tr,
             struct lineinfo *li){
  struct transithint *th = tr->ds.th;  /* Pointer to hint:                  */
  int rn;
  FILE *fp;  /* File pointer of info file:                                  */

  /* Decalre and initialize the union sign:                                 */
  /* sign.s contains the magic numbers of this machine's and TLI's:         */
  union {char sig[4]; int32_t s[2];} sign =
    {.s={0, ((0xff-'T')<<24)|((0xff-'L')<<16)|((0xff-'I')<<8)|(0xff)}};

  /* Get TLI file name from hint:                                           */
  if(!th->f_line){  /* Check if it was defined in hint                      */
    tr_output(TOUT_WARN, "No TLI file set.\n");
    tr->pi |= TRPI_READINFO;
    return -2;
  }
  /* Attempt to open the TLI file and make a pointer to it:                 */
  if((rn=fileexistopen(th->f_line, &fp)) != 1){
    tr_output(TOUT_ERROR,
      "Line info file '%s' is not found. "
      "fileexistopen() error code %i.\n", th->f_line, rn);
    return -1;
  }
  /* Set transit TLI file pointer and TLI file name:                        */
  tr->fp_line = fp;
  tr->f_line  = th->f_line;

  /* Read first four bytes, they should be either
  `(0xff-T)(0xff-L)(0xff-I)(0xff)'.  They are stored as integer.
  This checks whether the machine where the TLI file and the one this
  program is being run have the same endian order.                          */
  fread(sign.s, sizeof(int32_t), 1, fp);

  /* Read binary TLI:                                                       */
  if((rn=readtli_bin(fp, tr, li)) != 0){
    tr_output(TOUT_ERROR, "readtli_bin() return error code %i.\n", rn);
    return -6;
  }
  tr_output(TOUT_RESULT, "TLI file read from %g to %g microns.\n",
                             li->wi, li->wf);

  /* Declare linetransition struct and set wavelength and lower energy unit
     factors (As of TLI v5, always in microns and cm-1, respectively):      */
  struct line_transition *lt = &li->lt;
  lt->wfct = TLI_WAV_UNITS;
  lt->efct = TLI_E_UNITS;

  /* Close TLI file pointer, set progres indicator and return success:      */
  fclose(fp);
  tr->pi |= TRPI_READINFO;
  return 1;
}


/* FUNCTION:
  Read and store the line transition info (central wavelength, isotope
  ID, lowE, log(gf)) into lineinfo.  Return the number of lines read.

  Return: the number of records read on success, else:
          -1 unexpected EOF
          -2 file non-seekable
          -3 on non-integer number of structure records
          -4 First field is not valid while looking for starting point
          -5 One of the fields contained an invalid flaoating point         */
int readdatarng(struct transit *tr,   /* transit structure                  */
                struct lineinfo *li){ /* lineinfo structure                 */

  struct line_transition *lt = &li->lt;  /* line_transition structure       */
  FILE *fp;              /* Data file pointer                               */
  int niso,              /* Number of isotopes in line transition data      */
      nread,             /* Number of transitions to read for each isotope  */
      i;                 /* for-loop index                                  */
  PREC_NREC nlines,            /* Number of line transitions                */
            wl_loc, iso_loc,   /* Offsets for isoID, Elow, and gf data      */
            el_loc, gf_loc,    /* (in memory)                               */
            *isotran,          /* Number of transitions per isotope in TLI  */
            start,             /* Position of first LT for isotope in TLI   */
            offset=0;          /* Isotope offset (in number of transitions) */
  int rn;                /* Return IDs                                      */
  /* Indices of first and last transitions to be stored                     */
  PREC_NREC ifirst, ilast;

  /* Auxiliary variables to keep wavelength limits:                         */
  PREC_LNDATA iniw = 1.0/(tr->wns.f*tr->wns.fct) / TLI_WAV_UNITS;
  PREC_LNDATA finw = 1.0/(tr->wns.i*tr->wns.fct) / TLI_WAV_UNITS;

  /* Open line data file:                                                   */
  rn = fileexistopen(tr->f_line, &fp);
  if (rn == 0){  /* No file given */
    li->n_l = 0;
    return 0;
  }

  else if (rn != 1){
    tr_output(TOUT_ERROR, "Data file '%s' not found. "
      "fileexistopen() error code: %i.\n", tr->f_line, rn);
    return -1;
  }

  /* Check seekability:                                                     */
  if(fseek(fp, 0, SEEK_CUR)){
    tr_output(TOUT_ERROR, "File '%s' was not seekable.\n", tr->f_line);
    return -2;
  }

  /* Read total number of transitions in TLI file:                          */
  fseek(fp, li->endinfo, SEEK_SET);
  fread(&nlines, sizeof(PREC_NREC), 1, fp);
  tr_output(TOUT_INFO, "TLI has %d transition lines.\n", nlines);

  /* Number of isotopes in line transition data:                            */
  fread(&niso, sizeof(int), 1, fp);
  tr_output(TOUT_DEBUG, "TLI has %d isotopes.\n", niso);
  /* Number of transitions per isotope:                                     */
  isotran = calloc(niso, sizeof(PREC_NREC));
  fread(isotran, sizeof(PREC_NREC), niso, fp);
  for (i=0; i<niso; i++){
    tr_output(TOUT_DEBUG, "Ntransitions[%d]: %d.\n", i, isotran[i]);
  }

  /* Get current location of pointer:                                       */
  start = ftell(fp);
  li->n_l = 0;

  /* Allocation for line transition structures:                             */
  /* The size might be larger than needed, adjust at the end                */
  lt->gf    = (PREC_LNDATA *)calloc(nlines, sizeof(PREC_LNDATA));
  lt->wl    = (PREC_LNDATA *)calloc(nlines, sizeof(PREC_LNDATA));
  lt->elow  = (PREC_LNDATA *)calloc(nlines, sizeof(PREC_LNDATA));
  lt->isoid = (short       *)calloc(nlines, sizeof(short));
  /* Check for allocation errors:                                           */
  if(!lt->gf || !lt->wl || !lt->elow || !lt->isoid) {
    tr_output(TOUT_ERROR, "Couldn't allocate memory for "
      "linetran structure array of length %i, in function "
      "readdatarng.\n", nlines);
    exit(EXIT_FAILURE);
  }

  /* Starting location for wavelength, isoID, Elow, and gf data in file:    */
  wl_loc  = start;
  iso_loc = wl_loc  + nlines*sizeof(PREC_LNDATA);
  el_loc  = iso_loc + nlines*sizeof(short);
  gf_loc  = el_loc  + nlines*sizeof(PREC_LNDATA);

  for (i=0; i<niso; i++){
    /* Do binary search in units of TLI:                                    */
    datafileBS(fp, start, isotran[i], iniw, &ifirst, sizeof(PREC_LNDATA), 0);
    datafileBS(fp, start, isotran[i], finw, &ilast,  sizeof(PREC_LNDATA), 1);
    ifirst += offset;
    ilast  += offset;
    tr_output(TOUT_DEBUG, "Initial and final entries are: "
      "%li and %li.\n", ifirst, ilast);

    /* Number of transitions to read:                                       */
    nread = ilast - ifirst + 1;
    /* Move pointer to each section and read info:                          */
    /* Wavelength:                                                          */
    fseek(fp, ifirst*sizeof(PREC_LNDATA) + wl_loc,  SEEK_SET);
    fread(lt->wl+li->n_l,    sizeof(PREC_LNDATA), nread, fp);
    /* Isotope ID:                                                          */
    fseek(fp, ifirst*sizeof(short)       + iso_loc, SEEK_SET);
    fread(lt->isoid+li->n_l, sizeof(short),       nread, fp);
    /* Lower-state energy:                                                  */
    fseek(fp, ifirst*sizeof(PREC_LNDATA) + el_loc,  SEEK_SET);
    fread(lt->elow+li->n_l,  sizeof(PREC_LNDATA), nread, fp);
    /* gf:                                                                  */
    fseek(fp, ifirst*sizeof(PREC_LNDATA) + gf_loc,  SEEK_SET);
    fread(lt->gf+li->n_l,    sizeof(PREC_LNDATA), nread, fp);

    /* Count the number of lines:                                           */
    li->n_l += nread;
    /* Move the wl offset to next isotope:                                  */
    start  += isotran[i]*sizeof(double);
    offset += isotran[i];
  }

  /* Re-allocate arrays to their correct size:                              */
  lt->wl    = (PREC_LNDATA *)realloc(lt->wl,    li->n_l*sizeof(PREC_LNDATA));
  lt->isoid = (short       *)realloc(lt->isoid, li->n_l*sizeof(short));
  lt->elow  = (PREC_LNDATA *)realloc(lt->elow,  li->n_l*sizeof(PREC_LNDATA));
  lt->gf    = (PREC_LNDATA *)realloc(lt->gf,    li->n_l*sizeof(PREC_LNDATA));

  fclose(fp);               /* Close file                                   */
  tr->pi |= TRPI_READDATA;  /* Update progress indicator                    */
  return li->n_l;           /* Return the number of lines read              */
}


/* FUNCTION:
    Driver function to read TLI: read isotopes info, check
    and ranges, and read line transition information.
    Return: 0 on success.                                                   */
int
readlineinfo(struct transit *tr){
  struct transithint *th=tr->ds.th;
  static struct lineinfo li;
  static struct isotopes iso;
  long rn;  /* Sub-routines returned status */
  int filecheck;  /* Integer to check if opacity file exists */

  memset(&li,  0, sizeof(struct lineinfo));
  memset(&iso, 0, sizeof(struct isotopes));
  tr->ds.li  = &li;   /* lineinfo                                           */
  tr->ds.iso = &iso;  /* isotopes                                           */

  /* Set some defaults:                                                     */
  li.ni  = iso.n_i  = 0;  /* Number of isotopes                             */
  li.ndb = iso.n_db = 0;  /* Number of databases                            */
  /* Min and max allowed temperatures in TLI files:                         */
  li.tmin =     0.0;
  li.tmax = 70000.0;

  /* Read hinted info file:                                                 */
  tr_output(TOUT_INFO, "Reading info file '%s' ...\n", th->f_line);
  rn = readinfo_tli(tr, &li);
  tr_output(TOUT_INFO, "Done.\n\n");

  /* Check the remainder range of the hinted values
     related to line database reading:                                      */
  if (rn != -2){
    if((rn=checkrange(tr, &li)) < 0) {
      tr_output(TOUT_ERROR, "checkrange() returned error code %i.\n", rn);
      exit(EXIT_FAILURE);
    }
    /* Output status so far if the verbose level is enough:                 */
    if(rn>0 && verblevel>1)
      tr_output(TOUT_WARN, "checkrange() modified the suggested "
        "parameters, it returned code 0x%x.\n\n", rn);
  }

  /* Get the molecule index for the isotopes:                               */
  setimol(tr);

  /* Check for an opacity file:                                             */
  filecheck = access(th->f_opa, F_OK);
  /* Only read the TLI file if there is no opacity file                     */
  if(filecheck == -1 && rn != -2){
    /* Read data file:                                                      */
    tr_output(TOUT_INFO, "Reading data.\n");
    if((rn=readdatarng(tr, &li)) < 0) {
      tr_output(TOUT_ERROR, "readdatarng returned error code %li.\n", rn);
      exit(EXIT_FAILURE);
    }
    tr_output(TOUT_INFO, "Done.\n\n");
  }
  /* If there is an opacity file, update progress indicator so that
     program may continue:                                                  */
  else{
    tr_output(TOUT_INFO, "Skipping TLI reading.\n");
    tr->pi |= TRPI_READINFO;
    tr->pi |= TRPI_READDATA;
  }

  /* Scale factors:                                                         */
  double fct_to_microns = 1.0/tr->wns.fct/1e-4;
  /* Status so far:                                                         */
  tr_output(TOUT_INFO, "Status so far:\n"
    " * I read %li records from the datafile.\n"
    " * The wavelength range read was %.8g to %.8g microns.\n",
    li.n_l, 1.0/tr->wns.f*fct_to_microns,
    1.0/tr->wns.i*fct_to_microns);
  return 0;
}


/* FUNCTION:
   Frees lineinfo structure
   Return: 0 on success                                                     */
int
freemem_isotopes(struct isotopes *iso,
                 long *pi){
  int i;

  /* Free structures:                                                       */
  if (*pi &= TRPI_READBIN){
    for(i=0; i < iso->n_i; i++){
      free_isof(iso->isof+i);
      free_isov(iso->isov+i);
    }
    for(i=0; i < iso->n_db; i++)
      free_db(iso->db+i);

    /* Free arrays:                                                         */
    free(iso->isov);
    free(iso->isof);
    free(iso->db);
    free(iso->imol);
    free(iso->isoratio);
    /* Unset appropiate flags:                                              */
    *pi &= ~(TRPI_READBIN);
  }

  /* Unset flags:                                                           */
  *pi &= ~(TRPI_READINFO | TRPI_READDATA | TRPI_GETATM);
  return 0;
}


/* FUNCTION:
   Free lineinfo structure.
   Return: 0 on success                                                     */
int
freemem_lineinfo(struct lineinfo *li,
                 long *pi){
  int i;

  //transitprint(1,2, "%ld\n", *pi &= TRPI_READINFO);
  //transitprint(1,2, "%ld\n\n", *pi &= TRPI_READBIN);
  if (*pi &= TRPI_READBIN){
    /* Free isov, dbnoext and samp in li:                                     */
    free_isov(li->isov);
    free(li->isov);

    for(i=0; i<li->ndb; i++)
      free_dbnoext(li->db+i);
    free(li->db);

    /* Zero all the structure:                                                */
    memset(li, 0, sizeof(struct lineinfo));
    /* Unset appropiate flags:                                                */
    *pi &= ~(TRPI_READBIN);
  }

  /* Unset appropiate flags:                                                */
  *pi &= ~(TRPI_READINFO);
  return 0;
}


/* FUNCTION                                                                 */
int
freemem_linetransition(struct line_transition *lt,
                       long *pi){
  /* Free the four arrays of lt:                                            */
  free(lt->wl);
  free(lt->elow);
  free(lt->gf);
  free(lt->isoid);

  /* Unset appropiate flags:                                                */
  *pi &= ~TRPI_READDATA;
  return 0;
}


/* \fcnfh
   Saves line information */
void
saveline(FILE *fp,
         struct lineinfo *li){
}


#ifdef DBGREADLINEINFO
/* \fcnfh
   main function for debugging only */
int main(int argc, char **argv){
  struct transit tr;
  struct transithint th;
  struct lineinfo *li;
  int i,ti1,nbins,ans;
  PREC_LNDATA *ltgf;
  PREC_LNDATA *ltelow;
  PREC_LNDATA *ltwl,*twl;
  short *ltisoid,*tisoid;

  tr.ds.th = &th;
  th.na = 0;

  verblevel = 20;

  th.m = 0.001;
  th.na |= TRH_WM;
  char defile_line[] = "./res/lineread.tli";
  th.f_line = (char *)calloc(strlen(defile_line)+1, sizeof(char));
  strcpy(th.f_line, defile_line);

  th.na |= TRH_FL;

  nbins=20;
  Pprintf(2, "Number of bins[%i]?: ", nbins);
  if(Pgeti(0, &ti1, 6)>0)
    nbins = ti1;

  if((i=readlineinfo(&tr))!=0) {
    tr_output(TOUT_ERROR, "Error code: %i.\n", i);
    exit(EXIT_FAILURE);
  }
  tr_output(TOUT_DEBUG, "range: %.10g to %.10g.\n", tr.ds.li->wi, tr.ds.li->wf);
  li = tr.ds.li;
  ltgf = tr.ds->lt.gf;
  ltwl = tr.ds->lt.wl;
  ltisoid = tr.ds->lt.isoid;
  ltelow = tr.ds->lt.elow;

  ti1 = (int)(log10(li->n_l)+1);

  printf("Done reading the file.\n\n"
         "dbread_pands() test results:\n");
  printf("Chosen wavelength range was from %.10g to %.2f [nm].\n"
         " %*li lines read.\n"
         " Choosing %i equal-sized bins, the result is\n",
         li->wi, li->wf, ti1, li->n_l, nbins);

  long qb[tr.n_i];
  float szb = (li->wf-li->wi)/nbins;
  double endb;

  twl = ltwl;
  tisoid = ltisoid;
  if(!nbins)
    /* FINDME: is this a typo? */
    Pprintf(1, "  hmmm, you chose 0 bins!.\n");
  for(i=0; i<nbins; i++){
    memset(qb, 0, sizeof(*qb)*4);
    endb = li->wi+(i+1)*szb;
    //    PPprintf(1,2,"KK %g %f\n",lp->wl,endb);
    while(*twl<endb && twl-ltwl<li->n_l){
      qb[*tisoid++]++;
      twl++;
    }

    Pprintf(1, " %*i = %i + %i + %i + %i lines shorter than %.3f\n",
            ti1, qb[0]+qb[1]+qb[2]+qb[3], qb[0], qb[1], qb[2], qb[3], endb);
  }

  Pprintf(1, "\nWanna know the value of a single record?\n"
             "If so, write record number (range 0 - %i), else "
             "press ^C: ", li->n_l-1);

  while(Pgeti(0,&ans,(int)(log10(li->n_l))+1)>=0){
    if(ans<li->n_l&&ans>=0){
      Pprintf(1, "Wavelength: %.10g\n", ltwl[ans]);
      Pprintf(1, "Lower Energy Level: %.10g\nLog(gf): %.10g\n",
              ltelow[ans], ltgf[ans]);
      printf("Isotope Name: %s\n", tr.isof[ltisoid[ans]].n);
    }
    else
      Pprintf(1, "\nInvalid record number, so ...");

    Pprintf(1, "\nWanna know the value of another single record?\n"
               "If so, write the record number (range 0 - %i), else just "
               "press ^C: ", li->n_l-1);
  }

}

#undef checkprepost

#endif
