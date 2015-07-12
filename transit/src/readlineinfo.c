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

  transitDEBUG(30, verblevel, "BS: Start looking from %li in %li fields "
                              "for %f\n", offs, nfields, target);
  /* Binary search:                                                         */
  do{
    loc = (hi+lo)/2;                           /* Mid record's index        */
    fseek(fp, offs+reclength*loc, SEEK_SET);   /* Move pointer              */
    fread(&temp, sizeof(PREC_LNDATA), 1, fp);  /* Read value                */
    transitDEBUG(30, verblevel, "BS: found wl %.8f microns at position %li\n",
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
  transitprint(1, verblevel, "Binary search found wavelength: %.8f at "
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
  if(li->tli_ver != compattliversion)
    transiterror(TERR_SERIOUS,
                 "The version of the TLI file: %i (lineread v%i.%i) is not "
                 "compatible with this version of transit, which can only "
                 "read version %i.\n", li->tli_ver, li->lr_ver,
                 li->lr_rev, compattliversion);

  /* Read initial wavelength, final wavelength, and number of databases:    */
  fread(&iniw, sizeof(double), 1, fp);
  fread(&finw, sizeof(double), 1, fp);
  transitprint(1, verblevel, "Initial wavelength: %.2f (um)\n"
                             "Final   wavelength: %.2f (um)\n", iniw, finw);
  fread(&ndb, sizeof(unsigned short), 1, fp);
  transitprint(1, verblevel, "Number of databases: %d.\n", ndb);

  /* Allocate pointers according to the number of databases:                */
  iso->db       = (prop_db      *)calloc(ndb, sizeof(prop_db     ));
  li->db        = (prop_dbnoext *)calloc(ndb, sizeof(prop_dbnoext));
  /* Allocate isotope info.  Start with size 1, then reallocate as needed   */
  iso->isof     = (prop_isof    *)calloc(1,   sizeof(prop_isof));
  li->isov      = (prop_isov    *)calloc(1,   sizeof(prop_isov));
  iso->isoratio = (double       *)calloc(1,   sizeof(double));
  /* Min and max allowed temperatures in TLI files:                         */
  li->tmin      = -1.0;
  li->tmax      = 10000.0;

  /* Read info for each database:                                           */
  for(i=0; i<ndb; i++){
    /* Allocate and get DB's name:                                          */
    fread(&rs, sizeof(unsigned short), 1, fp);         /* Get DBname length */
    iso->db[i].n = (char *)calloc(rs+1, sizeof(char)); /* Allocate          */
    fread(iso->db[i].n, sizeof(char), rs, fp);         /* Read              */
    iso->db[i].n[rs] = '\0';
    transitprint(30, verblevel, "  DB name size: %d'\n", rs);
    /* Read and allocate the molecule's name:                               */
    fread(&rs, sizeof(unsigned short), 1, fp);
    iso->db[i].molname = (char *)calloc(rs+1, sizeof(char));
    fread(iso->db[i].molname, sizeof(char), rs, fp);
    iso->db[i].molname[rs] = '\0';

    transitprint(2,  verblevel, "Database (%d/%d) name: '%s' (%s molecule)\n",
                                i+1, ndb, iso->db[i].n, iso->db[i].molname);

    /* Get number of temperatures and isotopes:                             */
    fread(&nT,     sizeof(unsigned short), 1, fp);
    fread(&nDBiso, sizeof(unsigned short), 1, fp);
    transitprint(2, verblevel, "  Number of temperatures: %d\n"
                               "  Number of isotopes:     %d\n", nT, nDBiso);
    li->db[i].t  = nT;
    iso->db[i].i = nDBiso;

    /* Allocate for the different temperature points and read:              */
    T = li->db[i].T = (double *)calloc(nT, sizeof(double));
    fread(T, sizeof(double), nT, fp);
    transitprint(3, verblevel, "  Temperatures: [%6.1f, %6.1f, ..., %6.1f]\n",
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

    transitDEBUG(21, verblevel, "So far, cumIsotopes: %i, at databases: %i, "
                 "position %li.\n", correliso+nDBiso, i, ftell(fp));

    /* Display database general info:                                       */
    transitDEBUG(23, verblevel, "DB %i: '%s' has %i (%i) temperatures, "
                 "%i (%i) isotopes, and starts at cumulative isotope %i.\n",
                 iso->isof[correliso].d, iso->db[i].n,
                 li->db[i].t, nT,
                 iso->db[i].i, nDBiso, iso->db[i].s);


    /* Read the isotopes from this data base:                               */
    for (j=0; j<nDBiso; j++){
      /* Store isotopes'  DB index number:                                  */
      iso->isof[correliso].d = i;
      transitprint(10, verblevel,"  Correlative isotope number: %d", correliso);

      /* Read isotopes' name:                                               */
      fread(&rs, sizeof(unsigned short), 1, fp);
      iso->isof[correliso].n = (char *)calloc(rs+1, sizeof(char));
      fread(iso->isof[correliso].n, sizeof(char), rs, fp);
      iso->isof[correliso].n[rs] = '\0';
      transitprint(2,  verblevel, "  Isotope (%i/%i): '%s'\n", j+1, nDBiso,
                                  iso->isof[correliso].n);
      transitprint(30, verblevel, "   Isotope name size: %d'\n", rs);

      /* Read mass and isotopic ratio:                                      */
      fread(&iso->isof[correliso].m,   sizeof(double), 1, fp);
      fread((iso->isoratio+correliso), sizeof(double), 1, fp);
      transitprint(3,  verblevel, "    Mass:  %g u (%g gr)\n",
                          iso->isof[correliso].m, iso->isof[correliso].m*AMU);
      transitprint(3,  verblevel, "    Isotopic ratio: %.4g\n",
                          iso->isoratio[correliso]);
      transitDEBUG(30, verblevel, "    File position: %li.\n", ftell(fp));

      /* Read partition function:                                           */
      Z  = li->isov[correliso].z = li->isov[correliso-j].z + nT*j;
      fread(Z,  sizeof(double), nT, fp);
      li->isov[correliso].n = nT;

      transitprint(10, verblevel, "    Part Function:    [%.2e, %.2e, ..., "
                                  "%.2e]\n", Z[0],  Z[1],  Z[nT-1]);
      correliso++;
    }

    /* Update cumulative isotope count (index of first isitope in this DB): */
    iso->db[i].s = niso;
    niso += nDBiso;
  }

  transitprint(3, verblevel, "Cumulative number of isotopes per DB: [");
  for (i=0; i<ndb; i++)
    transitprint(3, verblevel,"%2d, ", iso->db[i].s);
  transitprint(3, verblevel, "\b\b].\n");
  transitprint(3, verblevel, "acum Iso: %2d.\n", niso);

  //transitprint(3, verblevel, "Iso ratio: %.5g %.5g %.5g %.5g\n", iso->isoratio[0], iso->isoratio[1], iso->isoratio[2], iso->isoratio[3]);

  /* Update structure values:                                               */
  li->ni = iso->n_i = niso; /* Number of isotopes                           */
  li->ndb = ndb;            /* Number of databases                          */
  /* Position of first line data (there's still one integer to be read):    */
  li->endinfo = ftell(fp) + sizeof(int);
  li->wi = iniw;            /* Initial wavelength                           */
  li->wf = finw;            /* Final wavelength                             */
  iso->n_db = ndb;          /* Number of databases                          */
  /* Allocate isotope's variable data                                       */
  iso->isov = (prop_isov *)calloc(iso->n_i, sizeof(prop_isov));

  return 0;
}


/* FUCTION:
   Set the molecule's index of the isotopes:                                */
int
setimol(struct transit *tr){
  struct molecules *mol = tr->ds.mol;
  struct isotopes  *iso = tr->ds.iso;
  int i; /* Auxiliary for-loop index                                        */

  iso->imol = (int *)calloc(iso->n_i, sizeof(int));
  iso->nmol = 0;
  for (i=0; i<iso->n_i; i++){
    //transitprint(1, verblevel, "Isotope %d is '%s', from DB %d: '%s'.\n",
    //                  i, iso->isof->n, iso->isof->d, iso->db[iso->isof->d].n);
    /* Search water molecule's index:                                       */
    iso->imol[i] = findstring(iso->db[iso->isof[i].d].molname, mol->name,
                              mol->nmol);
    transitprint(30, verblevel, "Isotope '%s', is mol %d: '%s'.\n",
                 iso->isof[i].n, iso->imol[i], mol->name[iso->imol[i]]);

    /* Find if current imol is already in imol array:                       */
    if (valueinarray(iso->imol, iso->imol[i], i) < 0){
      transitprint(20, verblevel, "Isotope %s (%d) is a new molecule (%s).\n",
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
  struct transithint *th = tr->ds.th;       /* transithint                  */
  prop_samp *msamp = &li->wavs;             /* transit wavelength sampling  */
  prop_samp *hsamp = &th->wavs;             /* hint    wavelength sampling  */
  PREC_LNDATA dbini = li->wi*TLI_WAV_UNITS, /* Minimum DB wavelength        */
              dbfin = li->wf*TLI_WAV_UNITS; /* Maximum DB wavelength        */
  double cm_to_micron = 1e4,                /* Conversion factor to microns */
         fct;

  /* FINDME: hack prints: */
  //transitprint(1, verblevel, " hsamp->f: %g,  hsamp->i: %g\n",
  //               hsamp->f, hsamp->i);
  //transitprint(1, verblevel, " hsamp->fct: %g\n", hsamp->fct);
  //transitprint(1, verblevel, " db i: %g,  db f: %g\n", dbini, dbfin);

  /* Initialize modified hints:                                             */
  msamp->n = -1;
  msamp->d = -1;
  msamp->v = NULL;
  msamp->fct = 1;

  /* Check that the hinted wavelength units factor is positive & non-zero:  */
  if(hsamp->fct <= 0)
    transiterror(TERR_SERIOUS, "User specified wavelength factor is "
                               "negative (%g).\n", hsamp->fct);

  /* Set lineinfo wavelength units factor equal to transithint factor:      */
  msamp->fct = hsamp->fct;

  /* transit lineinfo.wavs conversion factor to cgs:                        */
  fct = msamp->fct;

  transitDEBUG(10, verblevel,
               "Hinted initial and final wavelengths are %6g and %6g cm.\n"
               "The database max and min wavelengths are %6g and %6g cm.\n",
               hsamp->i*fct, hsamp->f*fct, dbini, dbfin);

  /* Set final wavelength:                                         */
  /* If invalid/not set hint final wavelength, default it to zero: */
  if(hsamp->f < 0){
    hsamp->f = 0;
    transiterror(TERR_WARNING,
                 "Incorrect upper wavelength limit in hint.  Default: setting "
                 "to %g before extraction.\n", hsamp->f*fct);
  }
  /* If hint is 0, set it to max db wavelength:                             */
  if(hsamp->f <= 0){
      msamp->f = dbfin/fct;
  }
  else{  /* Else, hinted f is a positive value:                             */
    transitDEBUG(20, verblevel, "dbini: %g  sampf: %g.\n",
                 dbini, hsamp->f);
    /* Check that it is not below the minimum value:                        */
    if(dbini > hsamp->f * fct){
      transiterror(TERR_SERIOUS|TERR_ALLOWCONT, "Final wavelength (%g * %g) "
                   "is smaller than minimum wavelength in database (%g).\n",
                   hsamp->f, fct, dbini);
      return -3;
    }
    /* Warn if it is above maximum value with information:                  */
    if(hsamp->f * fct > dbfin)
      transiterror(TERR_WARNING, "Final requested wavelength (%g microns) "
                   "is larger than the maximum informative value in database "
                   "(%g microns).\n", hsamp->f*fct * cm_to_micron,
                                      dbfin        * cm_to_micron);
    /* Set the final wavelength value:                                      */
    msamp->f = hsamp->f;
  }
  /* Set initial wavelength:                                                */
  /* If invalid value, default it to 0:                                     */
  if(hsamp->i < 0){
    hsamp->i = 0;
    transiterror(TERR_WARNING, "Setting hinted lower wavelength limit "
                 "before extraction as %g cgs. It was not user-hinted.\n",
                 hsamp->i*fct);
  }
  /* If default value, set it to min db wavelength:                         */
  if(hsamp->i<=0)
    msamp->i = dbini/fct;
  else{
    transitDEBUG(20, verblevel, "dbfin: %g  sampi: %g.\n",
                 dbfin, fct*hsamp->i);
    /* Check that it is not larger than the maximum db wavelength:          */
    if(dbfin < fct * hsamp->i){
      transiterror(TERR_SERIOUS|TERR_ALLOWCONT, "Initial wavelength (%g cm) "
                   "is larger than maximum wavelength in database (%g cm).\n",
                   fct*hsamp->i, dbfin);
      return -2;
    }
    if(fct * hsamp->i < dbini)
      transiterror(TERR_WARNING, "Initial requested wavelength (%g microns) "
                   "is smaller than the minimum informative value in database "
                   "(%g microns).\n", hsamp->i * fct * cm_to_micron,
                                      dbini          * cm_to_micron);
    msamp->i = hsamp->i;
  }

  /* Set progress indicator and return status:                     */
  tr->pi |= TRPI_CHKRNG;
  return res;
}


/* FUNCTION:
   Check TLI file exists.  Check that machine formating is compatible
   with lineread.  Determine if TLI is ASCII or binary.  Read either
   ASCII or binary TLI file. Declare line_transition.

  TD:  Checks on allocation errors.
  Return: 1 on success
         -1 unavailable file
         -2 Filename not hinted
         -3 TLI format not valid (missing magic bytes)
         -4 Improper TLI-ASCII input                              */
int
readinfo_tli(struct transit *tr,
             struct lineinfo *li){
  int rn;
  FILE *fp;  /* File pointer of info file: */

  /* Decalre and initialize the union sign:                         */
  /* sign.s contains the magic numbers of this machine's and TLI's: */
  union {char sig[4]; int32_t s[2];} sign =
    {.s={0, ((0xff-'T')<<24)|((0xff-'L')<<16)|((0xff-'I')<<8)|(0xff)}};
  char line[maxline+1];

  /* Pointer to hint: */
  struct transithint *th = tr->ds.th;

  /* Get TLI file name from hint:                              */
  if(!th->f_line){  /* Check if it was defined in hint         */
    transiterror(TERR_SERIOUS|TERR_ALLOWCONT, "Undefined TLI file name.\n");
    return -2;
  }
  /* Check that the file exists and make a pointer to read it: */
  if((rn=fileexistopen(th->f_line, &fp)) != 1){
    transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
                 "Line info file '%s' is not found. "
                 "fileexistopen() error code %i.\n", th->f_line, rn);
    return -1;
  }
  /* Set transit TLI file pointer and TLI file name: */
  tr->fp_line = fp;
  tr->f_line  = th->f_line;

  /* Read first four bytes, they should be either
  `(0xff-T)(0xff-L)(0xff-I)(0xff)' or '\#TLI'. They are stored as integer.
  This checks whether the machine where the TLI file and the one this
  program is being run have the same endian order.  If the first two are
  '\#TLI', then the first line should also start as '\#TLI-ascii'           */
  fread(sign.s, sizeof(int32_t), 1, fp);

  /* Determine if TLI is binary (asciiline=0) or ASCII (asciiline=1):       */
  li->asciiline = 0;
  transitDEBUG(13, verblevel, "Comparing %i and %i for Magic Number (len: "
                            "%li)\n", sign.s[0], sign.s[1], sizeof(sign.s[0]));

  if(sign.s[0] != sign.s[1]){
    /* Does it look like an ASCII TLI?, if so check it:                     */
    rn = strncasecmp(sign.sig, "#TLI", 4);  /* FINDME: strncasecmp */
    if(!rn){
      strcpy(line, "#TLI");
      fread(line+4, sizeof(char), 6, fp);
      rn = strncasecmp(line, "#TLI-ascii", 10);
    }
    /* If it wasn't a valid TLI, throw error and exit:                      */
    if(rn){
      transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
                   "The file '%s' has not a valid TLI format. It might be "
                   "because the machine were the file was created have "
                   "different endian order, which is incompatible.\n",
                   tr->f_line);
      return -3;
    }
    li->asciiline = 1;
    /* Ignore the rest of the first line: */
    fgetupto_err(line, maxline, fp, &linetoolong, tr->f_line, 1);
  }

  /* Read binary TLI:  */
  if((rn=readtli_bin(fp, tr, li))!=0){
    transiterror(TERR_CRITICAL|TERR_ALLOWCONT,
                 "readtli_bin() return error code %i.\n", rn);
    return -6;
  }
  transitprint(3, verblevel, "TLI file read from %g to %g microns.\n",
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
  int nlines,            /* Number of line transitions                      */
      niso,              /* Number of isotopes in line transition data      */
      nread,             /* Number of transitions to read for each isotope  */
      *isotran,          /* Number of transitions per isotope in TLI        */
      start,             /* Position of first LT for isotope in TLI         */
      i,                 /* for-loop index                                  */
      offset=0,          /* Isotope offset (in number of transitions)       */
      wl_loc, iso_loc,   /* Offsets for isoID, Elow, and gf data            */
      el_loc, gf_loc,    /* (in memory)                                     */
      rn;                /* Return IDs                                      */
  /* Indices of first and last transitions to be stored                     */
  long ifirst, ilast;

  /* Auxiliary variables to keep wavelength limits:                         */
  PREC_LNDATA iniw = li->wavs.i * li->wavs.fct / TLI_WAV_UNITS;
  PREC_LNDATA finw = li->wavs.f * li->wavs.fct / TLI_WAV_UNITS;
  PREC_LNDATA wltmp;   /* Auxiliary variable to store wavelength            */
  (void) wltmp;        /* Suppress warning for this unused variable         */

  /* Open line data file:                                                   */
  if((rn=fileexistopen(tr->f_line, &fp)) != 1){
    transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
                 "Data file '%s' not found.  fileexistopen() error "
                 "code: %i.\n", tr->f_line, rn);
    return -1;
  }

  /* Find starting point in datafile.  First with a binary search, then
     with a sequential search:                                              */

  /* Check seekability:                                                     */
  if(fseek(fp, 0, SEEK_CUR)){
    transiterror(TERR_CRITICAL|TERR_ALLOWCONT,
                 "File '%s' was not seekable.\n", tr->f_line);
    return -2;
  }

  /* Read total number of transitions in TLI file:                          */
  /* FINDME: May be better to put endinfo to the right position
             (avoid this  -sizeof(int)):                                    */
  fseek(fp, li->endinfo - sizeof(int), SEEK_SET);
  fread(&nlines, sizeof(int), 1, fp);
  transitprint(1, verblevel, "TLI has %d transition lines.\n", nlines);

  /* Number of isotopes in line transition data:                            */
  fread(&niso, sizeof(int), 1, fp);
  transitprint(10, verblevel, "TLI has %d isotopes.\n", niso);
  /* Number of transitions per isotope:                                     */
  isotran = calloc(niso, sizeof(int));
  fread(isotran, sizeof(int), niso, fp);
  for (i=0; i<niso; i++){
    transitprint(20, verblevel, "Ntransitions[%d]: %d.\n", i, isotran[i]);
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
  if(!lt->gf || !lt->wl || !lt->elow || !lt->isoid)
    transiterror(TERR_CRITICAL|TERR_ALLOC, "Couldn't allocate memory for "
                 "linetran structure array of length %i, in function "
                 "readdatarng.\n", nlines);

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
    transitprint(5, verblevel, "Initial and final entries are: "
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
    start += isotran[i]*sizeof(double);
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


  /* Read hinted info file:                                                 */
  transitprint(1, verblevel, "Reading info file '%s' ...\n", th->f_line);
  if((rn=readinfo_tli(tr, &li)) != 1)
    transiterror(TERR_SERIOUS, "readinfo_tli() returned an error "
                 "code %i.\n", rn);
  transitprint(1, verblevel, " Done.\n\n");

  /* Get the molecule index for the isotopes:                               */
  /* FINDME: Move this out of readline later.                               */
  rn = setimol(tr);

  /* Check the remainder range of the hinted values
     related to line database reading:                                      */
  if((rn=checkrange(tr, &li)) < 0)
    transiterror(TERR_SERIOUS, "checkrange() returned error code %i.\n", rn);
  /* Output status so far if the verbose level is enough:                   */
  if(rn>0 && verblevel>1)
    transiterror(TERR_WARNING, "checkrange() modified the suggested "
                               "parameters, it returned code 0x%x.\n\n", rn);

  /* Scale factors:                                                         */
  double fct = li.wavs.fct;
  double fct_to_microns = fct/1e-4;
  transitprint(2, verblevel, "The wavelength range to be used is %g to %g "
               "cm.\n", fct*tr->ds.li->wavs.i, fct*tr->ds.li->wavs.f);

  /* Check for an opacity file:                                             */
  filecheck = access(th->f_opa, F_OK);
  /* Only read the TLI file if there is no opacity file                     */
  if(filecheck == -1){
    /* Read data file:                                                      */
    transitprint(1, verblevel, "Reading data.\n");
    if((rn=readdatarng(tr, &li))<1)
      transiterror(TERR_SERIOUS, "readdatarng returned error code %li.\n", rn);
    transitprint(1, verblevel, "Done.\n\n");
  }
  /* If there is an opacity file, update progress indicator so that
     program may continue:                                                  */
  else{
    transitprint(1, verblevel, "Skipping TLI reading.\n");
    tr->pi |= TRPI_READINFO;
    tr->pi |= TRPI_READDATA;
  }

  /* Status so far:                                                         */
  transitprint(2, verblevel, "Status so far:\n"
               " * I read %li records from the datafile.\n"
               " * The wavelength range read was %.8g to %.8g microns.\n",
               li.n_l, li.wavs.i*fct_to_microns, li.wavs.f*fct_to_microns);

  transitDEBUG(21, verblevel,
               "Database min and max: %.10g(%.10g) and %.10g(%.10g)\n",
               li.wi, tr->ds.li->wi, li.wf, tr->ds.li->wf);
  return 0;
}


/* FUNCTION:
   Frees lineinfo structure
   Return: 0 on success                                      */
int
freemem_isotopes(struct isotopes *iso,
                 long *pi){
  int i;

  /* Free structures:                                         */
  for(i=0; i < iso->n_i; i++){      /* Allocated in readlineinfo */
    free_isof(iso->isof+i);
    free_isov(iso->isov+i);
  }
  for(i=0; i < iso->n_db; i++)
    free_db(iso->db+i);

  /* Free arrays:                                             */
  free(iso->isov);
  free(iso->isof);
  free(iso->db);
  free(iso->imol);
  free(iso->isoratio);

  /* Unset flags:                                             */
  *pi &= ~(TRPI_READINFO | TRPI_READDATA | TRPI_CHKRNG | TRPI_GETATM);
  return 0;
}


/* FUNCTION:
   Free lineinfo structure.
   Return: 0 on success                                                     */
int
freemem_lineinfo(struct lineinfo *li,
                 long *pi){
  int i;

  /* Free isov, dbnoext and samp in li:                                     */
  free_isov(li->isov);
  free(li->isov);

  for(i=0; i<li->ndb; i++)
    free_dbnoext(li->db+i);
  free(li->db);

  free_samp(&li->wavs);

  /* Zero all the structure:                                                */
  memset(li, 0, sizeof(struct lineinfo));

  /* Unset appropiate flags:                                                */
  *pi &= ~(TRPI_READINFO | TRPI_CHKRNG);
  return 0;
}

/* FUNCTION  */
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

  if((i=readlineinfo(&tr))!=0)
    transiterror(TERR_CRITICAL, "Error code: %i.\n", i);
  transitDEBUG(20, verblevel, "range: %.10g to %.10g.\n",
               tr.ds.li->wi, tr.ds.li->wf);
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
