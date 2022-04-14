// Copyright (C) 2015-2016 University of Central Florida. All rights reserved.
// Transit is under an open-source, reproducible-research license (see LICENSE).

/* Revision        March 19th,   2014 Jasmina Blecic
                   added eclipse_ray_solution                              */
/* Revision        April 26th,   2014 Jasmina Blecic
                   added new function for ray grid                         */

/* Required for compilation with shared memory libraries.                   */
#define _XOPEN_SOURCE 700

#include <time.h>
#include <transit.h>
#include <version_tr.h>
#include <math.h>
#include <procopt.h>
#include <libgen.h>
#include <string.h>
#include <assert.h>
#include <sys/stat.h>


/* Transit ray solutions:                                                   */
const static ray_solution *raysols[] = {&slantpath,   /* see slantpath.c    */
                                        &eclipsepath, /* see eclipse.c      */
                                        NULL};

#ifndef EXTRACFGFILES
#  define PREPEXTRACFGFILES ""
#endif

#ifdef NODOTDEFAULT
#  define DOTCFGFILE ""
#  define DOTCFGFILENM "NO DEFAULT FILE"
#else
#  define DOTCFGFILE "./.transitrc"
#  define DOTCFGFILENM DOTCFGFILE
#  ifdef EXTRACFGFILES
#    define PREPEXTRACFGFILES ","EXTRACFGFILES
#  endif
#endif

/* \fcnfh
   Process command line options, saving them in the hint structure.

   Return: 0 on success                                                */
int
processparameters(int argc,            /* Number of command-line args  */
                  char **argv,         /* Command-line arguments       */
                  struct transit *tr){ /* struct to store hinted pars  */
  enum param {       /* List of command-line arguments IDs: */
    CLA_DUMMY=128,
    CLA_ATMOSPHERE,
    CLA_LINEDB,
    CLA_RADLOW,
    CLA_RADHIGH,
    CLA_RADDELT,
    CLA_RADFCT,
    CLA_WAVLOW,
    CLA_WAVHIGH,
    CLA_WAVFCT,
    CLA_WAVNLOW,
    CLA_WAVNHIGH,
    CLA_WAVNDELT,
    CLA_WAVNOSAMP,
    CLA_WNFCT,
    CLA_ALLOWQ,
    CLA_GORBPAR,
    CLA_GORBPARFCT,
    CLA_TOOMUCH,
    CLA_OUTTOOMUCH,
    CLA_OUTSAMPLE,
    CLA_OUTSPEC,
    CLA_OUTINTENS,
    CLA_TAULEVEL,
    CLA_MODLEVEL,
    CLA_ETHRESH,
    CLA_CLOUD,
    CLA_CLOUDTOP,
    CLA_SCATTERING,
    CLA_TRANSPARENT,
    CLA_DETEXT,
    CLA_DETCIA,
    CLA_DETTAU,
    CLA_CSFILE,
    CLA_SAVEEXT,
    CLA_STARRAD,
    CLA_SOLUTION_TYPE,
    CLA_INTENS_GRID,
    CLA_OPACITYFILE,
    CLA_TEMPLOW,
    CLA_TEMPHIGH,
    CLA_TEMPDELT,
    CLA_MOLFILE,
    CLA_RPRESS,
    CLA_RRADIUS,
    CLA_GSURF,
    CLA_OPABREAK,
    CLA_OPASHARE,
    CLA_NDOP,
    CLA_NLOR,
    CLA_DMIN,
    CLA_DMAX,
    CLA_LMIN,
    CLA_LMAX,
    CLA_QSCALE,
    CLA_QMOL,
    CLA_SAVEFILES,
  };

  /* Generate the command-line option parser: */
  struct optdocs var_docs[]={
    /* General options:                       */
    {NULL,       0,   HELPTITLE,    NULL, NULL,
     "General Arguments:"},
    {"version",  'V', no_argument,  NULL, NULL,
     "Display Transit's version number."},
    {"help",     'h', no_argument,  NULL, NULL,
     "Display the list of command-line arguments."},
    /* FINDME: defaults not implemented                                     */
    //{"defaults", 'd', no_argument,  NULL, NULL,
    // "Prints default values of the different variable."},
    {"quiet",    'q', no_argument,  NULL, NULL,
     "Set the verbosity level to the minimum."},
    {"verb",     'v', required_argument,  "2", "verb",
     "Set the verbosity level (integer) to <verb>."},
    {"config_file",   'c', ADDPARAMFILE, NULL, "file",
     "Read command-line arguments from <file>."
     " '" DOTCFGFILENM PREPEXTRACFGFILES"'."},

    /* Input and output options:              */
    {NULL,          0,             HELPTITLE,         NULL, NULL,
     "INPUT/OUTPUT OPTIONS:"},
    {"atm",        CLA_ATMOSPHERE, required_argument, "NULL",  "atmfile",
     "File containing atmospheric info (Radius, pressure, temperature). A dash"
     " (-) indicates alternative input."},
    {"linedb",     CLA_LINEDB,     required_argument, NULL,
     "linedb", "File containing line information (TLI format, as given by "
               "'pylineread'."},
    {"outtoomuch", CLA_OUTTOOMUCH, required_argument, NULL, "filename",
     "Ouputs depth where toomuch optical depth has been attained as a function"
     " of wavelength."},
    {"outsample",  CLA_OUTSAMPLE,  required_argument, NULL, "filename",
     "Outputs sampling information. A dash (-) indicates standard input. By "
     "default there is no such output."},
    {"outspec",  CLA_OUTSPEC,  required_argument, "outspectrum", "filename",
     "Output modulation (transit) or flux (eclipse) spectrum file."},
    {"outintens",  CLA_OUTINTENS,  required_argument, NULL, "filename",
     "Outputs intensity information. A dash (-) indicates standard input. By "
     "default there is no such output."},
    {"molfile",    CLA_MOLFILE,    required_argument, "../inputs/molecules.dat",
     "filename", "Path to file with the molecular info."},
    {"savefiles", CLA_SAVEFILES, required_argument, NULL, "no",
     "Save output files with tau, molecular extinction, CIA, and total "
     "extinction, with the predetermined file names."},

    /* Radius options:                                                      */
    {NULL,       0,           HELPTITLE,         NULL, NULL,
     "RADIUS OPTIONS (defaults adopt the values from the atmosphere file):"},
    {"raddelt", CLA_RADDELT, required_argument, "-1", "spacing",
     "Radius spacing.  If set, resample the atmospheric layers to an "
     "equi-distant radius sampling."},
    {"radlow",  CLA_RADLOW,  required_argument, "0", "radius",
     "Lower radius.  If 0, use atmospheric data minimum."},
    {"radhigh", CLA_RADHIGH, required_argument, "0", "radius",
     "Higher radius.  If 0, use atmospheric data maximum."},
    {"radfct",  CLA_RADFCT,  required_argument, "0",  "factor",
     "Radius factor.  Multiplicating radius values by this gives centimeters. "
     "If 0, use atmosphere-file factor."},

    /* Atmosphere options:                                                  */
    {NULL,                0,            HELPTITLE,         NULL, NULL,
     "ATMOSPHERE OPTIONS:"},
    {"allowq",            CLA_ALLOWQ,   required_argument, "0.00001", "value",
     "Maximum allowed cumulative-abundance departure from 1.0."},
    {"refpress",          CLA_RPRESS,   required_argument, NULL, NULL,
     "Reference pressure of the planet's 'surface'."},
    {"refradius",         CLA_RRADIUS,  required_argument, NULL, NULL,
     "Reference altitude of the planet's 'surface'."},
    {"gsurf",             CLA_GSURF,    required_argument, NULL, NULL,
     "Surface gravity in cm/s^2."},
    {"qmol",              CLA_QMOL,     required_argument, NULL, NULL,
     "List of molecule names to modify their abundace with qscale."},
    {"qscale",            CLA_QSCALE,   required_argument, NULL, NULL,
     "log10-abundance scale factors for molecules in qmol."},

    /* Wavelength options:                                                  */
    {NULL,         0,             HELPTITLE,         NULL,       NULL,
     "WAVELENGTH OPTIONS (in fct units):"},
    {"wllow",     CLA_WAVLOW,    required_argument, NULL,        "wavel",
     "Lower wavelength. 0 if you want to use line data minimum."},
    {"wlhigh",    CLA_WAVHIGH,   required_argument, NULL,        "wavel",
     "Upper wavelength. 0 if you want to use line data maximum."},
    {"wlfct",     CLA_WAVFCT,    required_argument, "1e-4",      "factor",
     "Wavelength factor. Multiplicating wavelength values by this gives "
     "centimeters. If 0 or 1 then use centimeters."},

    /* Wavenumber options:                    */
    {NULL,         0,              HELPTITLE,         NULL, NULL,
     "WAVENUMBER OPTIONS (in fct units):"},
    {"wnlow",     CLA_WAVNLOW,    required_argument, NULL,  "waven",
     "Lower wavenumber. 0 if you want to use equivalent of the wavelength "
     "maximum."},
    {"wnhigh",    CLA_WAVNHIGH,   required_argument, NULL,  "waven",
     "Upper wavenumber. 0 if you want to use equivalent of the wavelength "
     "minimum."},
    {"wndelt",    CLA_WAVNDELT,   required_argument, "0",  "spacing",
     "Wavenumber spacing. 0 if you want to have the same number of points "
     "as in the wavelength sampling."},
    {"wnosamp",   CLA_WAVNOSAMP,  required_argument, "2160",  "integer",
     "Wavenumber oversampling. 0 if you want the same value as for the "
     "wavelengths."},
    {"wnfct",     CLA_WNFCT,      required_argument, "0",  "factor",
     "Output wavenumber factor. Multiplicating wavenumber values by this "
     "gives centimeters. If 0 then use wavelength's value. This only applies "
     "to output, internally wavenumbers will always be in cm-1."},

    /* Voigt profile options:                                               */
    {NULL, 0, HELPTITLE, NULL, NULL, "Voigt profile options:"},
    {"ndop",    CLA_NDOP, required_argument, "60",   "integer",
    "Number of Doppler-broadening width samples."},
    {"nlor",    CLA_NLOR, required_argument, "60",   "integer",
    "Number of Lorentz-broadening width samples."},
    {"dmin",    CLA_DMIN, required_argument, "1e-3", "float",
    "Minimum Doppler-broadening width (in cm-1)."},
    {"dmax",    CLA_DMAX, required_argument, "0.25", "float",
    "Maximum Doppler-broadening width (in cm-1)."},
    {"lmin",    CLA_LMIN, required_argument, "1e-4", "float",
    "Minimum Lorentz-broadening width (in cm-1)."},
    {"lmax",    CLA_LMAX, required_argument, "10.0", "float",
    "Maximum Lorentz-broadening width (in cm-1)."},
    {"nwidth",  'a',      required_argument, "20",   "number",
     "Number of the max-widths (the greater of Voigt or Doppler widths) "
     "that needs to be contained in a calculated profile."},

    /* Extinction calculation options:                                      */
    {NULL,         0,               HELPTITLE,         NULL,    NULL,
     "EXTINCTION CALCULATION OPTIONS:"},
    {"ethreshold", CLA_ETHRESH,   required_argument, "1e-8",    "ethreshold",
     "Minimum extinction-coefficient ratio (w.r.t. maximum in a layer) to "
     "consider in the calculation."},
    {"cloud",      CLA_CLOUD,      required_argument, NULL,
     "cloudtype,cloudext,cloudtop,cloudbot",
     "Gray-opacity layer with extinction linearly increasing from 0 at "
     "cloudtop, up cloudext at cloudbot.  Then keep a constant extinction "
     "until the bottom of the atmosphere.  cloudext has units of cm-1, "
     "cloudtop and cloudbot units are given by radfct."},
    {"cloudtop", CLA_CLOUDTOP, required_argument, NULL,
     "cloudtop", "Cloud-deck top pressure (log-bar)."},
    {"scattering", CLA_SCATTERING, required_argument, NULL,
     "scattering", "Scattering parameter.  "
     "Options: 'polar' for Villanueva et al. (2022) method; "
     "a float, giving the extinction in some log-units."},
    {"detailext",  CLA_DETEXT,      required_argument, NULL,
     "filename:wn1,wn2,...",
     "Save extinction at specified wavenumbers in filename."},
    {"detailcia",  CLA_DETCIA,      required_argument, NULL,
     "filename:wn1,wn2,...",
     "Save extinction due to CIA at specified wavenumbers in filename."},
    {"csfile",      CLA_CSFILE,      required_argument, NULL,    "filenames",
     "Use the indicated filenames for cross-section opacities (comma-"
     "separated list)."},
    {"saveext",    CLA_SAVEEXT,     required_argument, NULL,    "filename",
     "Save extinction array in this file which won't need to be recomputed "
     "if only the radius scale (scale height) changes."},

    /* Opacity grid options:                                                */
    {NULL,        0,            HELPTITLE,         NULL, NULL,
     "OPACITY GRID OPTIONS:"},
    {"opacityfile",    CLA_OPACITYFILE, required_argument,   NULL,  "filename",
     "Filename to read/save the opacity grid."},
    {"tlow",   CLA_TEMPLOW,    required_argument,  "500",  "temperature",
     "Lower temperature sample (in kelvin)."},
    {"thigh",  CLA_TEMPHIGH,   required_argument, "3000",  "temperature",
     "Upper temp"},
    {"tempdelt",  CLA_TEMPDELT,   required_argument,  "100.0",  "spacing",
     "Temperature sample spacing (in kelvin)."},
    {"justOpacity",      CLA_OPABREAK,  no_argument, NULL, NULL,
     "If set, End execution after the opacity-grid calculation."},
    {"shareOpacity",      CLA_OPASHARE,  no_argument, NULL, NULL,
     "If set, attempt to place the opacity grid into shared memory."},

    /* Resulting ray options:                 */
    {NULL,        0,            HELPTITLE,         NULL, NULL,
     "RESULTING RAY OPTIONS:"},
    {"solution",  's',          required_argument, "eclipse", "sol_name",
     "Name of the kind of output solution (eclipse or transit)."},
    {"toomuch",   CLA_TOOMUCH,  required_argument, "20", "optdepth",
     "If optical depth for a particular path is larger than optdepth, then do "
     "not proceed to lower radius."},
    {"taulevel",  CLA_TAULEVEL, required_argument, "1",  "integer",
     "Use constant (1) or variable (2) index of refraction for the transit "
     "ray path."},
    {"modlevel",  CLA_MODLEVEL, required_argument, "1",  "integer",
     "Do an integration of level <integer> to compute modulation. 1 doesn't "
     "consider limb darkening. -1 doesn't consider limb darkening and "
     "additionally only returns the moduated radius at which extinction "
     "becomes one."},
    {"detailtau", CLA_DETTAU,   required_argument, NULL, "filename:wn1,wn2,..",
     "Save optical depth at specified wavenumbers in filename"},

    /* Geometry options:                      */
    {NULL,          0,               HELPTITLE,         NULL,    NULL,
     "GEOMETRY PARAMETERS"},
    {"starrad",     CLA_STARRAD,     required_argument, "1.125", "radius_sun",
     "Stellar radius in solar radius."},
    {"gorbpar",    CLA_GORBPAR,     required_argument, NULL,
     "smaxis,time,incl,ecc,long_node,arg_per",
     "Orbital parameters. Use the above order. Default: 1, 0, 0, 0, 0, 0."},
    {"gorbparfct", CLA_GORBPARFCT,  required_argument, NULL,
     "unitsof:smaxis,time,incl,ecc,long_node,arg_per",
     "Units convertion factors to the cgs system of the orbital parameters. "
     "Same order of g-orbpar.  Default: AU, hours, deg, 1, deg, deg."},
    {"transparent", CLA_TRANSPARENT, no_argument,       NULL,    NULL,
     "If selected, the planet will have a maximum optical depth given by "
     "toomuch, it will never be totally opaque."},
    {"raygrid",      CLA_INTENS_GRID, required_argument, "0 20 40 60 80",
     NULL, "Intensity grid"},
    {NULL, 0, 0, NULL, NULL, NULL}
  };

  struct optcfg var_cfg;
  memset(&var_cfg, 0, sizeof(var_cfg));
  var_cfg.contact ="\n  Joseph Harrington <jh@physics.ucf.edu>\n"
                     "  Patricio Cubillos <pcubillos@fulbrightmail.org>\n"
                     "  Jasmina Blecic    <jasmina@physics.ucf.edu>";

  var_cfg.files   = DOTCFGFILE PREPEXTRACFGFILES;
  var_cfg.columns = 70;

  /* Declare transithint:                                                   */
  static struct transithint st_trh;
  memset(&st_trh, 0, sizeof(struct transithint));
  /* Assign pointer-to-transithint to transit data structure (ds):          */
  struct transithint *hints = tr->ds.th = &st_trh;
  /* Setup flags, verbose level, and wheter abundance is by mass or number: */
  st_trh.fl |= TRU_ATMASK1P | TRU_SAMPSPL | TRH_MASS;
  st_trh.verbnoise = 4;
  st_trh.mass = 1;
  st_trh.savefiles = 0;

  int rn,  /* optdocs' short option */
      i;   /* Auxilliary index      */
  char name[20];
  /* For interactive rad, wavelength, and wavenumber sample inputs: */

  int isYes, /* if savefiles is 'yes' */
       isNo; /* if savefiles is 'no' */

  /* Preset names for detailed output in transithint: */
  struct detailfld *det = &hints->det.tau;
  strcpy(det->name, "Optical depth");
  det = &hints->det.ext;
  strcpy(det->name, "Extinction");
  det = &hints->det.cia;
  strcpy(det->name, "CIA extinction");

  /* Start processing the command options:  */
  procopt_debug = 1;
  opterr = 0;
  while(1){
    /* Get an option:                       */
    rn = procopt(argc, argv, var_docs, &var_cfg, NULL);
    /* Exit statement:                      */
    if (rn==-1)
      break;

    /* Print info to screen if debugging:   */
    tr_output(TOUT_DEBUG, "Processing option '%c' / %d, argum: %s\n",
      rn, rn, optarg);
    /* Handle each case:                                                    */
    switch(rn){
    /* Cross-section data files:                                            */
    case CLA_CSFILE:
      hints->ncross  = nchar(optarg, ',') + 1;        /* Count files        */
      hints->csfile = splitnzero_alloc(optarg, ',');  /* Get file names     */
      break;
    case CLA_SAVEEXT:  /* Extinction array    */
      hints->save.ext = xstrdup(optarg);
      break;
    case CLA_OPACITYFILE:  /* Opacity array file                            */
      hints->f_opa = xstrdup(optarg);
      break;
    /* Set detailed fields files and wavenumbers: */
    case CLA_DETCIA:   /* Detailed CIA            */
      det = &hints->det.cia;
    case CLA_DETTAU:   /* Detailed tau            */
      if(rn == CLA_DETTAU)
        det = &hints->det.tau;
    case CLA_DETEXT:   /* Detailed extinction     */
      if(rn == CLA_DETEXT)
        det = &hints->det.ext;
      /* Get the filename:                                  */
      if (det->n)
        free(det->ref);
      rn = sizeof(det->file);
      if(strlen(optarg)<rn) rn = strlen(optarg);
      /* Get position of colon character, store it in 'i',
         and replace char with NULL:                        */
      i = 0;
      while(i<rn)                /* Find it                 */
        if(optarg[i++] == ':') break;
      optarg[i-1] = '\0';        /* Replace                 */
      strcpy(det->file, optarg); /* Copy name into det.file */
      /* Get wavenumbers, parse values into det.ref and how many in det.n: */
      det->n = getad(0, ',', optarg+i, &det->ref);
      /* Check for errors: */
      if(det->n<1 || i==rn-1) {
        tr_output(TOUT_ERROR | TOUT_BANNER,
                     "Bad format for detailed %s parameter, no valid "
                     "wavenumbers\n", det->name);
        exit(EXIT_FAILURE);
      }
      break;

    case CLA_ETHRESH:    /* Minimum extiction-coefficient threshold */
      hints->ethresh = atof(optarg);
      break;
    case 's':            /* Ray-solution type name     */
      hints->solname = (char *)realloc(hints->solname, strlen(optarg)+1);
      strcpy(hints->solname, optarg);
      break;
    case CLA_ATMOSPHERE: /* Atmosphere file name       */
      hints->f_atm = (char *)realloc(hints->f_atm, strlen(optarg)+1);
      strcpy(hints->f_atm, optarg);
      break;
    case CLA_LINEDB:     /* Line database file name    */
      hints->f_line = (char *)realloc(hints->f_line, strlen(optarg)+1);
      strcpy(hints->f_line, optarg);
      break;
    case CLA_MOLFILE:    /* Known molecular information                     */
      hints->f_molfile = (char *)realloc(hints->f_molfile, strlen(optarg)+1);
      strcpy(hints->f_molfile, optarg);
      break;
    case CLA_OUTSPEC:    /* Output spectrum file name                       */
      //if(hints->f_outspec) free_null(hints->f_outspec);
      hints->f_outspec = (char *)realloc(hints->f_outspec, strlen(optarg)+2);
      strcpy(hints->f_outspec, optarg);
      mkpath(hints->f_outspec, 0777);
      break;
    case CLA_OUTSAMPLE:  /* Sampling output file name  */
      if(hints->f_outsample) free_null(hints->f_outsample);
      hints->f_outsample = (char *)calloc(strlen(optarg)+1, sizeof(char));
      strcpy(hints->f_outsample, optarg);
      mkpath(hints->f_outsample, 0777);
      break;
    case CLA_OUTTOOMUCH: /* toomuch output file name   */
      if(hints->f_toomuch) free_null(hints->f_toomuch);
      if(*optarg != '\0'){
        hints->f_toomuch = (char *)calloc(strlen(optarg)+1, sizeof(char));
        strcpy(hints->f_toomuch, optarg);
        mkpath(hints->f_toomuch, 0777);
      }
      break;
    case CLA_OUTINTENS:  /* Intensity output file name  */
      if(hints->f_outintens) free_null(hints->f_outintens);
      hints->f_outintens = (char *)calloc(strlen(optarg)+2, sizeof(char));
      strcpy(hints->f_outintens, optarg);
      mkpath(hints->f_outintens, 0777);
      break;
    case CLA_SAVEFILES: /* output files with tau, extionction, CIA         */
      isYes = strncmp(optarg, "yes", 3)==0;
      isNo  = strncmp(optarg, "no" , 2)==0;
      if (isYes || isNo){
        hints->savefiles = isYes;
      }
      else{
        printf("\nAllowed arguments for savefiles are: 'yes' or 'no'\n\n");
        exit(EXIT_FAILURE);
      }
      break;
    case CLA_RPRESS:     /* Pressure reference level                        */
      hints->p0 = atof(optarg);
      break;
    case CLA_RRADIUS:    /* Radius reference level                          */
      hints->r0 = atof(optarg);
      break;
    case CLA_GSURF:      /* Surface gravity                                 */
      hints->gsurf = atof(optarg);
      break;

    case CLA_ALLOWQ:
      hints->allowrq = atof(optarg);
      break;
    case CLA_QMOL:
      hints->qmol = xstrdup(optarg);
      break;
    case CLA_QSCALE:
      hints->qscale = xstrdup(optarg);
      break;
    case CLA_OPABREAK: /* Bool: End after opacity calculation               */
      hints->opabreak = 1;
      break;
    case CLA_OPASHARE: /* Bool: Place opacity grid in shared memory         */
      hints->opashare = 1;
      break;

    /* Radius parameters:                                                   */
    case CLA_RADLOW:  /* Lower limit                                        */
      hints->rads.i = atof(optarg);
      break;
    case CLA_RADHIGH: /* Higher limit                                       */
      hints->rads.f = atof(optarg);
      break;
    case CLA_RADDELT: /* Spacing                                            */
      hints->rads.d = atof(optarg);
      break;
    case CLA_RADFCT: /* Units factor                                        */
      hints->rads.fct = atof(optarg);
      break;

    /* Wavelength parameters:                                               */
    case CLA_WAVLOW:    /* Lower limit                                      */
      hints->wavs.i = atof(optarg);
      break;
    case CLA_WAVHIGH:   /* Higher limit                                     */
      hints->wavs.f = atof(optarg);
      break;
    case CLA_WAVFCT:    /* Units factor                                     */
      hints->wavs.fct = atof(optarg);
      break;
    /* Wavenumber parameters:                                               */
    case CLA_WAVNLOW:    /* Lower limit                                     */
      hints->wns.i = atof(optarg);
      break;
    case CLA_WAVNHIGH:   /* Higher limit                                    */
      hints->wns.f = atof(optarg);
      break;
    case CLA_WAVNDELT:   /* Spacing                                         */
      hints->wns.d = atof(optarg);
      hints->wns.n = 0;
      hints->wns.v = NULL;
      break;
    case CLA_WAVNOSAMP:  /* Oversampling                                    */
      hints->wns.o = atof(optarg);
      break;
    case CLA_WNFCT:      /* Units factor                                    */
      hints->wns.fct = atof(optarg);
      break;

    /* Opacity's temperature grid parameters:                  */
    case CLA_TEMPLOW:
      hints->temp.i = atof(optarg);
      break;
    case CLA_TEMPHIGH:
      hints->temp.f = atof(optarg);
      break;
    case CLA_TEMPDELT:
      hints->temp.d = atof(optarg);
      break;

    case 'a':  /* Number of half-widths in a profile:            */
      hints->timesalpha = atof(optarg);
      break;
    case CLA_NDOP:
      hints->nDop = atoi(optarg);
      break;
    case CLA_NLOR:
      hints->nLor = atoi(optarg);
      break;
    case CLA_DMIN:
      hints->dmin = atof(optarg);
      break;
    case CLA_DMAX:
      hints->dmax = atof(optarg);
      break;
    case CLA_LMIN:
      hints->lmin = atof(optarg);
      break;
    case CLA_LMAX:
      hints->lmax = atof(optarg);
      break;

    case 'q':  /* Quiet run:                                                */
      verblevel = 1;
      break;

    case 'v':  /* Increase verbose level:                                   */
      verblevel = strtol(optarg, NULL, 10);
      break;

    case 'V':  /* Print version number and exit:                            */
      if(version_rc>0) snprintf(name, 20, "-rc%i", version_rc);
      else name[0] = '\0';
      printf("This is 'transit' version %i.%i%s\n\n", version, revision, name);
      exit(EXIT_SUCCESS);
      break;

    case 'd':  /* Show defaults. TBD                                        */
      break;
    case '?':
      rn = optopt;
      tr_output(TOUT_ERROR | TOUT_BANNER,
        "Unknown, unsupported, or missing parameter to option of "
        "code %i (%s) passed as argument, use '-h' to see the "
        "available options.\n", rn, (char)rn);
      exit(EXIT_FAILURE);
      break;
    default:   /* Ask for syntax help:                                      */
      tr_output(TOUT_ERROR | TOUT_BANNER,
        "Even though option of code %i (%c) had a valid structure "
        "element, it had no switch control statement.\n", rn, (char)rn);
      exit(EXIT_FAILURE);
      break;
    case 'h':  /* Print out doc-string help:                                */
      prochelp(EXIT_SUCCESS);
      break;

    case CLA_STARRAD:     /* Stellar radius                                 */
      hints->sg.starrad = atof(optarg);
      break;
    case CLA_GORBPAR:     /* Orbital parameters                             */
      getnd(6, ',', optarg,
            &hints->sg.smaxis, &hints->sg.time,  &hints->sg.incl,
            &hints->sg.ecc,    &hints->sg.lnode, &hints->sg.aper);
      break;
    case CLA_GORBPARFCT:  /* Orbital parameters units factors               */
      getnd(6, ',', optarg,
            &hints->sg.smaxisfct, &hints->sg.timefct,  &hints->sg.inclfct,
            &hints->sg.eccfct,    &hints->sg.lnodefct, &hints->sg.aperfct);
      break;
    case CLA_TRANSPARENT: /* Set maximum optical depth to toomuch           */
      hints->sg.transpplanet = 1;
      break;

    case CLA_TOOMUCH:    /* Maximum optical depth to make calculation       */
      hints->toomuch = atof(optarg);
      break;
    case CLA_TAULEVEL:   /* Optical depth integration level                 */
      hints->taulevel = atoi(optarg);
      break;
    case CLA_MODLEVEL:   /* Modulation integration level                    */
      hints->modlevel = atoi(optarg);
      break;

    case CLA_CLOUD:      /* Cloud arguments:                                */
      if(strstr(optarg, "ext") != NULL)
        hints->cl.flag = 1; // Constant extinction
      else if(strstr(optarg, "opa") != NULL)
        hints->cl.flag = 2; // Constant opacity
      optarg += 3;
      if(*optarg != ','  ||  optarg[1] == '\0') {
        tr_output(TOUT_ERROR, "Syntax error in option '--cloud', "
          "parameters need to be given as cloudtype,cloudext,cloudtop,cloudbot.\n");
        exit(EXIT_FAILURE);
      }

      hints->cl.cloudext = strtod(optarg+1, &optarg);
      if(*optarg != ','  ||  optarg[1] == '\0') {
        tr_output(TOUT_ERROR, "Syntax error in option '--cloud', "
          "parameters need to be given as cloudtype,cloudext,cloudtop,cloudbot.\n");
        exit(EXIT_FAILURE);
      }


      hints->cl.cloudtop = strtod(optarg+1, &optarg);
      if(*optarg != ','  ||  optarg[1] == '\0') {
        tr_output(TOUT_ERROR, "Syntax error in option '--cloud', "
          "parameters need to be given as cloudtype,cloudext,cloudtop,cloudbot.\n");
        exit(EXIT_FAILURE);
      }

      hints->cl.cloudbot = strtod(optarg+1, NULL);
      /* Safety check that top > bottom:                                    */
      if(hints->cl.cloudtop > hints->cl.cloudbot) {
        tr_output(TOUT_ERROR, "Syntax error in '--cloud', the cloud top "
                 "(%g) needs to be less than the cloud bottom (%g) .\n", 
                 hints->cl.cloudtop, hints->cl.cloudbot);
        exit(EXIT_FAILURE);
      }
      break;

    case CLA_CLOUDTOP:
      hints->cl.cloudtop = atof(optarg);
      // Ensure cloud extends to high enough pressure to absorb all radiation
      hints->cl.cloudbot = hints->cl.cloudtop + 10;
      hints->cl.cloudext = 100.0; // Large extinction so the cloud is opaque
      hints->cl.flag     = 1; // Constant extinction
      break;

    case CLA_SCATTERING:
      if(!strcmp(optarg, "polar"))
      {
        /* Polarizability method of Villanueva et al. (2022)                */
        hints->scattering_logext = 0.0;
        hints->scattering_flag = 2;
      }
      else
      {
        /* Approximation method of Lecavelier Des Etangs et al. (2008)      */
        hints->scattering_logext = atof(optarg);
        hints->scattering_flag = 1;
      }
      break;

    case CLA_INTENS_GRID:    /* Intensity grid                              */
      hints->angles = xstrdup(optarg);
      break;
    }
  }
  procopt_free();

  return 0;
}


/* \fcnfh
    Initialize transit ray solution sol. and determine if any of
     sol-> name match  hname.
    Return: 0 on success,
           -1 if no matching name was available                  */
int
acceptsoltype(ray_solution **sol,
              char *hname){
  int n=0;

  /* Search through each element in *sol, compare names to hname */
  while (raysols[n]){
    if(strcmp(hname, raysols[n]->name) == 0){
      *sol = (ray_solution *)raysols[n];
      return 0;
    }
    n++;
  }
  return -1;
}


/* \fcnfh
    Set output file names in transit (out, toomuch, and sample).
    Initialize transit.sol.  Set geometry and detailed output
    variables in transit.
    Return: 0 on success                                          */
int
acceptgenhints(struct transit *tr){
  /* Pointer to transithint:                       */
  struct transithint *th = tr->ds.th;

  /* Accept output spectrum file:                  */
  if(th->f_outspec)
    tr->f_outspec = th->f_outspec;
  else{ /* File not specified, use standard output */
    tr->f_outspec = (char *)calloc(2, sizeof(char));
    tr->f_outspec[0] = '-';
  }

  /* Accept toomuch, outsample, and intensity output files:    */
  tr->f_toomuch   = th->f_toomuch;
  tr->f_outsample = th->f_outsample;
  tr->f_outintens = th->f_outintens;
  /* FINDME: Should check if the file exists:                               */
  tr->f_molfile   = th->f_molfile;

  /* Initialize solution type, accept hinted solution if it's in list:      */
  if(acceptsoltype(&tr->sol, th->solname) != 0){
    tr_output(TOUT_ERROR, "Solution kind '%s' is invalid.\n"
      "Currently Accepted are:\n", th->solname);

    ray_solution **sol = (ray_solution **)raysols;
    while(*sol)
      tr_output(TOUT_ERROR, " %s\n", (*sol++)->name);
    exit(EXIT_FAILURE);
    /* FINDME: Fix this error message                                       */
  }

  /* Set hinted geometry hints:                                             */
  setgeomhint(tr);

  /* Accept hints for detailed output:                                      */
  static struct detailout det;
  memcpy(&det, &th->det, sizeof(struct detailout));
  tr->ds.det = &det;

  /* Accept line-profile arguments:                                         */
  /* Check that timesalpha (profile half width in units of Doppler/Lorentz
     broadening half widths) is > 1, and set it's value in transit:         */
  if(th->timesalpha < 1){
    tr_output(TOUT_ERROR,
      "Times of maximum width has to be greater than one: %i\n",
      th->timesalpha);
    return -1;
  }
  tr->timesalpha = th->timesalpha;

  if (th->ethresh <= 0){
    tr_output(TOUT_ERROR,
      "Extinction-coefficient threshold (%.3e) has to be positive.\n",
      th->ethresh);
    return -1;
  }

  /* Pass atmospheric flags into transit struct:                            */
  transitacceptflag(tr->fl, th->fl, TRU_ATMBITS); /* See transit.h          */

  /* Pass flag to break after the opacity-grid calculation:                 */
  tr->opabreak = th->opabreak;

  /* Pass flag to place opacity grid in shared memory:                      */
  tr->opashare = th->opashare;

  /* Set interpolation function flag:                                       */
  switch(tr->fl & TRU_SAMPBITS){
  case TRU_SAMPLIN:
    tr->interpflag = SAMP_LINEAR;
    break;
  case TRU_SAMPSPL:
    tr->interpflag = SAMP_SPLINE;
    break;
  default:
    tr_output(TOUT_ERROR, "Invalid sampling function specified.\n");
    exit(EXIT_FAILURE);
  }
  tr_output(TOUT_DEBUG,
    "transit interpolation flag: %li.\n", tr->interpflag);

  if (th->r0 < 0){
    tr_output(TOUT_ERROR,
      "Reference radius level (%g) must be positive.\n", th->r0);
    return -1;
  }
  tr->r0 = th->r0;

  if (th->p0 < 0){
    tr_output(TOUT_ERROR,
      "Reference pressure level (%g) must be positive.\n", th->p0);
    return -1;
  }
  tr->p0 = th->p0;

  if (th->gsurf < 0){
    tr_output(TOUT_ERROR,
      "Surface gravity (%g cm s^-2) must be positive.\n", th->gsurf);
    return -1;
  }
  tr->gsurf = th->gsurf;

  /* Read in the incident angles for eclipse geometry:                      */
  if (strcmp(tr->sol->name, "eclipse") == 0){
    parseArray(&tr->angles, &tr->ann, th->angles);
    /* FINDME: do some checks that the angles make sense                    */
  }
  if (th->qscale){
   parseArray(&tr->qscale, &tr->nqmol, th->qscale);
   if (countfields(th->qmol, ' ') != tr->nqmol) {
     tr_output(TOUT_ERROR, "qscale (%d) and qmol (%d) should have the "
        "same number of elements.\n", tr->nqmol, countfields(th->qmol,' '));
     exit(EXIT_FAILURE);
   }
  }
  else
    tr->nqmol = 0;

  /* Set cloud structure:                                                   */
  static struct extcloud cl;
  cl.cloudtop = th->cl.cloudtop;
  cl.cloudbot = th->cl.cloudbot;
  cl.cloudext = th->cl.cloudext;  /* Maximum cloud extinction               */
  cl.flag     = th->cl.flag;
  tr->ds.cl = &cl;
  /* Scattering extinction struct         */
  static struct extscat sc;
  sc.logext = th->scattering_logext;
  sc.flag = th->scattering_flag;
  tr->ds.sc = &sc;

  return 0;
}


/*\fcnfh
  Saves hints structure.                                                    */
void
savehint(FILE *out,
         struct transithint *hints){
  /* Save main structure:                                                   */
  fwrite(hints, sizeof(struct transithint), 1, out);

  /* Save strings:                                                          */
  savestr(out, hints->f_atm);
  savestr(out, hints->f_line);
  savestr(out, hints->f_outspec);
  savestr(out, hints->f_outintens);
  savestr(out, hints->f_toomuch);
  savestr(out, hints->f_outsample);
  savestr(out, hints->solname);
  for(int i=0; i<hints->ncross; i++){
    savestr(out, hints->csfile[i]);
  }

  /* Save sub-structures:                                                   */
  savesample_arr(out, &hints->rads);
  savesample_arr(out, &hints->wavs);
  savesample_arr(out, &hints->wns);
  savesample_arr(out, &hints->ips);

}

/* \fcnfh
   Restore hints structure, the structure needs to have been allocated
   before

   @returns 0 on success
            -1 if not all the expected information is read
            -2 if info read is wrong
            -3 if cannot allocate memory
             1 if information read was suspicious                           */
int
resthint(FILE *in,
         struct transithint *hint){
  int rn, res=0;
  /* Restore main structure:                                                */
  rn = fread(hint, sizeof(struct transithint), 1, in);
  if(rn<0)
    return rn;
  else
    res+=rn;

  /* Restore strings:                                                       */
  rn = reststr(in, &hint->f_atm);
  if(rn<0) return rn; else res += rn;
  rn = reststr(in, &hint->f_line);
  if(rn<0) return rn; else res += rn;
  rn = reststr(in, &hint->f_outspec);
  if(rn<0) return rn; else res += rn;
  rn = reststr(in, &hint->f_outintens);
  if(rn<0) return rn; else res += rn;
  rn = reststr(in, &hint->f_toomuch);
  if(rn<0) return rn; else res += rn;
  rn = reststr(in, &hint->f_outsample);
  if(rn<0) return rn; else res += rn;
  rn = reststr(in, &hint->solname);
  if(rn<0) return rn; else res += rn;
  for(int i=0; i<hint->ncross; i++){
    rn = reststr(in, hint->csfile+i);
    if(rn<0) return rn; else res += rn;
  }

  /* Restore sub-structures:                                                */
  restsample_arr(in, &hint->rads);
  if(rn<0) return rn; else res += rn;
  restsample_arr(in, &hint->wavs);
  if(rn<0) return rn; else res += rn;
  restsample_arr(in, &hint->wns);
  if(rn<0) return rn; else res += rn;
  restsample_arr(in, &hint->ips);
  if(rn<0) return rn; else res += rn;

  return res;
}

/* This function ensures that output path directories exist */
int
mkpath(char* path, mode_t mode){
  assert(path && *path);
  char* p;
  for (p=strchr(path+1, '/'); p; p=strchr(p+1, '/')){
    *p='\0';
    if (mkdir(path, mode)==-1){
      if (errno!=EEXIST){
        *p='/'; return -1;
      }
    }
    *p='/';
  }
  return 0;
}


void
printintro(){
  char rcname[20];
  if(version_rc>0)
    snprintf(rcname, 20, "-rc%i", version_rc);
  else
    rcname[0] = '\0';
  fprintf(stdout,
    "--------------------------------------------------\n"
    "                   TRANSIT v%i.%i%s\n"
    "--------------------------------------------------\n",
    version, revision, rcname);
  time_t tim = time(NULL);
  tr_output(TOUT_INFO, "Started on %s\n", ctime(&tim));
}


/* \fcnfh
   Frees hints structure.                                                   */
void
freemem_hints(struct transithint *h){
  /* Free strings which are copied into transit:                            */
  free(h->f_atm);
  free(h->f_line);
  free(h->f_outspec);
  free(h->f_outintens);
  free(h->f_toomuch);
  free(h->f_outsample);
  free(h->f_molfile);

  /* Free other strings:                                                    */
  free(h->solname);
  if (h->ncross){
    free(h->csfile[0]);
    free(h->csfile);
  }

  /* Free sub-structures:                                                   */
  freemem_samp(&h->rads);
  freemem_samp(&h->wavs);
  freemem_samp(&h->wns);
  freemem_samp(&h->ips);
  /* TBD: Free sve if it is ever enabled freesaves(&h->save);               */

  freemem_cloud(&h->cl);
  freemem_detailout(&h->det);
}

void
freemem_cloud(struct extcloud *c){
}

void
freemem_detailout(struct detailout *d){
  freemem_detailfld(&d->ext);
  freemem_detailfld(&d->tau);
  freemem_detailfld(&d->cia);
}

void
freemem_detailfld(struct detailfld *f){
  if (f->n)
    free(f->ref);
}
