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

Copyright (C) 2014 University of Central Florida.  All rights reserved.

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

/* Revision        March 19th,   2014 Jasmina Blecic
                   added eclipse_ray_solution                              */
/* Revision        April 26th,   2014 Jasmina Blecic
                   added new function for ray grid                         */

#include <time.h>
#include <transit.h>
#include <version_tr.h>
#include <math.h>
#include <pu/procopt.h>
#ifdef _USE_GSL
#include <gsl/gsl_version.h>
#endif /* _USE_GSL */

/* Transit ray solution:            */
const static transit_ray_solution *raysols[] = {&slantpath, NULL};
/* slantpath defined in slantpath.c */

/* Eclipse ray solution:            */
const static eclipse_ray_solution *eclipsesols[] = {&eclipsepath, NULL};
/* eclipsepath defined in eclipse.c */


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
    CLA_WAVLOW,
    CLA_WAVHIGH,
    CLA_WAVDELT,
    CLA_WAVOSAMP,
    CLA_WAVMARGIN,
    CLA_WAVNLOW,
    CLA_WAVNHIGH,
    CLA_WAVNDELT,
    CLA_WAVNOSAMP,
    CLA_WAVNMARGIN,
    CLA_ONEPT,
    CLA_ONEABUND,
    CLA_ONEINT,
    CLA_ONEEXTRA,
    CLA_NUMBERQ,
    CLA_ALLOWQ,
    CLA_EXTPERISO,
    CLA_NOEXTPERISO,
    CLA_GORBPAR,
    CLA_GORBPARFCT,
    CLA_GTIME,
    CLA_GTIMEFCT,
    CLA_GMASSRAD,
    CLA_GMASSRADFCT,
    CLA_OUTTAU,
    CLA_TOOMUCH,
    CLA_OUTTOOMUCH,
    CLA_RADFCT,
    CLA_WAVFCT,
    CLA_WNFCT,
    CLA_OUTSAMPLE,
    CLA_TAULEVEL,
    CLA_MODLEVEL,
    CLA_BLOWEX,
    CLA_TAUISO,
    CLA_MINELOW,
    CLA_ETHRESH,
    CLA_CLOUDRAD,
    CLA_CLOUDFCT,
    CLA_CLOUDE,
    CLA_TRANSPARENT,
    CLA_DETEXT,
    CLA_DETCIA,
    CLA_DETTAU,
    CLA_CIAFILE,
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
  };

  /* Generate the command-line option parser: */
  struct optdocs var_docs[]={
    /* General options:                       */
    {NULL,       0,   HELPTITLE,    NULL, NULL,
     "GENERAL OPTIONS:"},
    {"version",  'V', no_argument,  NULL, NULL,
     "Prints version number and exit."},
    {"help",     'h', no_argument,  NULL, NULL,
     "Prints list of possible parameters."},
    {"defaults", 'd', no_argument,  NULL, NULL,
     "Prints default values of the different variable."},
    {"quiet",    'q', no_argument,  NULL, NULL,
     "Decrease the verbose level to the minimum."},
    {"verb",     'v', required_argument,  "2", "Integer",
     "Set the verbosity level."},
    {"paramf",   'p', ADDPARAMFILE, NULL, "filename",
     "Use filename to read parameters in addition to default file(s):"
     " '" DOTCFGFILENM PREPEXTRACFGFILES"'."},

    /* Input and output options:              */
    {NULL,          0,             HELPTITLE,         NULL, NULL,
     "INPUT/OUTPUT OPTIONS:"},
    {"output",     'o',            required_argument, "-",  "outfile",
     "Change output file name, a dash (-) directs to standard output."},
    {"atm",        CLA_ATMOSPHERE, required_argument, "-",  "atmfile",
     "File containing atmospheric info (Radius, pressure, temperature). A dash"
     " (-) indicates alternative input."},
    {"linedb",     CLA_LINEDB,     required_argument, "./res/lineread.tli",
     "linedb", "File containing line information (TLI format, as given by"
               "'lineread'."},
    {"outtoomuch", CLA_OUTTOOMUCH, required_argument, NULL, "filename",
     "Ouputs depth where toomuch optical depth has been attained as a function"
     " of wavelength."},
    {"outsample",  CLA_OUTSAMPLE,  required_argument, NULL, "filename",
     "Outputs sampling information. A dash (-) indicates standard input. By "
     "default there is no such output."},
    {"molfile",    CLA_MOLFILE,    required_argument, "../inputs/molecules.dat",
     "filename", "Path to file with the molecular info."},

    /* Radius options:                                                      */
    {NULL,       0,           HELPTITLE,         NULL, NULL,
     "RADIUS OPTIONS (0-valued defaults would mean to use the values given by "
     "the atmosphere file):"},
    {"raddelt", CLA_RADDELT, required_argument, "-1", "spacing",
     "Radius spacing.  If set, resample the atmospheric layers to a "
     "radius-equidistant sampling array."},
    {"radlow",  CLA_RADLOW,  required_argument, "0", "radius",
     "Lower radius.  If 0, use atmospheric data minimum."},
    {"radhigh", CLA_RADHIGH, required_argument, "0", "radius",
     "Higher radius.  If 0, use atmospheric data maximum."},
    {"radfct",  CLA_RADFCT,  required_argument, "0",  "factor",
     "Radius factor. Multiplicating radius values by this gives centimeters. "
     "If 0, use atmosphere-file factor."},

    /* Atmosphere options:                                                  */
    {NULL,                0,            HELPTITLE,         NULL, NULL,
     "ATMPOSPHERE OPTIONS:"},
    {"numberabund",      CLA_NUMBERQ,  no_argument,       NULL, NULL,
     "Boolean: 0 if the abundances are by number, 1 if by mass."},
    {"allowq",            CLA_ALLOWQ,   required_argument, "0.00001", "value",
     "Maximum allowed cumulative-abundance departure from 1.0."},
    {"refpress",          CLA_RPRESS,   required_argument, NULL, NULL,
     "Reference pressure of the planet's 'surface'."},
    {"refradius",         CLA_RRADIUS,  required_argument, NULL, NULL,
     "Reference altitude of the planet's 'surface'."},
    {"gsurf",             CLA_GSURF,    required_argument, NULL, NULL,
     "Surface gravity in cm/s^2."},

    /* Wavelength options:                                                  */
    {NULL,         0,             HELPTITLE,         NULL,       NULL,
     "WAVELENGTH OPTIONS (all in fct units):"},
    {"wllow",     CLA_WAVLOW,    required_argument, "0",        "wavel",
     "Lower wavelength. 0 if you want to use line data minimum."},
    {"wlhigh",    CLA_WAVHIGH,   required_argument, "0",        "wavel",
     "Upper wavelength. 0 if you want to use line data maximum."},
    {"wldelt",    CLA_WAVDELT,   required_argument, "0.00002",
     "spacing",  "Wavelength spacing. It cannot be 0 or less."},
    {"wlosamp",   CLA_WAVOSAMP,  required_argument, "100",      "integer",
     "Wavelength oversampling. It cannot be 0 or less."},
    {"wlfct",     CLA_WAVFCT,    required_argument, "1",        "factor",
     "Wavelength factor. Multiplicating wavelength values by this gives "
     "centimeters. If 0 or 1 then use centimeters."},
    {"wlmarg",    CLA_WAVMARGIN, required_argument, "0.00000", "boundary",
     "Not trustable range at boundary of line databases. Also transitions "
     "this much away from the requested range will be considered."},

    /* Wavenumber options:                    */
    {NULL,         0,              HELPTITLE,         NULL, NULL,
     "WAVENUMBER OPTIONS (all in cm-1):"},
    //{"wavenumber", 'n',            no_argument,       NULL, NULL,  
    // "Interactively input wavenumber parameters."},
    {"wnlow",     CLA_WAVNLOW,    required_argument, "0",  "waven",
     "Lower wavenumber. 0 if you want to use equivalent of the wavelength "
     "maximum."},
    {"wnhigh",    CLA_WAVNHIGH,   required_argument, "0",  "waven",
     "Upper wavenumber. 0 if you want to use equivalent of the wavelength "
     "minimum."},
    {"wndelt",    CLA_WAVNDELT,   required_argument, "0",  "spacing",
     "Wavenumber spacing. 0 if you want to have the same number of points "
     "as in the wavelength sampling."},
    {"wnosamp",   CLA_WAVNOSAMP,  required_argument, "0",  "integer",
     "Wavenumber oversampling. 0 if you want the same value as for the "
     "wavelengths."},
    {"wnfct",     CLA_WNFCT,      required_argument, "0",  "factor",
     "Output wavenumber factor. Multiplicating wavenumber values by this "
     "gives centimeters. If 0 then use wavelength's value. This only applies "
     "to output, internally wavenumbers will always be in cm-1."},
    {"wnmarg",    CLA_WAVNMARGIN, required_argument, "0",  "boundary",
     "Not trustable range in cm-1 at boundaries. Transitions this much away "
     "from the requested range will be considered. Use the maximum of the "
     "wavelength boundaries if this value is 0."},

    /* Extinction calculation options:        */
    {NULL,         0,               HELPTITLE,         NULL,    NULL,
     "EXTINCTION CALCULATION OPTIONS:"},
    {"finebin",    'f',             required_argument, "5",     "integer",
     "Number of fine-bins to calculate the Voigt function."},
    {"nwidth",     'a',             required_argument, "50",    "number",
     "Number of the max-widths (the greater of Voigt or Doppler widths) that "
     "need to be contained in a calculated profile."},
    {"maxratio",   'u',             required_argument, "0.001", "uncert",
     "Maximum allowed uncertainty in doppler width before recalculating "
     "profile."},
    {"periso",    CLA_EXTPERISO,   no_argument,       NULL,    NULL,
     "Calculate extinction per isotope (allows to display the contribution "
     "from different isotopes, but consumes more memory."},
    {"blowex",     CLA_BLOWEX,      required_argument, "1",     "factor",
     "Blow extinction by factor before computing tau. No physical"
     "significance (use only for debugging)."},
    {"minelow",    CLA_MINELOW,     required_argument, "0",     "low-energy",
     "Lowest limit of low energy to consider (in cm-1)."},
    {"ethreshold", CLA_ETHRESH,   required_argument, "1e-8",    "ethreshold",
     "Minimum extinction-coefficient ratio (w.r.t. maximum in a layer) to "
     "consider in the calculation."},
    {"cloudrad",   CLA_CLOUDRAD,    required_argument, NULL, "radup,raddown",
     "Make a cloud appear linearly from radup to raddown. Use '--cloudfct' "
     "units; if not defined, use radfct."},
    {"cloudfct",   CLA_CLOUDFCT,    required_argument, NULL,    "factor",
     "Cloud radius values specified by '--cloudrad' will be multiplied by "
     "this to convert to cgs units."},
    {"cloudext",   CLA_CLOUDE,      required_argument, NULL,    "extinction",
     "Maximum extinction of the cloud, which opacity will linearly increase "
     "from 'radup' to 'raddown'."},
    {"detailext",  CLA_DETEXT,      required_argument, NULL,
     "filename:wn1,wn2,...",
     "Save extinction at specified wavenumbers in filename."},
    {"detailcia",  CLA_DETCIA,      required_argument, NULL,
     "filename:wn1,wn2,...",
     "Save extinction due to CIA at specified wavenumbers in filename."},
    {"cia",        CLA_CIAFILE,     required_argument, NULL,    "filenames",
     "Use the indicated filenames for CIA opacities, it is a comma"
     "separated list."},
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

    /* Resulting ray options:                 */
    {NULL,        0,            HELPTITLE,         NULL, NULL,
     "RESULTING RAY OPTIONS:"},
    {"solution",  's',          required_argument, "Slant Path", "sol_name",
     "Name of the kind of output solution ('slant path'is currently the only "
     "availabale alternative)."},
    {"toomuch",   CLA_TOOMUCH,  required_argument, "20", "optdepth",
     "If optical depth for a particular path is larger than optdepth, then do "
     "not proceed to lower radius."},
    {"tauiso",    CLA_TAUISO,   required_argument, "0",  "isoid",
     "Compute tau only for isotope indexed in isoid (index which can actually "
     "be different from what you expect)."},
    /* FINDME: Explain me */
    {"outtau",    CLA_OUTTAU,   required_argument, "0",  "#radius",
     "Output is optical depth instead of modulation. It will be asked which "
     "radius to plot."},
    {"taulevel",  CLA_TAULEVEL, required_argument, "1",  "integer",
     "Calculate the lightray path with a constant (1) or variable (2) "
     "index of refraction."},
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
    {"solutiontype", CLA_SOLUTION_TYPE, required_argument, "eclipse",
     NULL, "Ray solution type (eclipse or transit)."},
    {"raygrid",      CLA_INTENS_GRID, required_argument, "0, 20, 40, 60, 80",
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
  /* TBD: have this option user selectable */
  st_trh.tauiso = 0;

  int rn,  /* optdocs' short option */ 
      i;   /* Auxilliary index      */
  prop_samp *samp = NULL;
  char name[20], rc, *lp; /* FINDME */
  /* For interactive rad, wavelength, and wavenumber sample inputs: */
  char *sampv[] = {"Initial", "Final", "Spacing", "Oversampling integer for"};
  double rf; /* FINDME */

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
    transitDEBUG(21, verblevel, "Processing option '%c', argum: %s\n",
                 rn, optarg);
    /* Handle each case: */
    switch(rn){
    case CLA_CIAFILE:
      hints->ncia    = nchar(optarg, ',') + 1;        /* Count files */
      hints->ciafile = splitnzero_alloc(optarg, ','); /* Get files   */
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
      if(det->n<1 || i==rn-1)
        transiterror(TERR_SERIOUS,
                     "Bad format for detailed %s parameter, no valid "
                     "wavenumbers\n", det->name);
      break;

    case CLA_MINELOW:    /* Lowest limit of low energy */
      hints->minelow = atof(optarg);
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
    case 'o':            /* Output file name           */
      hints->f_out = (char *)realloc(hints->f_out, strlen(optarg)+2);
      strcpy(hints->f_out, optarg);
      break;
    case CLA_OUTSAMPLE:  /* Sampling output file name  */
      if(hints->f_outsample) free_null(hints->f_outsample);
      hints->f_outsample = (char *)calloc(strlen(optarg)+1, sizeof(char));
      strcpy(hints->f_outsample, optarg);
      break;
    case CLA_OUTTOOMUCH: /* toomuch output file name   */
      if(hints->f_toomuch) free_null(hints->f_toomuch);
      if(*optarg != '\0'){
        hints->f_toomuch = (char *)calloc(strlen(optarg)+1, sizeof(char));
        strcpy(hints->f_toomuch, optarg);
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

    case CLA_OUTTAU:
      if(atoi(optarg))
        hints->fl |= TRU_OUTTAU;
      hints->ot = atoi(optarg)-1;
      break;
    /* Interactive radius, wavelength, and wavenumber inputs: */
    case 'r':
      samp = &hints->rads;
      fprintpad(stderr, 1, "In units of planetary radius ...\n");
      strcpy(name, "radius");
    case 'w':
      if(rn=='w'){
        samp = &hints->wavs;
        /* FINDME: should it be CGS units? */
        fprintpad(stderr, 1, "In nanometers ...\n");
        strcpy(name, "wavelength");
      }
    case 'n':
      if(rn=='n'){
        samp = &hints->wns;
        fprintpad(stderr, 1, "In cm-1 ...\n");
        strcpy(name, "wavenumber");
      }
      /* Ask for initial, final, spacing, and oversampling: */
      for(i=0; i<4; i++){
        if(i==3 && samp==&hints->rads) /* Skip radius oversampling */
          break;
        while(rn){
          /* Read from standard input (screen): */
          fprintf(stderr, "- %s %s: ", sampv[i], name);
          /* FINDME: This will probably break */
          samp->i = readd(stdin, &rc);
          if(!rc)
            break;
          switch(rc){
            /*
              case '?':
              fprintpad(stderr,1,"%s\n",var_docs[longidx].doc);
            */
          default:
            fprintf(stderr, "Try again.\n");
          }
        }
        rn = 1;
        /*        longidx++;*/
      }
      break;

    case CLA_ALLOWQ:
      hints->allowrq = atof(optarg);
      break;
    /* FINDME: This is not even an option! */
    case CLA_NUMBERQ: /* Bool: abundances by number (0), or by mass (1) */
      hints->mass = 0;
      break;
    case CLA_ONEPT: /* Single P-T calculation: */
      /* getnd defined in /pu/src/iomisc.c */
      if((rn=getnd(3, ',', optarg, &hints->onept.p, &hints->onept.t, &rf))!=3){
        if(rn>0)
          transiterror(TERR_SERIOUS,
                       "At least one of the values given for the floats "
                       "pressure (%g), temperature (%g), or integer number of "
                       "extra isotopes (%g), was not a correct value.\n",
                       hints->onept.p, hints->onept.t, rf);
        else
          transiterror(TERR_SERIOUS,
                       "There were %i comma-separated fields instead of 3 \n"
                       "for '--onept' option", -rn);
      }
      hints->onept.ne = (int)rf; /* Number of extra isotopes: */
      if(rf != hints->onept.ne)  /* rf is not an integer      */
        transiterror(TERR_SERIOUS,
                     "A non-integer (%g) number of extra isotopes was given "
                     "with the option --onept\n", rf);
      hints->onept.one = 1;
      break;

    case CLA_ONEABUND:
      /* Get number and values of input abundances: */
      if((hints->onept.nq = getad(0, ',', optarg, &hints->onept.q))<1)
        transiterror(TERR_SERIOUS,
                     "None of the given isotope abundances were "
                     "accepted %s\n", optarg);
      transitprint(2, verblevel, "%i abundance isotopes were correctly "
                                 "given: %s\n", hints->onept.nq, optarg);
      break;

    case CLA_ONEEXTRA:
      /* Count instances of ',' and replace with '\0': */
      rn = ncharchg(optarg, ',', '\0')+1;
      if(hints->onept.n){
        free(hints->onept.n);
        free(hints->onept.n[0]);
        free(hints->onept.m);
      }
      /* Allocate pointer to array of extra isotopes names and masses: */
      hints->onept.n    = (char     **)calloc(rn,            sizeof(char *));
      /* FINDME: Where maxeisoname comes from? */
      hints->onept.n[0] = (char      *)calloc(rn*maxeisoname,sizeof(char));
      hints->onept.m    = (PREC_ZREC *)calloc(rn,            sizeof(PREC_ZREC));
      for(i=0; i<rn; i++){
        /* Set pointer at the correct position to write the name: */
        hints->onept.n[i] = hints->onept.n[0] + i*maxeisoname;
        /* Get mass, set pointer lp right after:                  */
        hints->onept.m[i] = strtod(optarg, &lp);
        if(lp==optarg)
          transiterror(TERR_SERIOUS,
                       "Bad format in the field #%i of --oneextra. It doesn't "
                       "have a valid value for mass. The field should be "
                       "<mass1><name1> with only an optional dash betweeen "
                       "the mass and name:\n %s\n", i+1, optarg);
        /* If there is a dash, ignore it:   */
        if(*lp=='-') lp++;
        /* Get name:                        */
        strncpy(hints->onept.n[i], lp, maxeisoname);
        optarg = lp;
        /* Set last name character to '\0': */
        hints->onept.n[i][maxeisoname-1] = '\0';
        /* FINDME: Explain this while: */
        while(*optarg++);
        if(lp==optarg)
          transiterror(TERR_SERIOUS,
                       "Bad format in the field #%i of --oneextra. It doesn't "
                       "have a valid isotope name. The field should be "
                       "<mass1><name1> with only an optional dash betweeen "
                       "the mass and name:\n %s\n", i+1, optarg);
      }
      /* Number of extra mass-name isotope pairs: */
      hints->onept.nm = rn;
      break;

    case CLA_ONEINT:
      /* Set flags ... FINDME: try to understand */
      hints->fl = (hints->fl & ~TRU_ATM1PBITS) | TRU_ATMASK1P;
      break;
    /* Radius parameters:             */
    case CLA_RADLOW:  /* Lower limit  */
      hints->rads.i = atof(optarg);
      break;
    case CLA_RADHIGH: /* Higher limit */
      hints->rads.f = atof(optarg);
      break;
    case CLA_RADDELT: /* Spacing      */
      hints->rads.d = atof(optarg);
      break;
    case CLA_RADFCT: /* Units factor  */
      hints->rads.fct = atof(optarg);
      break;
    /* Wavelength parameters:                        */
    case CLA_WAVLOW:    /* Lower limit               */
      hints->wavs.i = atof(optarg);
      break;
    case CLA_WAVHIGH:   /* Higher limit              */
      hints->wavs.f = atof(optarg);
      break;
    case CLA_WAVDELT:  /* Spacing                    */
      hints->wavs.d = atof(optarg);
      if(hints->wavs.d<=0)
        transiterror(TERR_SERIOUS,
                     "Wavelength spacing has to be greater than zero, "
                     "instead of %g.\n", hints->wavs.d);
      hints->wavs.n = 0;
      hints->wavs.v = NULL;
      break;
    case CLA_WAVFCT:    /* Units factor              */
      hints->wavs.fct = atof(optarg);
      break;
    case CLA_WAVOSAMP:  /* Oversampling              */
      hints->wavs.o = atof(optarg);
      break;
    case CLA_WAVMARGIN: /* Trustable margin at edges */
      hints->margin = atof(optarg);
      break;
    /* Wavenumber parameters:                 */
    case CLA_WAVNLOW:    /* Lower limit       */
      hints->wns.i = atof(optarg);
      break;
    case CLA_WAVNHIGH:   /* Higher limit      */
      hints->wns.f = atof(optarg);
      break;
    case CLA_WAVNDELT:   /* Spacing           */
      hints->wns.d = atof(optarg);
      hints->wns.n = 0;
      hints->wns.v = NULL;
      break;
    case CLA_WAVNOSAMP:  /* Oversampling      */
      hints->wns.o = atof(optarg);
      break;
    case CLA_WAVNMARGIN: /* Wavenumber margin */
      hints->wnm = atof(optarg);
      break;
    case CLA_WNFCT:      /* Units factor      */
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

    case 'u':  /* Maximum accepted Doppler-width ratio         */
      hints->maxratio_doppler = atof(optarg);
      break;
    case 'f':  /* Number of fine-binning for Voigt calculation */
      hints->voigtfine = atoi(optarg);
      break;
    case 'a':  /* Number of half-widths in a profile            */
      hints->timesalpha = atof(optarg);
      break;

    case 'q':  /* Quiet run                                    */
      verblevel = 0;
      break;

    case 'v':  /* Increase verbose level                       */
      verblevel = strtol(optarg, NULL, 10);
      break;

    case 'V': /* Print version number and exit                 */
      if(version_rc>0) snprintf(name, 20, "-rc%i", version_rc);
      else name[0] = '\0';
      printf("This is 'transit' version %i.%i%s\n\n", version, revision, name);
      exit(EXIT_SUCCESS);
      break;

    case 'd':  /* Show defaults. TBD */
      break;
    case '?':
      rn = optopt;
      transiterror(TERR_SERIOUS,
                   "Unknown, unsupported, or missing parameter to option of "
                   "code %i (%s) passed as argument, use '-h' to see "
                   "accepted options.\n", rn, (char)rn);
      break;
    default:   /* Ask for syntax help */
      transiterror(TERR_CRITICAL,
                   "Even though option of code %i(%c) had a valid structure "
                   "element, it had no switch control statement. Code need to "
                   "be revised.\n", rn, (char)rn);
      break;
    case 'h': /* Print out doc-string help */
      prochelp(EXIT_SUCCESS);
      break;
    case CLA_EXTPERISO:   /* Calculate extinction per isotope     */
      hints->fl |= TRU_EXTINPERISO;
      break;
    case CLA_NOEXTPERISO: /* FINDME: Redundant, deprecate         */
      hints->fl &= ~TRU_EXTINPERISO;
      break;
    case CLA_BLOWEX:      /* Multiplicating factor for extinction */
      hints->blowex = atof(optarg);
      break;

    case CLA_STARRAD:     /* Stellar radius                       */
      hints->sg.starrad = atof(optarg);
      break;
    case CLA_GORBPAR:     /* Orbital parameters                   */
      getnd(6, ',', optarg,
            &hints->sg.smaxis, &hints->sg.time,  &hints->sg.incl, 
            &hints->sg.ecc,    &hints->sg.lnode, &hints->sg.aper);
      break;
    case CLA_GORBPARFCT:  /* Orbital parameters units factors     */
      getnd(6, ',', optarg,
            &hints->sg.smaxisfct, &hints->sg.timefct,  &hints->sg.inclfct,
            &hints->sg.eccfct,    &hints->sg.lnodefct, &hints->sg.aperfct);
      break;
    case CLA_TRANSPARENT: /* Set maximum optical depth to toomuch */
      hints->sg.transpplanet = 1;
      break;

    case CLA_TOOMUCH:    /* Maximum optical depth to make calculation */
      hints->toomuch = atof(optarg);
      break;
    case CLA_TAUISO:     /* Single isotope index to calculate tau */
      hints->tauiso = atoi(optarg);
      break;
    case CLA_TAULEVEL:   /* Optical depth integration level       */
      hints->taulevel = atoi(optarg);
      break;
    case CLA_MODLEVEL:   /* Modulation integration level          */
      hints->modlevel = atoi(optarg);
      break;

    case CLA_CLOUDRAD:   /* Lower and higher limits of a cloud    */
      hints->cl.rini = strtod(optarg, &optarg);
      if(*optarg!=',' || optarg[1]=='\0')
        transiterror(TERR_SERIOUS,
                     "Syntax error in option '--cloudrad', parameters need "
                     "to be radup,raddown.\n");
      hints->cl.rfin = strtod(optarg+1, NULL);
      if(hints->cl.rini<hints->cl.rfin ||
         (hints->cl.rfin<=0 && hints->cl.rini!=0))
        transiterror(TERR_SERIOUS,
                     "Syntax error in option '--cloudrad', radup(%g) needs to "
                     "be bigger than raddown (%g) and both greater than "
                     "zero.\n", hints->cl.rini, hints->cl.rfin);
      break;
    case CLA_CLOUDFCT:  /* Cloud limits units factor              */
      hints->cl.rfct = atof(optarg);
      break;
    case CLA_CLOUDE:    /* Maximum cloud opacity                  */
      hints->cl.maxe = atof(optarg);
      break;
    case CLA_SOLUTION_TYPE:     /* Ray-solution type              */
      hints->path = eclipse;
      if(strncasecmp(optarg,"transit",7)==0)
        hints->path = transit;
      break;
    case CLA_INTENS_GRID:      /* Intensity grid                  */
      parseAngles(hints, optarg);
      break;
    }
  }
  procopt_free();

  return 0;
}

/* FINDME: Move to pu */
/* \fcnfh
   Parser for angle parameter from configuration file
   Reads the angles and write them in the transithint variable
   Returns: zero on success                                      */

int
parseAngles(struct transithint *hints, char *slant){
  /* Declares a copy variable of the text from the
     configuration file so we can work on a copy                */
  char angles_copy[100];

  /* Defines a delimiter */
  const char s[2] = ",";

  /* Takes characters (tokens) between delimiters               */
  char *token;

  /* Declares floating variable for strings                     */
  double theta;

  /* Declares variable for number of angles                     */
  long int angNum = 0;

  /* Declares result of the conversion from string to floats    */
  int res;

  /* Makes a copy of angles read from cfg file                  */
  strcpy(angles_copy, slant);

  /* Gets first token                                           */
  token = strtok(angles_copy, s);

  /* Walks through other tokens                                 */
  while(token != NULL){
    /* Converts strings to floats                              */
    res = sscanf(token, "%lf", &theta);

    if(res==1){
      /* Fills out hints variable angles                       */
      hints->angles[angNum] = theta;
      angNum++;
    }
    token = strtok(NULL, s);
  }
  /* Fill out hints variable ann                                */
  hints->ann = angNum;

  return 0;
}


/* \fcnfh
    Initialize transit ray solution sol. and determine if any of
     sol-> name match  hname.
    Return: 0 on success,
           -1 if no matching name was available                  */
int
acceptsoltype(transit_ray_solution **sol, 
              char *hname){
  /* raysols is defined at the beginning of this file            */
  *sol = (transit_ray_solution *)raysols[0];
  int len;  /* Length of hname                                   */
  len = strlen(hname);

  /* Search through each element in *sol, compare names to hname */
  while(*sol){
    if(strncasecmp(hname, (*sol)->name, len)==0)
      return 0;
    sol++;
  }

  return -1;
}


/* \fcnfh
    Initialize eclipse ray solution solution and determine
    if any of ecl-> name match hname.
    Return: 0 on success,
           -1 if no matching name was available                  */
int
accepteclipsetype(eclipse_ray_solution **ecl,
                  char *hname){
  *ecl = (eclipse_ray_solution *)eclipsesols[0];
  int len;

  len = strlen(hname);

  while(*ecl){
    if(strncasecmp(hname,(*ecl)->name,len)==0)
      return 0;
    ecl++;
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

  /* Accept output file:                           */
  if(th->f_out)
    tr->f_out = th->f_out;
  else{ /* File not specified, use standard output */
    tr->f_out    = (char *)calloc(2, sizeof(char));
    tr->f_out[0] = '-';
  }

  /* Accept toomuch and outsample output files:    */
  tr->f_toomuch   = th->f_toomuch;
  tr->f_outsample = th->f_outsample;
  /* FINDME: Should check if the file exists:                               */
  tr->f_molfile   = th->f_molfile;

  /* Initialize solution-type, accept hinted ray-solution
     if it's name exists:                          */
  int noSolName = acceptsoltype(&tr->sol,     th->solname) != 0;
  int noEclName = accepteclipsetype(&tr->ecl, th->solname) != 0;
  if(noSolName && noEclName){
    transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
                 "Solution kind '%s' is invalid.\n"
                 "Currently Accepted are:\n", th->solname);

    transit_ray_solution **sol = (transit_ray_solution **)raysols;
    while(*sol)
      transiterror(TERR_SERIOUS | TERR_NOPREAMBLE | TERR_ALLOWCONT,
                   " %s\n", (*sol++)->name);

    eclipse_ray_solution **ecl = (eclipse_ray_solution **)eclipsesols;
    while(*ecl)
      transiterror(TERR_SERIOUS | TERR_NOPREAMBLE | TERR_ALLOWCONT,
                   " %s\n", (*ecl++)->name);
    exit(EXIT_FAILURE);
  }

  /* Set hinted geometry hints:                    */
  setgeomhint(tr);

  /* Accept hints for detailed output:             */
  static struct detailout det;
  memcpy(&det, &th->det, sizeof(struct detailout));
  tr->ds.det = &det;

  /* Accept line-profile arguments:                                         */

  /* Check that voigtfine (oversampling factor of the profile binning)
     is > 1, and set it's value in transit:                                 */
  if(th->voigtfine < 1){
    transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
                 "Fine binning of Voigt function has to be positive: %i.\n",
                 th->voigtfine);
    return -1;
  }
  tr->voigtfine = th->voigtfine;

  /* Check that timesalpha (profile half width in units of Doppler/Lorentz
     broadening half widths) is > 1, and set it's value in transit:         */
  if(th->timesalpha < 1){
    transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
                 "Times of maximum width has to be greater than one: %i\n",
                 th->voigtfine);
    return -1;
  }
  tr->timesalpha = th->timesalpha;

  /* maxratio is the maximum allowed ratio change before recalculating
     profile array:                                                         */
  if(th->maxratio_doppler<0){
    transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
                 "Maximum allowed Doppler width ratio change. Has to "
                 "be 0 or positive (%g).\n", th->maxratio_doppler);
    return -1;
  }

  if(th->minelow < 0){
    transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
                 "Minimum line-transition low energy limit has to be positive:"
                 " %i\n", th->minelow);
    return -1;
  }
  tr->minelow = th->minelow;

  if (th->ethresh <= 0){
    transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
                 "Extinction-coefficient threshold (%.3e) has to be "
                 "positive.\n", th->ethresh);
    return -1;
  }

  /* Pass atmospheric flags into transit struct:                            */
  transitacceptflag(tr->fl, th->fl, TRU_ATMBITS); /* See transit.h          */

  /* Set interpolation function flag:                                       */
  switch(tr->fl & TRU_SAMPBITS){
  case TRU_SAMPLIN:
    tr->interpflag = SAMP_LINEAR;
    break;
  case TRU_SAMPSPL:
    tr->interpflag = SAMP_SPLINE;
    break;
  default:
    transiterror(TERR_SERIOUS, "Invalid sampling function specified.\n");
  }
  transitprint(1, verblevel, "transit interpolation flag: %li.\n",
                             tr->interpflag);

  if (th->r0 < 0){
    transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
                 "Reference radius level (%g) must be positive.\n", th->r0);
    return -1;
  }
  tr->r0 = th->r0;

  if (th->p0 < 0){
    transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
                 "Reference pressure level (%g) must be positive.\n", th->p0);
    return -1;
  }
  tr->p0 = th->p0;

  if (th->gsurf < 0){
    transiterror(TERR_SERIOUS|TERR_ALLOWCONT,
                "Surface gravity (%g cm s^-2) must be positive.\n", th->gsurf);
    return -1;
  }
  tr->gsurf = th->gsurf;
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
  savestr(out, hints->f_out);
  savestr(out, hints->f_toomuch);
  savestr(out, hints->f_outsample);
  savestr(out, hints->solname);
  for(int i=0; i<hints->ncia; i++){
    savestr(out, hints->ciafile[i]);
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
  rn = reststr(in, &hint->f_out);
  if(rn<0) return rn; else res += rn;
  rn = reststr(in, &hint->f_toomuch);
  if(rn<0) return rn; else res += rn;
  rn = reststr(in, &hint->f_outsample);
  if(rn<0) return rn; else res += rn;
  rn = reststr(in, &hint->solname);
  if(rn<0) return rn; else res += rn;
  for(int i=0; i<hint->ncia; i++){
    rn = reststr(in, hint->ciafile+i);
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


void
printintro(){
  char rcname[20];
  if(version_rc>0)
    snprintf(rcname, 20, "-rc%i", version_rc);
  else
    rcname[0] = '\0';
  transitprint(1, verblevel,
               "-----------------------------------------------\n"
               "                TRANSIT v%i.%i%s\n"
               "-----------------------------------------------\n",
               version, revision, rcname);
  time_t tim = time(NULL);
  transitprint(2, verblevel, "Started on %s\n", ctime(&tim));
}


/* \fcnfh
   Frees hints structure.                                                   */
void
freemem_hints(struct transithint *h){
  /* Free strings which are copied into transit:                            */
  free(h->f_atm);
  free(h->f_line);
  free(h->f_out);
  free(h->f_toomuch);
  free(h->f_outsample);
  free(h->f_molfile);

  /* Free other strings:                                                    */
  free(h->solname);
  if (h->ncia){
    free(h->ciafile[0]);
    free(h->ciafile);
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
