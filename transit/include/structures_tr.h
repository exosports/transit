// Copyright (C) 2015-2016 University of Central Florida. All rights reserved.
// Transit is under an open-source, reproducible-research license (see LICENSE).

#ifndef _TRANSIT_STRUCTURES_H
#define _TRANSIT_STRUCTURES_H

/*  Structures  */

/* Forward declarations:                                                    */
struct transit;
struct geometry;

/* Structure definitions:                                                   */
typedef struct {    /* Sampling struct                                      */
  PREC_NREC n;      /* Number of elements                                   */
  PREC_RES d;       /* Spacing                                              */
  PREC_RES i;       /* Initial value                                        */
  PREC_RES f;       /* Final value                                          */
  int o;            /* Oversampling                                         */
  PREC_RES *v;      /* Values of the sampling                               */
  double fct;       /* v units factor to cgs                                */
} prop_samp;


typedef struct {    /* Isotope's variable (per layer) information:          */
  unsigned int n;   /* Arrays' length                                       */
  double *z;        /* Partition function [radius or temp]                  */
} prop_isov;


typedef struct {    /* Isotope's fixed information:  */
  int d;            /* Database to which they belong */
  char *n;          /* Isotope name                  */
  PREC_ZREC m;      /* Isotope mass                  */
} prop_isof;


typedef struct{    /* Molecule's information: */
  int n;           /* Number of elements      */
  PREC_ATM *d;     /* Density   [n]           */
  PREC_ATM *q;     /* Abundance [n]           */
} prop_mol;


typedef struct {    /* Atmospheric conditions:          */
  double *mm;       /* Mean molecular mass [rad]        */
  PREC_ATM *p;      /* Pressure    [rad]                */
  PREC_ATM *t;      /* Temperature [rad]                */
  PREC_ATM pfct;    /* p units factor to cgs (dyne/cm2) */
  PREC_ATM tfct;    /* t units factor to cgs (Kelvin)   */
} prop_atm;


typedef struct {    /* Database properties:                                 */
  char *n;          /* Database name                                        */
  char *molname;    /* Molecule name                                        */
  unsigned int i;   /* Number of isotopes                                   */
  int s;            /* Cumulative first isotope's index                     */
} prop_db;


typedef struct {    /* One item per database                                */
  unsigned int t;   /* Number of temperatures */ /* FINDME: rename to nt    */
  double *T;        /* Temperatures                                         */
} prop_dbnoext;


typedef struct {             /* Ray solution parameters:                    */
  const char *name;          /* Ray solution name                           */
  const char *file;          /* Ray solution source file                    */
  const short monospace;     /* Request equispaced inpact parameter?        */
  PREC_RES (*optdepth)       /* Extinction-coefficient integrator function  */
       (struct transit *tr,
        PREC_RES b,          /*  Height of ray path                         */
        PREC_RES *ex);       /*  Extinction array [rad]                     */
  PREC_RES (*spectrum)       /* Optical-depth integrator function           */
        (struct transit *tr,
         PREC_RES *tau,      /*  Optical depth                              */
         PREC_RES w,         /*  Wavenumber value                           */
         long last,          /*  index where tau exceeded toomuch           */
         PREC_RES toomuch,   /*  Cutoff optical depth                       */
         prop_samp *r);      /*  Impact parameter or layers' radius         */
} ray_solution;


struct atm_isoprop{    /* Proportional-abundance isotopic parameters: */
  double f;            /* Fractional abundance                        */
  double m;            /* Isotope mass                                */
  int eq;              /* Isotope index from transit.ds.isotopes      */
  char n[maxeisoname]; /* Isotope name                                */
  char t[maxeisoname]; /* Molecule name                               */
};


struct line_transition{  /* Line transition parameters:                     */
  PREC_LNDATA *wl;       /* Wavelength                                      */
  PREC_LNDATA *elow;     /* Lower-state energy                              */
  PREC_LNDATA *gf;       /* gf value                                        */
  short *isoid;          /* Isotope ID (Assumed to be in range)             */
  double wfct;           /* wl units factor to cgs                          */
  double efct;           /* elow units factor to cgs                        */
};


struct lineinfo{             /* Line information parameters:                */
  struct line_transition lt; /* Line transitions                            */
  unsigned short tli_ver;    /* TLI version                                 */
  unsigned short lr_ver;     /* lineread version                            */
  unsigned short lr_rev;     /* lineread revision                           */
  double wi, wf;             /* Initial and final wavelength in database    */
  long endinfo;              /* Position at the end of the info part
                                of the info file                            */
  int ni;                    /* Number of isotopes                          */
  int ndb;                   /* Number of databases                         */
  prop_isov *isov;           /* Variable isotope information (w/temp) [iso] */
  prop_dbnoext *db;          /* Temperature info from databases [DB]        */
  double tmin, tmax;         /* Min and max allowed TLI temperatures        */
  PREC_NREC n_l;             /* Number of lines in database                 */
};


struct atm_data{     /* Atmospheric file parameters:                    */
  int n_aiso;        /* Number of molecules in atmosphere file          */
  prop_samp rads;    /* Radius sampling                                 */
  prop_atm atm;      /* Atmospheric properties                          */
  prop_mol *molec;   /* Molecular information [n_aiso]                  */
  double *mm;        /* Mean molecular mass [rad]                       */
  char *info;        /* Optional atmosphere file information or label   */
  _Bool mass;        /* Abundances in isov by mass (1) of by number (0) */
  int begline;       /* Line of first radius dependent info             */
  long begpos;       /* Position of first radius dependent info         */
};


struct extinction{
  PREC_RES **e;      /* Extinction value [rad][wav]                         */
  int vf;            /* Number of fine-bins of the Voigt function           */
  float ta;          /* Number of alphas that have to be contained in
                        the profile                                         */
  _Bool *computed;   /* Whether the extinction at the given radius was
                        computed [rad]                                      */
  double ethresh;    /* Lower extinction-coefficient threshold              */
};


struct opacityhint{
  int master_PID;       /* Process that will write the shared memory        */
  int num_attached;     /* Count of processes attached to the main segment  */
  volatile long status; /* Flags concerning the state of the shared memory  */
  long Nwave, Ntemp, Nlayer, Nmol; /* Info necessary for sizing shared mem  */
};


struct opacity{
  PREC_RES ****o;         /* Opacity grid [temp][iso][rad][wav]             */
  PREC_VOIGT ***profile;  /* Voigt profiles [nDop][nLor][2*profsize+1]      */
  PREC_NREC **profsize;   /* Half-size of Voigt profiles [nDop][nLor]       */
  double *aDop,           /* Sample of Doppler widths [nDop]                */
         *aLor;           /* Sample of Lorentz widths [nLor]                */
  PREC_RES *temp,         /* Opacity-grid temperature array                 */
           *press,        /* Opacity-grid pressure array                    */
           *wns;          /* Opacity-grid wavenumber array                  */
  PREC_ATM **ziso;        /* Partition function per isotope [niso][Ntemp]   */
  int *molID;             /* Opacity-grid molecule ID array                 */
  long Nwave, Ntemp, Nlayer, Nmol, /* Number of elements in opacity grid    */
      nDop, nLor;         /* Number of Doppler and Lorentz-width samples    */
  int hintID;             /* Shared memory ID of the hint segment           */
  struct opacityhint *hint; /* Information about the shared memory          */
  int mainID;             /* Shared memory ID of the main segment           */
  void *mainaddr;         /* Shared memory address of the main segment      */
};


struct idxref{
  PREC_RES *n;   /* Index of refraction [rad]                               */
};


#if 0
struct savefiles {
  char *ext;   /* saves extinction */
  char *tau;   /* after tau() savefile */
  char *modulation;  /* after modulation() savefile */
};
#endif


struct optdepth{
  PREC_RES **t;     /* Optical depth [wn][ip]                               */
  long *last;       /* Index of the lowest impact parameter, where the
                       optical depth is greater than 'toomuch' (assumed that
                       optical depth increases inwards)  [wn]               */
  double toomuch;   /* Optical depth values greater than this won't be
                       calculated: the extinction is assumed to be zero.    */
};


struct grid{
  PREC_RES **a;      /* Intensity grid, 2D, [an][wnn]                       */
};


struct geometry{
  float smaxis;       /* Semimajor axis                                     */
  double smaxisfct;   /* 'smaxis' times this gives cgs units.               */
  double time;        /* this value is 0 when in the middle of the eclipse  */
  double timefct;     /* 'time' times this gives cgs units                  */
  float incl;         /* inclination of the planetary orbit with respect
                         to the observer, 90 degrees is edge on             */
  float inclfct;      /* Units to convert inclination to radians            */
  double ecc;         /* eccentricty                                        */
  double eccfct;      /* eccentricity's units                               */
  double lnode;       /* longitud of the ascending node                     */
  double lnodefct;    /* longitud of the ascending node units               */
  double aper;        /* argument of the pericenter                         */
  double aperfct;     /* argument of the pericenter units                   */

  double starmass;    /* Mass of the star                                   */
  double starmassfct; /* 'starmass' times this gives cgs units.             */

  double starrad;     /* Star's radius                                      */
  double starradfct;  /* 'starrad' times this gives cgs units.              */

  double x, y;        /* Coordinates of the center of the planet with
                         respect to the star. 'fct' to convert to cgs is
                         found in rads.fct. These fields are not hinted.    */

  _Bool transpplanet; /* If true, set maximum optical depth to toomuch      */
};


struct isotopes{
  prop_isof *isof;    /* Fixed isotope information      [n_i]               */
  prop_isov *isov;    /* Variable isotope information   [n_i]               */
  double *isoratio;   /* Isotopic abundance ratio       [n_i]               */
  int *imol;          /* Molecule index for this isotope[n_i]               */
  prop_db *db;        /* Database's info [n_db]                             */
  int n_db,           /* Number of databases                                */
      n_i,            /* Number of isotopes                                 */
      nmol;           /* Number of different molecules having a line list   */
};

struct molecules{
  int nmol;           /* Number of molecules                                */
  prop_mol *molec;    /* Molecular properties                               */
  char **name;        /* Molecules' names                                   */
  PREC_ZREC *mass;    /* Molecules' masses                                  */
  PREC_ZREC *radius;  /* Molecules' radii                                   */
  PREC_ZREC *pol;     /* Molecules' polarizability                          */
  int *ID;            /* Molecule universal ID                              */
};


struct outputray{
  PREC_RES *o;     /* Output as seen before interaction with telescope      */
};

struct extcloud{
  double cloudext;  /* Maximum opacity in [cm-1] or [cm2/g]                 */
  double cloudtop;  /* Pressure at which clouds start                       */
  double cloudbot;  /* Pressure at which clouds end                         */
  int flag;         /* Sets the cloud model                                 */
};

struct extscat{
  double logext; /* Scattering slope                                        */
  int flag;      /* Sets the scattering model                               */
};

struct saves{
  char *ext;   /* Extinction grid                                           */
};

/* Struct to store requested ext, tau, or cia detailed information:  */
struct detailfld{
  int n;         /* Number of requested wavenumber samples */
  PREC_RES *ref; /* Array of wavenumbers requested         */
  char file[80]; /* Output filename                        */
  char name[30]; /* Name of field                          */
};


struct detailout{
  struct detailfld ext, tau, cia;
};


struct cross{
  int nfiles;         /* Number of cross-section files                      */
  PREC_CS **e;        /* Extinction from all CS sources [nwn][nrad]         */
  PREC_CS ***cs;      /* Tabulated CS extinction     [nfiles][nwave][ntemp] */
  PREC_CS **wn;       /* Tabulated wavenumber arrays [nfiles][nwave]        */
  PREC_CS **temp;     /* Tabulated temperatures      [nfiles][ntemp]        */
  int *nwave;         /* Number of wavenumbers       [nfiles]               */
  int *ntemp;         /* Number of temperatures      [nfiles]               */
  int *nspec;         /* Number of species           [nfiles]               */
  int **mol;          /* Species' ID for each file   [nfiles][2]            */
  double tmin, tmax;  /* CIAs minimum and maximum temperatures              */
};


/* Structure with user hinted data that should go to the 'struct
   transit' upon approval                                                   */
struct transithint{  
  char *f_atm,          /* Atmosphere filename                              */
       *f_line,         /* TLI filename                                     */
       *f_opa,          /* Opacity filename                                 */
       *f_outspec,      /* Output spectrum filename                         */
       *f_toomuch,      /* Output toomuch filename                          */
       *f_outsample,    /* Output sample filename                           */
       *f_outintens,    /* Output intensity filename                        */
       *f_molfile;      /* Known molecular info filename                    */
  PREC_NREC ot;         /* Radius index at which to print output from tau   */
  prop_samp rads, ips,  /* Sampling properties of radius, impact parameter, */
       wavs, wns, temp; /*   wavelength, wavenumber, and temperature        */
  char *angles;         /* String with incident angles (for eclipse)        */
  char *qmol, *qscale;  /* String with species scale factors                */
  float allowrq;        /* How much less than one is accepted, and no warning
                           is issued if abundances don't ad up to that      */
  float timesalpha;     /* Number of alphas that have to be contained in a
                           calculated profile, one side only                */
  int voigtfine;        /* Fine-binning for Voigt function in kapwl(), if
                           accepted it goes to tr.ds.op.vf                  */
  int nDop, nLor;       /* Number of broadening width samples               */
  float dmin, dmax, lmin, lmax; /* Broadening-width samples boundaries      */
  int verbnoise;        /* Noisiest verbose level in a non debugging run    */ 
  _Bool mass;           /* Whether the abundances read by getatm are by
                           mass or number                                   */
  _Bool opabreak;       /* Break after opacity calculation flag             */
  _Bool opashare;       /* Attempt to place opacity grid in shared memory.  */
  long fl;              /* flags                                            */
  _Bool userefraction;  /* Whether to use variable refraction               */
  _Bool savefiles;      /* Whether to save files                            */
  double p0, r0;        /* Pressure and radius reference level              */
  double gsurf;         /* Surface gravity                                  */

  double toomuch;       /* Optical depth values greater than this won't be
                           calculated: the extinction is assumed to be zero */
  int taulevel;         /* Tau integration level of precision               */
  int modlevel;         /* Modulation integration level of precision        */
  char *solname;        /* Name of the type of solution                     */
  struct geometry sg;   /* System geometry                                  */
  struct saves save;    /* Saves indicator of program stats                 */

  struct extcloud cl;
  double scattering_logext;
  int scattering_flag;
  struct detailout det;

  double ethresh;       /* Lower extinction-coefficient threshold           */
  char **csfile;
  int ncross;

};

/* Main data structure                                                      */
struct transit{  
  char *f_atm,       /* Atmosphere filename                                 */
       *f_line,      /* TLI filename                                        */
       *f_opa,       /* Opacity filename                                    */
       *f_outspec,   /* Output spectrum filename                            */
       *f_toomuch,   /* Output toomuch filename                             */
       *f_outsample, /* Output sample filename                              */
       *f_outintens, /* Output intensity filename                           */
       *f_molfile;   /* Known molecular info filename                       */
  PREC_NREC ot;      /* Radius index at which to print output from tau      */

  FILE *fp_atm, *fp_opa, *fp_out, *fp_line; /* Pointers to files            */
  float allowrq;    /* How much less than one is accepted, so that no warning
                       is issued if abundances don't ad up to that          */
  PREC_RES telres;  /* Telescope resolution                                 */
  long int angleIndex; /* Index of the current angle                        */
  prop_samp rads, ips, /* Sampling properties of radius, impact parameter,  */
      owns,            /* oversampled wavenumber,                           */
      wavs, wns, temp; /* wavelength, wavenumber, and temperature           */
  prop_atm atm;      /* Sampled atmospheric data                            */
  _Bool opabreak;    /* Break after opacity calculation                     */
  _Bool opashare;    /* Attempt to place opacity grid in shared memory.     */
  int ndivs,         /* Number of exact divisors of the oversampling factor */
     *odivs;         /* Exact divisors of the oversampling factor           */
  int voigtfine;     /* Number of fine-bins of the Voigt function           */
  float timesalpha;  /* Broadening profile width in number of Doppler or
                        Lorentz half width                                  */
  double p0, r0;     /* Pressure and radius reference level                 */
  double gsurf;      /* Surface gravity                                     */
  int ann;           /* Number of angles                                    */
  double *angles;    /* Array of incident angles for eclipse geometry       */
  int nqmol;         /* Number of species scale factors                     */
  double *qscale;    /* Species scale factors                               */
  int *qmol;         /* Species with scale factors                          */
  int taulevel;      /* Tau integration level of precision                  */
  int modlevel;      /* Modulation integration level of precision           */

  long fl;           /* flags                                               */
  long interpflag;   /* Interpolation flag                                  */
  long pi;           /* progress indicator                                  */

  ray_solution *sol; /* Transit solution type                               */
  PREC_RES *outpret; /* Output dependent on wavelength only as it travels
                        to Earth before telescope                           */

  struct saves save; /* Saves indicator of program stats                    */

  struct {           /* Data structures pointers, this is data that is not
                        required for the final computation                  */
    struct transithint *th;
    struct lineinfo    *li;
    struct atm_data    *at;
    struct extinction  *ex;
    struct opacity     *op;
    struct grid        *intens;
    struct optdepth    *tau;
    struct idxref      *ir;
    struct geometry    *sg;
#if 0
    struct savefiles   *sf;
#endif
    struct isotopes    *iso;
    struct molecules   *mol;
    struct outputray   *out;
    struct extcloud    *cl;
    struct extscat     *sc;
    struct detailout   *det;
    struct cross       *cross;
  }ds;
};

#endif /* _TRANSIT_STRUCTURES_H */
