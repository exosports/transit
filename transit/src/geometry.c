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

#include <transit.h>

/*\fcnfh
  Set transit geometry variables (tr.ds.sg) from hint or default variables
  Return: 0 on success                                                  */
int
setgeomhint(struct transit *tr){
  /* Declare geometry structure and set values to 0: */
  static struct geometry st_sg;
  memset(&st_sg, 0, sizeof(struct geometry));

  struct geometry *sg = tr->ds.sg = &st_sg;
  struct geometry *hg = &tr->ds.th->sg;

  /* Copy hint transpplanet into transit: */
  sg->transpplanet = hg->transpplanet;

  /* If factor values are correctly hinted then use them, otherwise
     assume hard coded values. If you change the following lines, don't
     forget to change the help page in transit.h that the semimajor axis
     is in AU; timing is in hours from mid eclipse; star's mass and
     radius are in solar units.                                          */
  sg->smaxisfct = hg->smaxisfct>0 ? hg->smaxisfct : AU;
  sg->timefct   = hg->timefct>0   ? hg->timefct   : HOUR;
  sg->eccfct    = hg->eccfct>0    ? hg->eccfct    : 1;
  sg->inclfct   = hg->inclfct>0   ? hg->inclfct   : DEGREES;
  sg->aperfct   = hg->aperfct>0   ? hg->aperfct   : DEGREES;
  sg->lnode     = hg->lnodefct>0  ? hg->lnodefct  : DEGREES;
  /* Stellar parameters:                                                 */
  sg->starmassfct = hg->starmassfct>0 ? hg->starmassfct : SUNMASS;
  sg->starradfct  = hg->starradfct>0  ? hg->starradfct  : SUNRADIUS;

  /* If values are not hinted, then use hard coded values:               */
  sg->smaxis = hg->smaxis>0 ? hg->smaxis : AU/sg->smaxisfct;
  sg->time   = hg->time>0   ? hg->time   : 0;
  sg->ecc    = hg->ecc>0    ? hg->ecc    : 0;
  sg->incl   = hg->incl>0   ? hg->incl   : 0;
  sg->aper   = hg->aper>0   ? hg->aper   : 0;
  sg->lnode  = hg->lnode>0  ? hg->lnode  : 0;
  /* Stellar parameters:                                                 */
  sg->starmass = hg->starmass>0 ? hg->starmass:1.101*SUNMASS/sg->starmass;
  sg->starrad  = hg->starrad>0  ? hg->starrad :1.125*SUNRADIUS/sg->starradfct;

  /* Set progressindicator and return:                                   */
  tr->pi |= TRPI_GEOMETRYHINT;
  return 0;
}


/* \fcnfh
   Set X and Y geometry variables

   Return: 0 on success     */
int
setgeom(struct geometry *sg, /* geometry structure                           */
        double time,         /* Time for which to set planet position. If
                                HUGE_VAL, use the default geometry value     */
        long *flags){        /* Progress indicator flag                      */

  transitcheckcalled(*flags, "setgeom", 1, "setgeomhint", TRPI_GEOMETRYHINT);

  /* TBD: include argument of pericenter and longitud of the node in the
     computations, right now both values are used as zero.                */
  /* Auxiliary variables in cgs */
  double smaxis = sg->smaxis*sg->smaxisfct; /* semi-major axis */
  double ecc    = sg->ecc*sg->eccfct;       /* eccentricity    */
  double incl   = sg->incl*sg->inclfct;     /* inclination     */
  // double aper  = sg->aper*sg->aperfct;
  // double lnode = sg->lnode*sg->lnodefct;
  double t = (time<HUGE_VAL ? time:sg->time)*sg->timefct; /* observation time */
  double mass = sg->starmass*sg->starmassfct; /* Stellar mass */

  double prec2 = 1e-12;  /* Square eccentric anomaly precision limit */
  double n = sqrt(GGRAV*mass/(smaxis*smaxis*smaxis)); /* Mean motion */
  double Ea = HUGE_VAL,  /* Approximation to eccentric anomaly       */
         E  = n*time;    /* Eccentric anomaly                        */
  /* Solve Kepler equation to get the eccentric anomaly at time t:   */
  while((E-Ea)*(E-Ea) > prec2){
    Ea = E;
    E  = n*t + ecc*sin(E);
  }

  double M     = smaxis*(1-ecc*ecc);
  double Delta = smaxis*(1-ecc*cos(E));
  double cosv = (M-Delta)/(ecc*Delta);  /* a * cosine of true anomaly */
  double sini = sin(incl);
  double cosi = cos(incl);

  /* Set X and Y: */
  sg->x = Delta*sqrt(cosi*cosi-cosv*cosv);
  sg->y = Delta*sini;

  /* Set progress indicator and return: */
  *flags |= TRPI_GEOMETRY;
  return 0;
}


/* \fcnfh
   Evaluate if x**2 + y**2 < radius**2.
   Any normalized variation of the star, probably limb darkening
   Return: Normalized value
   FINDME: It's not really doing what's described                      */
inline PREC_RES
starvariation(double x,
	      double y,
	      double radius){
  if(x*x + y*y > radius*radius)
    return 0;

  return 1;
}
