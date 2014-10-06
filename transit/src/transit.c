/*
 * transit.c   - Models the modulation produced by a Planet crossing in
 *               front of the star. Main component of the Transit program.
 *
 * Copyright (C) 2003-2006 Patricio Rojo (pato@astro.cornell.edu)
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of version 2 of the GNU General 
 * Public License as published by the Free Software Foundation.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 * 02111-1307, USA.
 */

/* Revision        January 23rd, 2014 Jasmina Blecic
                   implemented eclipse                                     */
/* Revision        March 19th,   2014 Jasmina Blecic
                   implemented switch eclipse/transit                      */
/* Revision        April 26th,   2014 Jasmina Blecic
                   implemented intensity grid and flux                     */

/* TBD: calloc checks */

#include <transit.h>

/* \fcnfh                                                                  */
int main(int argc,      /* Number of variables                             */
         char **argv){  /* Variables                                       */

  /* FINDME: What is s_lt?                                                 */
  /* Initialization of data structure's pointers. Note that s_lt is not
     assigned because the array that is going to point has not been
     initialized yet.                                                      */
  struct transit transit;
  long itr=0;
  struct timeval tv;
  double t0=0.0;

  // /* Intermedio */
  // double press[82] = {3.48800000e+03,   3.08100000e+03,   2.71900000e+03, 2.39800000e+03,
  // 2.11500000e+03,   1.86500000e+03,   1.64500000e+03, 1.45100000e+03,
  // 1.27900000e+03,   1.12700000e+03,   9.93400000e+02, 8.75700000e+02,
  // 7.72100000e+02,   6.81000000e+02,   6.00900000e+02, 5.30600000e+02,   
  // 4.68700000e+02,   4.14200000e+02,   3.66300000e+02, 2.87200000e+02,
  // 2.54500000e+02,   2.00600000e+02,   1.78200000e+02, 1.58500000e+02,
  // 1.24200000e+02,   1.07800000e+02,   7.78700000e+01, 6.48800000e+01,
  // 5.34000000e+01,   4.34600000e+01,   3.50100000e+01, 2.25000000e+01,
  // 1.00000000e+01,   8.50000000e+00,   7.25000000e+00, 6.00000000e+00,
  // 4.75000000e+00,   3.50000000e+00,   2.25000000e+00, 1.00000000e+00,
  // 8.50000000e-01,   7.25000000e-01,   6.00000000e-01, 4.75000000e-01,
  // 3.50000000e-01,   2.25000000e-01,   1.00000000e-01, 8.50000000e-02,
  // 7.25000000e-02,   6.00000000e-02,   4.75000000e-02, 3.50000000e-02,
  // 2.25000000e-02,   1.00000000e-02,   8.50000000e-03, 7.25000000e-03,
  // 6.00000000e-03,   4.75000000e-03,   3.50000000e-03, 2.25000000e-03,
  // 1.00000000e-03,   8.50000000e-04,   7.25000000e-04, 6.00000000e-04,
  // 4.75000000e-04,   3.50000000e-04,   2.25000000e-04, 1.00000000e-04,
  // 8.50000000e-05,   7.25000000e-05,   6.00000000e-05, 4.75000000e-05,
  // 3.50000000e-05,   2.25000000e-05,   1.00000000e-05, 8.50000000e-06,
  // 7.25000000e-06,   6.00000000e-06,   4.75000000e-06, 3.50000000e-06,
  // 2.25000000e-06,   1.00000000e-06};

  // double temp[82] = {5276.19885763,  5135.52171269, 4996.56400992,  4859.67345826,
  //          4725.60624726,  4594.07087687,  4465.56057154, 4339.80360765,
  //          4216.08115684,  4094.74624706,  3976.43759165, 3860.87093829,
  //          3748.13485899,  3638.33043109,  3531.46206944, 3427.72284147,
  //          3326.79043866,  3228.71188471,  3133.60458723, 2951.76341386,
  //          2865.01458616,  2699.92692996,  2620.92622678, 2544.74316773,
  //          2391.74982836,  2306.87359909,  2138.13073955, 2039.83630646,
  //          1969.26677184,  1904.17248806,  1856.12640729, 1791.60938166,
  //          1758.99546   ,  1749.25163068,  1745.06893914, 1738.47329663,
  //          1734.19075206,  1725.65313619,  1716.41297207, 1683.21187516,
  //          1674.13951857,  1665.68063409,  1653.52899661, 1634.39067658,
  //          1605.811824  ,  1547.2910364 ,  1432.0174548 , 1393.44970589,
  //          1377.87891486,  1341.14487078,  1314.76017886, 1269.01704875,
  //          1226.61094637,  1159.9259367 ,  1144.29379242, 1138.40893895,
  //          1123.6064931 ,  1114.91257672,  1096.3489126 , 1081.07115092,
  //          1049.74496698,  1041.49001874,  1038.03736434, 1029.34737827,
  //          1023.3851272 ,  1011.05034207,   998.85438538, 972.77877526,
  //           964.65434658,   962.08540007,   953.38936222, 948.41701059,
  //           935.49498022,   923.74550514,   894.42542711, 885.10559872,
  //           882.4138626 ,   872.13928837,   867.11801587, 852.73758811,
  //           837.63002197,   817.245427};

  // double mu[82] = {2.339,2.339, 2.339, 2.339, 2.339, 2.339, 2.339, 2.339,
  //       2.339, 2.339, 2.339, 2.339, 2.339, 2.339, 2.339, 2.339,
  //       2.339, 2.339, 2.339, 2.339, 2.339, 2.339, 2.339, 2.339,
  //       2.339, 2.339, 2.339, 2.339, 2.339, 2.339, 2.339, 2.339,
  //       2.339, 2.339, 2.339, 2.339, 2.339, 2.339, 2.339, 2.339,
  //       2.339, 2.339, 2.339, 2.339, 2.339, 2.339, 2.339, 2.339,
  //       2.339, 2.339, 2.34, 2.34, 2.34, 2.34, 2.34, 2.339, 2.339,
  //       2.339, 2.339, 2.339, 2.339, 2.338, 2.338, 2.338, 2.338,
  //       2.338, 2.338, 2.338, 2.338, 2.338, 2.338, 2.338, 2.337,
  //       2.337, 2.337, 2.337, 2.336, 2.336, 2.335, 2.335, 2.334,
  //       2.334};
  // double *radio = (double *)calloc(82, sizeof(double));
  // int rp = radpress(939.177, 1.0, 0.0, temp, mu, press, radio, 82, 1e5);
  // transitprint(1,2, "Radius = [");
  // for (rp=0; rp<82; rp++){
  //   transitprint(1,2, "%.2f, ", radio[rp]);
  // }
  // transitprint(1,2,"\b\b]\n");
  // /* Intermedio */
  // return 1;

  memset(&transit, 0, sizeof(struct transit));
  verblevel=2;

  /* Process the command line arguments:                                   */
  fw(processparameters, !=0, argc, argv, &transit);
  t0 = timecheck(verblevel, itr,  0, "processparameters", tv, t0);

  /* Accept all general hints:                                             */
  fw(acceptgenhints, !=0, &transit);
  /* Presentation:                                                         */
  printintro();

  /* No program warnings if verblevel is 0 or 1:                           */
  if(verblevel<2)
    transit_nowarn = 1;

  /* Make wavenumber binning:                                              */
  fw(makewnsample0, <0, &transit);
  t0 = timecheck(verblevel, itr,  1, "makewnsample0", tv, t0);
  if(fw_status>0)
    transitprint(7, verblevel,
                 "makewnsample() modified some of the hinted "
                 "parameters according to returned flag: 0x%lx.\n",
                 fw_status);

  /* Read Atmosphere information:                                          */
  fw(getatm, !=0, &transit);
  t0 = timecheck(verblevel, itr,  2, "getatm", tv, t0);

  /* Read line info:                                                       */
  fw(readlineinfo, !=0, &transit);
  t0 = timecheck(verblevel, itr,  3, "readlineinfo", tv, t0);

  /* Hack: add an index to f_out filename to get different files
     and compare them:                                                     */
  int fout_len = (int)strlen(transit.f_out);  /* Length of tr.f_out        */
  char *dot    = strchr(transit.f_out, '.');  /* Search for '.' char       */
  int dot_pos  = fout_len - (int)strlen(dot); /* Position of dot           */
  char fout[fout_len+1];         /* Copy of tr.f_out                       */
  char str_iter[1];              /* String of iteration number             */
  strcpy(fout, transit.f_out);
  strcpy(fout+dot_pos+1, dot);
  strncpy(fout+dot_pos, "0", 1);

  /* Make radius binning and interpolate data to new value:                */
  fw(makeradsample, <0, &transit);
  t0 = timecheck(verblevel, itr,  4, "makeradsample", tv, t0);
  if(fw_status>0)
    transitprint(7, verblevel,
                 "makeradsample() modified some of the hinted "
                 "parameters according to returned flag: 0x%lx.\n",
                 fw_status);

  /* Calculate opacity grid:                                               */
  fw(opacity, <0, &transit);
  t0 = timecheck(verblevel, itr,  5, "opacity", tv, t0);

  /* Compute sampling of impact parameter:                                 */
  fw(makeipsample, <0, &transit);
  t0 = timecheck(verblevel, itr,  6, "makeipsample", tv, t0);
  if(fw_status>0)
    transitprint(7, verblevel,
                 "makeipsample() modified some of the hinted "
                 "parameters according to returned flag: 0x%lx.\n",
                 fw_status);
 
  /* Print sampling info:                                                  */
  fw(outsample, !=0, &transit);
  t0 = timecheck(verblevel, itr,  7, "outsample", tv, t0);

  /* EDIT: The loop should enclose getatm (getatm will change in the
     future, but for the moment we will leave it as it is).                */
  for (itr=0; itr<1; itr++){
    t0 = timecheck(verblevel, itr,  8, "Start loop", tv, t0);

    /* Initialize CIA:                                                     */
    fw(interpolatecia, !=0, &transit);
    t0 = timecheck(verblevel, itr,  9, "interpolatecia", tv, t0);
 
    /* Compute index of refraction:                                        */
    fw(idxrefrac, !=0, &transit);
    t0 = timecheck(verblevel, itr,  10, "idxrefrac", tv, t0);
 
    /* Calculate extinction coefficient:                                   */
    fw(extwn, !=0, &transit);
    t0 = timecheck(verblevel, itr, 11, "extwn", tv, t0);
 
    sprintf(str_iter, "%li", itr);
    strncpy(fout+dot_pos, str_iter, 1);
    strcpy(transit.f_out, fout);

    /* Ray solutions choice:                                               */
    RaySol path=transit.ds.th->path;

    /* Calculates optical depth for eclipse                                */
    if(path == eclipse){
      transitprint(1,verblevel, "\nCalculating eclipse:\n");

      /* Angle number                                                      */
      struct transithint *th = transit.ds.th;
      long int an = th->ann;

      /* Sets intensity grid:                                              */
      fw(intens_grid, !=0, &transit);
      for(int i = 0; i < an; i++){
        /* Fills out angle index                                           */
        transit.angleIndex = i;

        fw(tau_eclipse, !=0, &transit);
        t0 = timecheck(verblevel, itr, 12, "tau eclipse", tv, t0);
  
        /* Calculates eclipse intensity:                                   */
        /* In cgs units erg/s/sr/cm                                        */
        fw(emergent_intens, !=0, &transit);
        t0 = timecheck(verblevel, itr, 13, "emergent intensity", tv, t0);
      }

      /* Calculates flux  erg/s/cm                                         */
      fw(flux, !=0, &transit);
      t0 = timecheck(verblevel, itr, 14, "flux", tv, t0);

      /* Free no longer needed memory                                      */
      freemem_intensityGrid(transit.ds.intens, &transit.pi);
    }

    /* Calculate optical depth for transit:                                */
    else{
      transitprint(1,verblevel, "\nCalculating transit:\n");
      fw(tau, !=0, &transit);
      t0 = timecheck(verblevel, itr, 12, "tau", tv, t0); 

     /* Calculates transit modulation:                                     */
      fw(modulation, !=0, &transit);
      t0 = timecheck(verblevel, itr, 13, "modulation", tv, t0);
   }
 
    /* Free no longer needed memory                                        */
    freemem_idexrefrac(transit.ds.ir,        &transit.pi);
    freemem_extinction(transit.ds.ex,        &transit.pi);
    freemem_tau(transit.ds.tau,              &transit.pi);

    free(transit.save.ext);
    freemem_cia      (transit.ds.cia, &transit.pi);
    freemem_outputray(transit.ds.out, &transit.pi);
    t0 = timecheck(verblevel, itr, 14, "THE END", tv, t0);
    transitprint(1, verblevel, "----------------------------\n");
  }
  freemem_isotopes(     transit.ds.iso, &transit.pi);
  freemem_molecules(    transit.ds.mol, &transit.pi);
  freemem_atmosphere(   transit.ds.at,  &transit.pi);
  freemem_lineinfotrans(transit.ds.li,  &transit.pi);
  freemem_transit(&transit);

  return EXIT_SUCCESS;
}


/* \fcnfh
   Frees transit structure.                 */
void
freemem_transit(struct transit *tr){
  freemem_hints(tr->ds.th);

  freemem_samp(&tr->rads);
  freemem_samp(&tr->wns);
  freemem_samp(&tr->ips);
  free_atm(&tr->atm);

  free(tr->outpret);
  /* TBD: Free saves once it is enabled
  freemem_saves();                          */
}
