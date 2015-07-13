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

Joseph Harrington, Patricio Cubillos, and Jasmina Blecic
UCF PSB 441
4111 Libra Drive
Orlando, FL 32816-2385
USA

Thank you for using transit!
******************************* END LICENSE ******************************/

#include <spline.h>

/* FUNCTION
   Calculates the z array for cubic splines. Called by splinterp and
   spline_init                                                              */
inline double *
tri(double *a,
    double *d,
    double *c,
    double *b,
    double *e,
    long n){

  int i, j;
  double xmult;

  for (i=1; i<=n; i++){
    xmult = a[i-1]/d[i-1];
    d[i] = d[i] - xmult*c[i-1];
    b[i] = b[i] - xmult*b[i-1];
  }

  e[n-1] = b[n-1]/d[n-1];

  for (j=n; j >= 1; j--){
    e[j] = (b[j-1]-c[j-1]*e[j+1])/d[j-1];
  }

  e[0] = 0;
  e[n+1] = 0;
  return e;
}


/* FUNCTION
   Calculate the y values for a cubic spline                                */
inline double *
spline3(double *xi,
        double *yi,
        double *x,
        double *z,
        double *h,
        double *y,
        long nx,
        long N){

  int i, j=0, n;  /* Indices                                                  */
  double dx;    /* Distance from interpolation x-coordinate to
                   previous x-coordinate in the source array                */
  double amult, bmult, cmult;  /* Coefficients of the cubic polynomial      */

  y[0] = yi[0];

  /* Calculate the spline interpolation y-value for each x-coordinate
     in the requested array                                                 */
  for (n=0; n<=nx; n++){
    for (i=0; i<N-1; i++){
      if (xi[i] <= x[n]){
        j = i;
      }
    }
    amult = (z[j+1] - z[j]) / (6*h[j]);
    bmult = z[j] / 2;
    cmult = (yi[j+1] - yi[j]) / h[j] - h[j]/6 * (z[j+1] + 2 * z[j]);
    dx = x[n] - xi[j];
    /* Calculate y-value from cubic polynomial:                             */
    y[n] = yi[j] + dx*cmult + dx*dx*bmult + dx*dx*dx*amult;
  }

  return y;
}


/* FUNCTION
   Cubic spline interpolation.  Given points described by arrays
   xi and yi, interpolates to coordinates of xout and puts them in yout.   */
double *
splinterp(long N,         /* Length of xi                                  */
          double *xi,
          double *yi,
          long nx,
          double *xout,
          double *yout){

  double *b;  /* dy/dx                                                      */
  double *h;  /* Spacing between x-coordinates                              */
  double *k;  /* Array related to derivative of spline                      */
  double *a;  /* Shifted h array                                            */
  double *d;  /* Midstep array used in calculations                         */
  double *z;  /* Array created by tri() function                            */
  double *c;  /* FINDME: unnecessary. Remove                                */
  double *y;  /* Array created by spline3() function                        */
  double *e;
  int i;

  /* Account for the endpoints. The code is written to calculate nx
     points not including the final endpoint. This line makes it so
     the return has the desired number of values                            */
  nx -= 1;

  /* Allocate all arrays to be used:                                        */
  b = calloc(N,    sizeof(double));
  h = calloc(N,    sizeof(double));
  k = calloc(N,    sizeof(double));
  a = calloc(N,    sizeof(double));
  d = calloc(N,    sizeof(double));
  z = calloc(N+1,  sizeof(double));
  c = calloc(N,    sizeof(double));
  e = calloc(N,    sizeof(double));
  y = calloc(nx+1, sizeof(double));


  /* Calculate first entry of array                                         */
  b[0] = (yi[1] - yi[0]) / (xi[1] - xi[0]);

  /* The following loops fill out the arrays                                */
  for (i=0; i<=N-2; i++){
    b[i+1] = -b[i] + 2*(yi[i+1] - yi[i]) / (xi[i+1] - xi[i]);
  }

  for (i=0; i<=N-2; i++){
    h[i] = xi[i+1] - xi[i];
  }

  for (i=0; i<=N-3; i++){
    k[i] = 6*( (yi[i+2] - yi[i+1])/h[i+1] - (yi[i+1] - yi[i])/h[i] );
  }

  for (i=0; i<=N-3; i++){
    a[i] = h[i+1];
  }

  for (i=0; i<=N-3; i++){
    d[i] = 2*(h[i] + h[i+1]);
  }

  /* FINDME: remove this                                                    */
  c = a;

  /* Calculate z array                                                      */
  z = tri(a, d, c, k, e, N-2);

  /* Calculate output array                                                 */
  yout = spline3(xi, yi, xout, z, h, y, nx, N);

  /* FINDME: Free arrays.                                                   */
  return yout;
}


double
splinterp_pt(double *z,
             long N,
             double *x,
             double *y,
             double xout,
             double yout){

  int index;    /* Index of x such that: x[index] <= xout < x[index+1]      */

  /* Indices of bounds:                                                     */
  int first=0,
      last =N;

  double x_hi, x_lo;
  double y_hi, y_lo;

  double h;
  double dy;

  double dx;
  double a, b, c;

  index = (first + last)/2;  /* Middle index for binary search              */
  /* Binary search to find index:                                           */
  while(first <= last){
    if(x[index] < xout && x[index + 1] > xout){
      break;
    }
    else if(x[index] < xout){
      first = index + 1;
      index = (first + last)/2;
    }
    else if(x[index] > xout){
      last  = index - 1;
      index = (first + last)/2;
    }
    else{
      break;
    }
  }

  /* x- and y-values which mark bounds of the desired value:                */
  x_lo = x[index  ];
  x_hi = x[index+1];
  y_lo = y[index  ];
  y_hi = y[index+1];

  /* Calculate range of area in question:                                   */
  /* FINDME: Rename h as dx                                                 */
  h  = x_hi - x_lo;
  dy = y_hi - y_lo;

  /* If the requested x falls on a given x, no interpolation is necessary   */
  if (x[index] == xout)
    yout = y[index];
  /* As long as the input makes sense, make calculations */
  else if(h > 0){
    dx = xout - x_lo;
    a = (z[index+1] - z[index])/(6*h);
    b = z[index]/2;
    c = dy/h - h/6 * (z[index+1] + 2*z[index]);
    yout = y_lo + dx*(c + dx*(b + dx*a));
  }
  else  /* If given a non-ascending x-array, return 0 */
    yout = 0;

  return yout;
}


double *
spline_init(double *z,
            double *x,
            double *y,
            long N){
  /* See splinterp() for description of arrays */
  double *a;
  double *d;
  double *b;
  double *h;
  double *k;
  int i;

  a = calloc(N, sizeof(double));
  d = calloc(N, sizeof(double));
  b = calloc(N, sizeof(double));
  h = calloc(N-1, sizeof(double));
  k = calloc(N-1, sizeof(double));

  /* Calculate arrays for given x, y arrays (same as in splinterp() )       */
  for(i=0;i<N-1;i++)
    h[i] = x[i+1]-x[i];

  b[0] = (y[1] - y[0]) / (x[1] - x[0]);

  for(i=0;i<N-1;i++)
    b[i+1] = -b[i] + 2*(y[i+1] - y[i]) / (x[i+1] - x[i]);

  for(i=0;i<N-2;i++)
    d[i] = 2*(h[i] + h[i+1]);

  for(i=0;i<N-2;i++)
    a[i] = h[i+1];

  for(i=0;i<N-1;i++)
    k[i] = 6*((y[i+2] - y[i+1])/h[i+1] - (y[i+1] - y[i])/h[i]);

  /* Use calculated arrays and call tri() to get z array:                   */
  z = tri(a, d, a, k, z, N-2);

  /* FINDME: Free arrays                                                    */
  return z;
}
