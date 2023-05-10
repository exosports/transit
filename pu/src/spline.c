// Copyright (C) 2015-2016 University of Central Florida. All rights reserved.
// Transit is under an open-source, reproducible-research license (see LICENSE).

#include <spline.h>
#include <iomisc.h>
#include <math.h>

/* FUNCTION
   Solve the Natural Cubic Spline Tridiagonal System.
   See: Kincaid & Cheney (2008): Numerical Mathematics & Computing,
        Algorithm 1 (page 391).                                             */
inline void
tri(double *h,  /* Input: spacing between x-coordinates for y               */
    double *y,  /* Input: Y array of samples to interpolate from            */
    double *z,  /* Output: Second derivative of the spline                  */
    long n){    /* Number of y samples                                      */

  int i;  /* Auxiliary for-loop indices                                     */
  double *b, *u, *v;
  b = calloc(n-1, sizeof(double));
  u = calloc(n-1, sizeof(double));
  v = calloc(n-1, sizeof(double));

  /* Step 1:                                                                */
  for (i=0; i<n-1; i++){
    b[i] = (y[i+1] - y[i]) / h[i];
  }

  /* Step 2:                                                                */
  u[1] = 2 * (h[1] + h[0]);
  v[1] = 6 * (b[1] - b[0]);

  for (i=2; i<n-1; i++){
    u[i] = 2*(h[i] + h[i-1]) - h[i-1]*h[i-1]/u[i-1];
    v[i] = 6*(b[i] - b[i-1]) - v[i-1]*h[i-1]/u[i-1];
  }

  /* Step 3:                                                                */
  z[0] = z[n-1] = 0;  /* First and last are zeroes                          */
  for (i=n-2; i > 0; i--){
    z[i] = (v[i] - h[i]*z[i+1]) / u[i];
  }

  free(b);
  free(u);
  free(v);
  return;
}


/* FUNCTION
   Evaluate the cubic spline  y(x) interpolating from yi(xi).
   See: Kincaid & Cheney (2008): Numerical Mathematics & Computing,
        Equation 12, (page 392).                                            */
inline void
spline3(double *xi,  /* Input X array to interpolate from                   */
        double *yi,  /* Input Y array to interpolate from                   */
        double *x,   /* Input array with X coordinates where to interpolate */
        double *z,   /* Second derivatives of yi at xi                      */
        double *h,   /* Spacing between xi points                           */
        double *y,   /* Output: spline interpolated values                  */
        long nx,     /* Length of x                                         */
        long N){     /* Length of xi                                        */

  int i, n;  /* Indices                                                     */
  double B;

  /* Calculate the spline interpolation y-value for each x-coordinate
     in the requested array                                                 */
  for (n=0; n<nx; n++){
    /* If the array is sorted:
    while (xi[i+1] < x[n]){
      i++;
    }                                                                       */
    /* Else, do a binary search:                                            */
    i = binsearchapprox(xi, x[n], 0, N-1);
    /* Enforce: xi[i] <= x[n] (except if xi[N-1] == x[n]):                  */
    if (i == N-1 || x[n] < xi[i]){
      i--;
    }

    /* Factor for linear coefficient:                                       */
    B = (yi[i+1] - yi[i]) / h[i] - h[i]/6 * (z[i+1] + 2 * z[i]);

    /* Calculate y-value from cubic polynomial:                             */
    y[n] = yi[i] +    (x[n] - xi[i])    * B         +
                   pow(x[n] - xi[i], 2) * 0.5*z[i]  +
                   pow(x[n] - xi[i], 3) * (z[i+1] - z[i]) / (6*h[i]);
  }
  return;
}


/* FUNCTION
   Cubic spline interpolation.  Given points described by arrays
   xi and yi, interpolates to coordinates of xout and puts them in yout.   */
void
splinterp(long N,         /* Length of xi                                   */
          double *xi,
          double *yi,
          long nx,        /* Length of xout                                 */
          double *xout,
          double *yout){

  double *h;  /* Spacing between xi-coordinates                             */
  double *z;  /* Array created by tri() function                            */
  int i;

  /* Allocate all arrays to be used:                                        */
  h = calloc(N-1, sizeof(double));
  z = calloc(N,   sizeof(double));

  for (i=0; i<N-1; i++){
    h[i] = xi[i+1] - xi[i];
  }

  /* Calculate z array:                                                     */
  tri(h, yi, z, N);

  /* Calculate output array                                                 */
  spline3(xi, yi, xout, z, h, yout, nx, N);

  /* Free arrays:                                                           */
  free(h);
  free(z);

  return;
}


double
splinterp_pt(double *z,
             long N,
             double *x,
             double *y,
             double xout){

  double yout;  /* Output interpolated value for xout                       */
  int index;    /* Index of x such that: x[index] <= xout < x[index+1]      */

  double x_hi, x_lo;
  double y_hi, y_lo;

  double h;
  double dy;

  double dx;
  double a, b, c;

  /* Binary search to find index:                                           */
  index = binsearchapprox(x, xout, 0, N-1);
  /* Enforce: x[i] <= xout (except if x[N-1] == xout):                      */
  if (index == N-1 || xout < x[index]){
    index--;
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
    b = 0.5*z[index];
    c = dy/h - h/6 * (z[index+1] + 2*z[index]);
    yout = y_lo + dx*(c + dx*(b + dx*a));
  }
  else  /* If given a non-ascending x-array, return 0                       */
    yout = 0;

  return yout;
}


void
spline_init(double *z,  /* Output: Second derivative of spline              */
            double *x,
            double *y,
            long N){
  /* See splinterp() for description of arrays:                             */
  double *h;
  int i;
  h = calloc(N-1, sizeof(double));

  /* Calculate arrays for given x, y arrays (same as in splinterp() )       */
  for(i=0; i<N-1; i++)
    h[i] = x[i+1] - x[i];

  /* Call tri() to get z array:                                             */
  tri(h, y, z, N);

  /* Free arrays:                                                           */
  free(h);
  return;
}
