// Copyright (C) 2015-2016 University of Central Florida. All rights reserved.
// Transit is under an open-source, reproducible-research license (see LICENSE).

#include <math.h>
#include <numerical.h>

/* \fcnfh
   Binary search for index such that arr[index] <= val < arr[index+1]

   Return: index, on success, else:
           -1 if val    < arr[i]
           -2 if arr[f] < val
           -3 if only one element in array and it was not the looked one 
           -4 if indices are invalid:  f < i, or i < 0
           -5 if the value is the last element in arr: val = arr[f]    */
inline int
binsearchie(double *arr,  /* Array of values                           */
            long i,       /* Index of first element to consider in arr */
            long f,       /* Index of final element to consider        */
            double val){  /* Value to look for in the array            */
  long m;

  /* Exceptional cases: */
  if(arr[i]>val)
    return -1;
  if(arr[f]<val)
    return -2;
  if(arr[f]==val)
    return -5;
  if(i==f && arr[i]!=val)
    return -3;
  if(f<i && i<0)
    return -4;

  /* Binary search:     */
  while(f-i>1){
    m = (f+i)>>1;
    if(arr[m]>val)
      f = m;
    else
      i = m;
  }

  return i;
}


/* \fcnfh
   Search for index such that val is between arr[index] exclusive and
   arr[index+1] inclusive.

   @returns index that is looked for
            -1 if value is before first element 'i' of array
            -2 if value is after last element 'f' of array
            -3 if only one element in array and it was not the looked
               one 
            -4 if 'f' is equal or smaller than 'i', or 'i' is less than
               0 
            -6 if it is first index of the array
*/
inline int
binsearchei(double *arr,  /* Array of length of at least 'f'              */
            long i,       /* initial search index, cannot be negative     */
            long f,       /* final search index (or array length minus 1) */
            double val){  /* number to look for in the array              */
  long m;

  if(arr[i]>val)
    return -1;
  if(arr[i]==val)
    return -6;
  if(arr[f]<val)
    return -2;
  if(i==f&&arr[i]!=val)
    return -3;
  if(f<i&&i<0)
    return -4;

  while(f-i>1){
    m=(f+i)>>1;
    if(arr[m]<val)
      i=m;
    else
      f=m;
  }

  return i;
}

/*\fcnfh
  binsearch() defaults to an inclusive, exclusive search

  Return: output from binsearchie                          */
int 
binsearch(double *arr,  /* Array to search in              */ 
          long i,       /* Index of first element to check */
          long f,       /* Index of final element to check */
          double val){  /* Value to compare                */
  return binsearchie(arr, i, f, val);
}


/* \fcnfh
   Integrate using simpson method and trapezoid at one interval if even
   number, it requires an equispaced x-grid

   @returns value of the integration
*/
double
integ_trasim(double dx,
             double *y,
             long n)
{
  double restrap=0,res=0;
  long i;

  if(n<2){
    fprintf(stderr,
            "%s:: integ_trasim: At least 2 points are required to perform\n"
            "integration\n"
            ,__FILE__);
    exit(EXIT_FAILURE);
  }

  //use trapezoidal in the last even element
  if(!(n&1)){
    n--;
    restrap=.5*dx*(y[n]+y[n-1]);
  }

  //if there is enough elements do a Simpson integral
  if(n>2){
    //add the middle values, start with the odd elements which will be
    //multiplied by 4
    n--;
    for(i=1;i<n;i+=2)
      res+=y[i];
    res*=2;

    //now the even elements to be multiplied by 2
    for(i=2;i<n;i+=2)
      res+=y[i];
    res*=2;

    //now the borders
    res+=y[0]+y[n];
  }

  //finish'em
  return res*dx/3.0+restrap;
}


double
integ_trapz(double *x,  /* Independent variable                             */
            double *y,  /* Function to integrate                            */
            long n){    /* Number of data points                            */

  double res=0; /* Integral value                                           */
  long i;

  if(n < 2){
    fprintf(stderr, "%s:: integ_trasim: At least 2 points are required "
                    "to perform the integration\n", __FILE__);
    exit(EXIT_FAILURE);
  }
  /* Trapezoidal-rule sum:                                                  */
  for(i=0; i < n-1; i++){
    res += (x[i+1] - x[i]) * (y[i+1] + y[i]);
  }
  return 0.5*res;
}


/* \fcnfh
   Interpolates a parabola in three points and return requested
   value. X-array must need equispaced. This function doesn't check for
   less than required or non-equispaced elements.

   @returns value interpolated
*/
double
interp_parab(double *x,   /* x-array with at least 3 equispaced elements */
             double *y,   /* y-array with at least 3 elements            */
             double xr){  /* requested x-value to interpolate            */

  const double dx = x[1] - x[0];
  const double x0 = x[0] / dx;
  const double my = y[0] + y[2] - 2*y[1];
  const double a  = my / (2.0 * dx * dx);
  const double b  = (y[2] - y[1] - (x0 + 1.5) * my) / dx;
  const double c  = y[0] + x0 * ( y[2] - 4*y[1] + 3*y[0] + x0 * my ) / 2.0;

  return xr * xr * a  +  xr * b  +  c;
}


/* FUNCTION
   Interpolates a line in two points and return requested
   value. This function doesn't check for less than required elements.
   Returns: value interpolated                                              */
double
interp_line(double *x,   /* x-array with at least 2 elements                */
            double *y,   /* y-array with at least 2 elements                */
            double xr){  /* requested x-value to interpolate                */

  const double dx = x[1] - x[0];
  const double m  = (y[1] - y[0]) / dx;

  return y[0] + (xr - x[0]) * m;
}


/* FUNCTION
   return $x^n$, is a faster version of pow() that only works if n is
   integer
   Returns: result                                                          */
double 
powi(double x,
     int n){

  double y;
  _Bool negn=n<0;

  y = 1;
  if (negn)
    n *= -1;

  for(; n>0; --n){
    while((n&1)==0){
      x *= x;
      n >>= 1;
    }
    y = y*x;
  }

  if (negn)
    y = 1.0/y;

  return y;
}


/* FUNCTION
   Compares up to the 'prec' most significative digits
   Returns: true if both numbers are the same to the required accuracy      */
_Bool
fixedcmp(double d1,
         double d2,
         int prec){
  if (prec > 8){
    fprintf(stderr, "fixedcmp:: Sorry, but requested precision can be 8 at"
                    " most. Not %i. STOPPING.\n", prec);
    exit(EXIT_FAILURE);
  }

  double l10 = log(d1)/log(10);
  int expv;

  double prec10 = powi(10, prec);
  d1 *= prec10;
  d2 *= prec10;

  if (l10 < 0)
    expv = -l10 + 0.999999999999;
  else
    expv = -l10;

  prec10 = powi(10, expv);
  d1 *= prec10;
  d2 *= prec10;

  l10 = d1-d2;

  return l10*l10 < 1.0;
}


/* Resample an array by a integer factor of points                        */
int
resample(double *input,  /* Input array                                   */
         double *out,    /* Output array                                  */
         int n,          /* Number of elements in the input array         */
         int scale){     /* Resampling factor                             */
  /* Simple resample to low-res array by taking corresponding values
     from hi-res array.                                                   */
  int j;
  int m = 1 + (n-1)/scale;   /* Number of points in the resampled array   */
  for (j=0; j<m; j++){
    out[j] = input[scale*j];
  }
  return 0;
}


/* Downsample an array by a integer factor of points                        */
int
downsample(double *input,  /* Input array                                   */
           double *out,    /* Output array                                  */
           int n,          /* Number of elements in th einput array         */
           int scale){     /* Resampling factor                             */

  /* - If the scaling factor (f) is an odd number, this routine simply
       performs an averaging of the f points around the resampled value.
     - If f is even, then the boundary points are weighted by one half and
       distributed into the two resampled points.
     - The x coordinates of the first and last values of the input and output
       arrays will coincide.

     The integral area under the curves (the input and output arrays) is
     nearly conserved.

     For example, if the input array is:
       I = {I0 I1 I2 I3 I4 I5 I6}
     The output for a scaling factor f=3 is:
       O0 = (     I0 + I1) / [0.5(f+1)]
       O1 = (I2 + I3 + I4) / [    f   ]
       O2 = (I5 + I6     ) / [0.5(f+1)]
     The output for a scaling factor f=2 is:
       O0 = (         I0 + 0.5 I1) / [0.5(f+1)]
       O1 = (0.5 I1 + I2 + 0.5 I3) / [    f+1 ]
       O2 = (0.5 I3 + I4 + 0.5 I5) / [    f+1 ]
       O3 = (0.5 I5 + I6         ) / [0.5(f+1)]                             */

  int i, j;   /* for-loop indices */
  /* Number of points in the downsampled array:                             */
  int m = 1 + (n-1)/scale;
  /* Kernel size:                                                           */
  int ks = 2*(scale/2) + 1;

  /* Odd/even flag:                                                         */
  int even = 1;
  if (scale % 2 != 0)
    even = 0;

  /* First point:                                                           */
  out[0] = 0.0;
  for (i=0; i<ks/2+1; i++)
    out[0] += input[i];
  if (even == 1)
    out[0] -= 0.5*input[ks/2];
  out[0] /= 0.5*(scale+1);

  /* Down-sampling:                                                         */
  for (j=1; j<m-1; j++){
    out[j] = 0.0;
    for (i=-ks/2; i < ks/2+1; i++){
      out[j] += input[scale*j + i];
    }
    if (even == 1)
      out[j] -= 0.5*(input[scale*j-ks/2] + input[scale*j+ks/2]);
    out[j] /= scale;
  }

  /* Last point:                                                            */
  out[m-1] = 0.0;
  for (i=n-1-ks/2; i<n; i++)
    out[m-1] += input[i];
  if (even == 1)
    out[m-1] -= 0.5*input[n-ks/2];
  out[m-1] /= 0.5*(scale+1);

  return 0;
}


/* FUNCTION
   Calculate the differentials for a Simpson-rule integration.

   Parameters:
   -----------
   h: 1D double ndarray
      Intervals between the X-axis samples.

   Returns:
   --------
   hsum: 1D double ndarray
      Sums of interval pairs.
   hratio: 1D double ndarray
      Ratio of consecutive intervals.
   hfactor: 1D double ndarray
      Factor interval.

   Notes:
   ------
   - If there are even samples, skip the first interval.
   - hsum    = [h0+h1, h2+h3, h4+h5, ...]
   - hratio  = [h1/h0, h3/h2, h5/h4, ...]
   - hfactor = [hsum0*hsum0/h0*h1, hsum1*hsum1/h2*h3, ...]                  */
void
geth(double *h,
     double *hsum,
     double *hratio,
     double *hfactor,
     int n){

  int size;            /* Size of output array                              */
  int i, j, even=0;    /* Auxilliary for-loop indices                       */

  /* Empty array case:                                                      */
  if (n==0){
    hsum    = 0;
    hratio  = 0;
    hfactor = 0;
    return;
  }

  /* Check for even number of samples (odd number of intervals):            */
  even = n%2;
  /* Calculate size of h arrays                                             */
  size = (n-1)/2;

  if(even)
    even = 0;
  else
    even = 1;

  /* Calculate hsum, hratio, hfactor                                        */
  for (i=0; i<size; i++){
    j = 2*i + even;
    hsum   [i] = h[j  ] + h[j+1];
    hratio [i] = h[j+1] / h[j  ];
    hfactor[i] = hsum[i] * hsum[i] / (h[j] * h[j+1]);
  }
}


/* FUNCTION
Wrapper for Simpson-rule integration.

Parameters:
-----------
y: 1D double ndarray
   Function to integrate.
h: 1D double ndarray
   Intervals between function evaluations.
hsum: 1D double ndarray
   Sums of interval pairs.
hratio: 1D double ndarray
   Ratio of consecutive intervals.
hfactor: 1D double ndarray
   Factor interval.

Returns:
--------
integ: Float
   Integral of y over intervals h using the Simpson rule.

Notes:
------
- If there are even samples, use a trapezoidal integration for
  the first interval.
- See geth for formula for hsum, hratio, and hfactor");                     */
double
simps(double *y,
      double *h,
      double *hsum,
      double *hratio,
      double *hfactor,
      int n){
  int even;
  double integ=0;
  even = n%2 == 0;

  /* Simple case, nothing to integrate:                                     */
  if (n == 1)
    return 0.0;
  /* Simple case, do a trapezoidal integration:                             */
  if (n == 2)
    return h[0] * (y[0] + y[1]) / 2;

  /* Do Simpson integration (skip first if even):                           */
  integ = simpson(y, hsum, hratio, hfactor, n);

  /* Add trapezoidal rule for first interval if n is even:                  */
  if (even){
    integ += h[0] * (y[0] + y[1]) / 2;
  }

  return integ;
}


/* FUNCTION
   Make a spacing array for integration                                     */
void
makeh(double *x,
      double *h,
      int n){
  int i;
  /* Calculate spacing between each point:                                  */
  for(i=0; i<n-1; i++){
    h[i] = x[i+1] - x[i];
  }
}


/* FUNCTION
   Perform Simpson integration calculation                                  */
inline double
simpson(double *y,         /* Values of function to integrate               */
        double *hsum,      /* See geth function                             */
        double *hratio,    /* See geth function                             */
        double *hfactor,   /* See geth function                             */
        int n){            /* Number of elements in y                       */

  /* Do the final sum for the Simpson integration algorithm.  Based on
     the Python implementation:
     github.com/scipy/scipy/blob/v0.15.1/scipy/integrate/quadrature.py      */

  int i=0,           /* for-loop index                                      */
      j;             /* Array index for each interval                       */
  double res = 0.0;  /* The results                                         */

  /* Add contribution from each interval:                                   */
  for (i=0; i < (n-1)/2; i++){
    /* Skip first value of y if there's an even number of samples:          */
    j = 2*i + (n%2==0);
    res += (y[j  ] * (2.0 - hratio[i])     +
            y[j+1] * hfactor[i]            +
            y[j+2] * (2.0 - 1.0/hratio[i]) ) * hsum[i];
  }

  return res/6.0;
}
