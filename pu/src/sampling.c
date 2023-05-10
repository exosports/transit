// Copyright (C) 2015-2016 University of Central Florida. All rights reserved.
// Transit is under an open-source, reproducible-research license (see LICENSE).

#include <sampling.h>

/* \fcnfh
   Finds value from dataset x,y to output yout
   Returns: 1 on success                                                    */
int 
lineinterpol(int ndat,         /* Length of x,y array                       */
             double *x,        /* x-axis                                    */
             double *y,        /* y-axis                                    */
             int n,            /* length of output arrays (indx,t,yout)     */
             long *indx,       /* Integer x-axis reference index            */
             float *t,         /* Value [0..1] to indicate how far from the
                                   value at 'indx' to the next value        */
             double *yout,     /* Interpolated output values                */
             double *dbgout){  /* Not used in linear interpolation          */
  int i, idx;

  for(i=0;i<n;i++){
    idx=indx[i];
    if(idx>=ndat||idx<0){
      fprintf(stderr,
              "%s:: index #%i(%i) is out of dataset\n"
              "range (0 - %i) in lineinterpol()\n"
              ,__FILE__,i,idx,ndat-1);
      exit(EXIT_FAILURE);
    }
    if(idx==ndat-1){
      if(t[i]!=0){
        fprintf(stderr,
                "%s:: Trying to extrapolate in last index(%i) while\n"
                "in lineinterpol(). t=%f\n"
                ,__FILE__,idx,t[i]);
        exit(EXIT_FAILURE);
      }
      yout[i]=y[idx];
    }
    else
      yout[i]=y[idx]+t[i]*(y[idx+1]-y[idx]);
  }

  return 1;
}


/*\fcnfh
  It calculates the 2nd derivatives of the function to interpolate
  storing them in 'D'. Requires the 'y' values and the x-axis spacing in
  'h' of the 'n' points dataset.
*/
inline void
natcubsplinecoef(long n,        /* Number of data points */
                 double *x,        /* y-axis values */
                 double *y,        /* y-axis values */
                 double *h,        /* x-axis spacing */
                 double *D)        /* 2nd-derivatives values */
{
  long i,n1;
  double u,v,w;

  //'n'-1 will be used a lot, so giving it to a new variable
  n1=n-1;
  //Initialization of array D as first derivative
  u=y[0];
  for(i=0;i<n1;i++){
    v=y[i+1];
    h[i]=x[i+1]-x[i];
    D[i+1]=(v-u)/h[i];
    u=v;
  }

  //Natural spline condition. D[0]=0.
  //U=h[i-1]/P[i-1], V=h[i-1], W=h[i], P[i] stored in h[i], where P[i]
  //denotes diagonal coefficient in the Gaussian elimination. h[i] is
  //only restored in \lin{hrest}
  D[0]=u=0;
  w=h[0];
  for(i=1;i<n1;i++){
    v=w;
    w=h[i];
    h[i]=(v+w)*2-u*v;
    D[i]=D[i+1]-D[i]-u*D[i-1];
    u=w/h[i];
  }

  D[n1]=0;
  //Back substitution and restoration of 'h'.
  for(i=n1-1;i>0;i--){
    w=x[i+1]-x[i];
    D[i]=(6*D[i]-w*D[i+1])/h[i];
    //\linelabel{hrest}
    h[i]=w;
  }
}


/* \fcnfh
   Interpolates in arrays (x,y), the value corresponding to 'refx'. 'n'
   is the arrays' length. Uses requested interpolation.

   @returns the interpolated y-value
 */
inline double
interp(double refx,   /* Reference x point     */
       double *x,     /* x-axis                */
       double *y,     /* y-axis                */
       long n,        /* lenght of array       */
       int intkind){  /* Kind of interpolation */

  switch(intkind){
  case INTERP_LINEAR:
    return lineinterp(refx, x, y, n);
  default:
    break;
  }
  fprintf(stderr,
          "interp(): Requested interpolation %i, is non existent.\n", intkind);
  exit(EXIT_FAILURE);
}



/* \fcnfh
   Linearly interpolates in arrays (x,y), the value corresponding to
   'refx'. 'n' is the arrays' length

   @returns the interpolated y-value
 */
double
lineinterp(double refx,                /* Reference x point */
           double *x,                /* x-axis */
           double *y,                /* y-axis */
           long n)                /* lenght of array */
{
  double yout;
  _Bool ascend=(x[1]-x[0]>0);
  n--;

  //Initial range is correct?
  if((ascend && (refx<*x || refx>x[n])) || (!ascend && (refx>*x || refx<x[n]))){
    fprintf(stderr,
            "ascend: %i, %i, %i; noascend %i, %i\n"
            "refx:%.10g  x: %.10g - %.10g\n"
            ,ascend,refx<*x,refx>x[n],refx>*x,refx<x[n]
            ,refx,x[0],x[n]);
    fprintf(stderr,
            "interpline():: refx(%.8g) is outside x range (%.8g - %.8g)\n"
            ,refx,*x,x[n]);
    exit(EXIT_FAILURE);
  }

  //check that enough data points are given
  if(n<2){
    fprintf(stderr,
            "interpline():: Array of length less than 2.\n");
    exit(EXIT_FAILURE);
  }

  //go through the arrays
  while(1){
    //Is final range correct?
    if(!n){
      if(*x==refx)
        return *y;
      fprintf(stderr,
              "interp():: Last x-value (%g) is smaller than refx(%g)\n"
              ,*x,refx);
      exit(EXIT_FAILURE);
    }
    //found!
    if((ascend && x[1]>refx) || x[1]<refx)
      break;
    n--;
    x++;
    y++;
  }

  yout=*y+(refx-*x)*(y[1]-*y)/(x[1]-*x);
  return yout;
}


//\delfh
#ifdef DBGSPLINE
#define na 5
int main(int argc, char **argv)
{
  //  double x[na]={0,1,2,3,4,5,6,7,8,9};
  //  double y[na]={0,1,2,3,4,5,6,7,8,9};
  //  double y[na]={2,1,2,3,2,1,2,3,5,9};
  double x[na]={-3,-1,0,3,4};
  double y[na]={7,11,26,56,29};
  int nn=101;
  double s,tt,hh;
  int i,ii;

  if(argc>1)
    nn=atoi(argv[1]);

  int subn=(nn-1)/(na-1)+0.5;
  nn=subn*(na-1) +1;
  long idx[nn];
  float t[nn];
  double out[nn];
  double D[na];

  for(i=0;i<nn;i++){
    idx[i]=i/subn;
    t[i]=(i%subn) / (float)subn;
  }

  natcubspline(na,x,y,nn,idx,t,out,D);
  /*
  for(i=0;i<na;i++){
    printf("%g\t%g\t%g\t%g\n",x[i],y[i],h[i],D[i]/2);
  }
  */

  printf("#x      i        t              S             S'           S\"/2          S\"'/3\n");
  for(ii=0;ii<nn-1;ii++){
    i=idx[ii];
    tt=t[ii];
    hh=x[i+1]-x[i];
    s=(y[i+1]-y[i])/hh + hh*hh/6.0*(D[i]*(-2+6*tt/hh-3*tt*tt/hh)
                                    +D[i+1]*(3*tt*tt/hh-1));

    printf("%-8.3g%-4i%6.2f%15.5g%15.5g%15.5g%15.5g\n"
           ,x[i]+tt*hh,i,tt,out[ii],
           s,(D[i]*(1-tt)+D[i+1]*tt)/2.0
           ,(D[i+1]-D[i])/3.0);
  }
  i=idx[ii];
  tt=t[ii];
  hh=x[i+1]-x[i];

  printf("%-8.3g%-4i%6.2f%15.5g\n"
         ,x[i],i,tt,out[ii]);

  
}
#endif /* DBGSPLINE */

#ifdef DBGSAMPLING
int main(int argc, char **argv){
  //  double x[na]={0,1,2,3,4,5,6,7,8,9};

  long fl=TRU_SAMPLIN;


#define ndat 8
  double x [ndat] = {1, 3, 5, 5.5, 7, 9.3, 11, 14};
  double y [ndat] = {2, 1, 0, 3.2, 3, 4.9, -1, 58};
  double y2[ndat] = {1, 2, 1, 4.2, 5, 4.3, -8, 28};
  int i;
}
#endif


//\deluh
