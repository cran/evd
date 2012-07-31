#include "header.h"

void nlgpd(double *data, int *n, double *loc, double *scale, 
    double *shape, double *dns)
{
  int i;
  double *dvec, eps;
  
  dvec = (double *)R_alloc(*n, sizeof(double));
  eps = R_pow(DOUBLE_EPS, 0.3);

  if(*scale <= 0) {
     *dns = 1e6;
     return;
  }

  for(i=0;i<*n;i++)  {
    data[i] = (data[i] - *loc) / *scale;
    if(fabs(*shape) <= eps) 
      dvec[i] = log(1 / *scale) - data[i];
    else {
      data[i] = 1 + *shape * data[i];
      if(data[i] <= 0) {
        *dns = 1e6;
        return;
      }
      dvec[i] = log(1 / *scale) - (1 / *shape + 1) * log(data[i]);
    }
  }

  for(i=0;i<*n;i++) 
    *dns = *dns - dvec[i];
}

void nlpp(double *exceed, int *nhigh, double *loc, double *scale, 
    double *shape, double *thresh, double *nop, double *dns)
{
  int i;
  double *dvec, d2, u, eps;

  dvec = (double *)R_alloc(*nhigh, sizeof(double));
  eps = R_pow(DOUBLE_EPS, 0.3);

  if(*scale <= 0) {
     *dns = 1e6;
     return;
  }

  for(i=0;i<*nhigh;i++)  {
    exceed[i] = (exceed[i] - *loc) / *scale;
    if(fabs(*shape) <= eps) 
      dvec[i] = log(1 / *scale) - exceed[i];
    else {
      exceed[i] = 1 + *shape * exceed[i];
      if(exceed[i] <= 0) {
        *dns = 1e6;
        return;
      }
    dvec[i] = log(1 / *scale) - (1 / *shape + 1) * log(exceed[i]);
    }
  }

  u = (*thresh - *loc) / *scale;
  if(fabs(*shape) <= eps) 
    d2 = - *nop * exp(-u);
  else {
    u = 1 + *shape * u;
    if(u <= 0 && *shape > 0) {
      *dns = 1e6;
      return;
    }
    if(u <= 0 && *shape < 0) d2 = 0;
    else d2 = - *nop * R_pow(u, -1 / *shape);
  }

  *dns = -d2;
  for(i=0;i<*nhigh;i++) 
    *dns = *dns - dvec[i];
}

void clusters(double *high, double *high2, int *n, int *r, 
              int *rlow, double *clstrs)
{
  int i,j,rr;
  int incl = 0, clind = 0, shigh = 0, shigh2 = 0;

  for(i=0;i<*n;i++)  {
    if(high[i] && incl) {
       clstrs[i + 0 * *n] = clind;
    }
    if(high[i] && !incl) {
      incl = 1;
      clstrs[i + 1 * *n] = 1;
      clind++;
      clstrs[i + 0 * *n] = clind;
    }
    if(!high[i] && incl) {
      if(*r > *n-i) rr = *n-i;
      else rr = *r;
      for(j=i;j<(i+rr);j++) {
        shigh = shigh + high[j];
      }
      if(*rlow > *n-i) rr = *n-i;
      else rr = *rlow;
      for(j=i;j<(i+rr);j++) {
        shigh2 = shigh2 + high2[j];
      }
      if(!shigh || !shigh2) {
        incl = 0;
        clstrs[i - 1 + 2 * *n] = 1;
      }
      else clstrs[i + 0 * *n] = clind;
      shigh = shigh2 = 0;
    }
  }
  if(incl) clstrs[*n - 1 + 2 * *n] = 1;
}

