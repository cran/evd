#include "header.h"

/* Conditional copulas condition on 2nd margin. */

double ccbvlog(double m1, double m2, double oldm1, double dep)
{
  double tm1,tm2,idep,u,v,fval;

  tm1 = -log(m1);
  tm2 = -log(m2);

  idep = 1/dep;
  u = R_pow(tm1, idep) + R_pow(tm2, idep);
  v = R_pow(u, dep);
  fval = exp(-v) * (1 / m2) * R_pow(tm2, idep-1) * R_pow(u, dep-1) - oldm1;
  return fval;
}

double ccbvalog(double m1, double m2, double oldm1, double dep, double asy1, 
                double asy2)
{
  double tm1,tm2,idep,u,v,fval;

  tm1 = -log(m1);
  tm2 = -log(m2);

  idep = 1/dep;
  u = R_pow(asy1*tm1, idep) + R_pow(asy2*tm2, idep);
  v = (1-asy1)*tm1 + (1-asy2)*tm2 + R_pow(u, dep);
  fval = exp(-v) * (1 / m2) * (1 - asy2 + R_pow(asy2, idep) * 
    R_pow(tm2, idep-1) * R_pow(u, dep-1)) - oldm1;
  return fval;
}

double ccbvhr(double m1, double m2, double oldm1, double dep)
{
  double tm1,tm2,v,idep,fval;

  tm1 = -log(m1);
  tm2 = -log(m2);
  idep = 1 / dep;
  v = tm2 * pnorm(idep + (log(tm2) - log(tm1)) * dep/2, 0, 1, 1, 0) +
    tm1 * pnorm(idep + (log(tm1) - log(tm2)) * dep/2, 0, 1, 1, 0);
  fval = pnorm(idep + (log(tm2) - log(tm1)) * dep/2, 0, 1, 1, 0) *
    exp(-v) / m2 - oldm1;
  return fval;
}

double ccbvneglog(double m1, double m2, double oldm1, double dep)
{
  double tm1,tm2,v,idep,fval;

  tm1 = -log(m1);
  tm2 = -log(m2);
  idep = 1 / dep;
  v = R_pow((R_pow(tm2,-dep) + R_pow(tm1,-dep)),-idep);
  fval = exp(v) * m1 * (1-R_pow(1 + R_pow(tm2/tm1,dep), -1-idep)) - oldm1;
  return fval;
}

double ccbvaneglog(double m1, double m2, double oldm1, double dep, 
                   double asy1, double asy2)
{
  double tm1,tm2,v,idep,fval;

  tm1 = -log(m1);
  tm2 = -log(m2);
  idep = 1 / dep;
  v = R_pow(asy1 * tm2, -dep) + R_pow(asy2 * tm1, -dep);
  fval = exp(R_pow(v, -idep)) * m1 * (1 - R_pow(asy1, -dep) * 
    R_pow(tm2, -dep-1) * R_pow(v, -idep-1)) - oldm1;
  return fval;
}

double ccbvbilog(double m1, double m2, double oldm1, double alpha, 
                double beta)
{
  int i;
  double tm1,tm2,v,fval;
  double delta,eps,llim,midpt,ilen,lval,midval,uval;

  tm1 = -log(m1);
  tm2 = -log(m2);

  delta = eps = R_pow(DOUBLE_EPS, 0.75);
  llim = 0;
  ilen = 1;
  lval = (1 - alpha) * tm1;
  uval = (beta - 1) * tm2;
  if(!(sign(lval) != sign(uval))) 
    error("values at end points are not of opposite sign");
  for(i=0;i<DOUBLE_DIGITS;i++) {
    ilen = ilen/2;
    midpt = llim + ilen;
    midval = (1 - alpha) * tm1 * R_pow(1 - midpt, beta) - 
             (1 - beta) * tm2 * R_pow(midpt, alpha);
    if(fabs(midval) < eps || fabs(ilen) < delta) 
      break;
    if(sign(lval) != sign(midval)) {
      uval = midval;
    }
    else {
      llim = midpt;
      lval = midval;
    }
  if(i == DOUBLE_DIGITS-1) 
    error("numerical problem in root finding algorithm");
  }

  v = tm1 * R_pow(midpt, 1 - alpha) + tm2 * R_pow(1 - midpt, 1 - beta);
  fval = exp(-v) * (1 / m2) * R_pow(1 - midpt, 1 - beta) - oldm1;
  return fval;
}

double ccbvnegbilog(double m1, double m2, double oldm1, double alpha, 
                   double beta)
{
  int i;
  double tm1,tm2,v,fval;
  double delta,eps,llim,midpt,ilen,lval,midval,uval;

  tm1 = -log(m1);
  tm2 = -log(m2);

  delta = eps = R_pow(DOUBLE_EPS, 0.75);
  llim = 0;
  ilen = 1;
  lval = - (1 + beta) * tm2;
  uval = (1 + alpha) * tm1;
  if(!(sign(lval) != sign(uval))) 
    error("values at end points are not of opposite sign1");
  for(i=0;i<DOUBLE_DIGITS;i++) {
    ilen = ilen/2;
    midpt = llim + ilen;
    midval = (1 + alpha) * tm1 * R_pow(midpt, alpha) - 
             (1 + beta) * tm2 * R_pow(1 - midpt, beta);
    if(fabs(midval) < eps || fabs(ilen) < delta) 
      break;
    if(sign(lval) != sign(midval)) {
      uval = midval;
    }
    else {
      llim = midpt;
      lval = midval;
    }
  if(i == DOUBLE_DIGITS-1) 
    error("numerical problem in root finding algorithm");
  }

  v = - tm1 - tm2 + tm1 * R_pow(midpt, 1 + alpha) + 
    tm2 * R_pow(1 - midpt, 1 + beta);
  fval = exp(v) * (1 / m2) * (1 - R_pow(1 - midpt, 1 + beta)) - oldm1;
  return fval;
}

double ccbvct(double m1, double m2, double oldm1, double alpha, double beta)
{
  double tm1,tm2,u,v,fval;

  tm1 = -log(m1);
  tm2 = -log(m2);

  u = alpha * tm2 / (alpha * tm2 + beta * tm1);
  v = tm1 * pbeta(u, alpha + 1, beta, 0, 0) + 
    tm2 * pbeta(u, alpha, beta + 1, 1, 0);
  fval = exp(-v) * (1 / m2) * pbeta(u, alpha, beta + 1, 1, 0) - oldm1;
  return fval;
}

double ccbvamix(double m1, double m2, double oldm1, double alpha, double beta)
{
  double tm1,tm2,tm1a,v,v2,fval;

  tm1 = -log(m1);
  tm2 = -log(m2);
  tm1a = tm1/(tm1 + tm2);

  v = tm1 + tm2 - tm1 * ((alpha + beta) - alpha * tm1a - beta * tm1a * tm1a);
  v2 = 1 - alpha * tm1a * tm1a - 2 * beta * tm1a * tm1a * tm1a; 
  fval = exp(-v) * (1 / m2) * v2 - oldm1;

  return fval;
}

/*
   Calculates conditional copula for any model, conditioning on
   the margin `cnd'.
*/

void ccop(double *m1, double *m2, int *cnd, double *dep, double *asy1, double *asy2, double *alpha, double *beta, int *n, int *model, double *ccop)
{
  int i;

  switch(*model) {
  case 1:
    for(i=0;i<*n;i++) {
      if(*cnd == 2) ccop[i] = ccbvlog(m1[i], m2[i], 0, *dep);
      else  ccop[i] = ccbvlog(m2[i], m1[i], 0, *dep);
    }
    break;
  case 2:
    for(i=0;i<*n;i++) {
      if(*cnd == 2) ccop[i] = ccbvalog(m1[i], m2[i], 0, *dep, *asy1, *asy2);
      else ccop[i] = ccbvalog(m2[i], m1[i], 0, *dep, *asy2, *asy1); 
    }
    break;
  case 3:
    for(i=0;i<*n;i++) {
      if(*cnd == 2) ccop[i] = ccbvhr(m1[i], m2[i], 0, *dep);
      else ccop[i] = ccbvhr(m2[i], m1[i], 0, *dep);
    }
    break;
  case 4:
    for(i=0;i<*n;i++) {
      if(*cnd == 2) ccop[i] = ccbvneglog(m1[i], m2[i], 0, *dep);
      else ccop[i] = ccbvneglog(m2[i], m1[i], 0, *dep);
    }
    break;  
  case 5:
    for(i=0;i<*n;i++) {
      if(*cnd == 2) ccop[i] = ccbvaneglog(m1[i], m2[i], 0, *dep, *asy1, *asy2);
      else ccop[i] = ccbvaneglog(m2[i], m1[i], 0, *dep, *asy2, *asy1);
    }
    break;
  case 6:
    for(i=0;i<*n;i++) {
      if(*cnd == 2) ccop[i] = ccbvbilog(m1[i], m2[i], 0, *alpha, *beta);
      else ccop[i] = ccbvbilog(m2[i], m1[i], 0, *beta, *alpha);
    }
    break;
  case 7:
    for(i=0;i<*n;i++) {
      if(*cnd == 2) ccop[i] = ccbvnegbilog(m1[i], m2[i], 0, *alpha, *beta);
      else ccop[i] = ccbvnegbilog(m2[i], m1[i], 0, *beta, *alpha);
    }
    break;
  case 8:
    for(i=0;i<*n;i++) {
      if(*cnd == 2) ccop[i] = ccbvct(m1[i], m2[i], 0, *alpha, *beta);
      else ccop[i] = ccbvct(m2[i], m1[i], 0, *beta, *alpha);
    }
    break;
    case 9:
    for(i=0;i<*n;i++) {
      if(*cnd == 2) ccop[i] = ccbvamix(m1[i], m2[i], 0, *alpha, *beta);
      else ccop[i] = ccbvamix(m2[i], m1[i], 0, *alpha + 3 * *beta, - *beta);
    }
    break;
  default:
     error("no copula found for this model");
  }
}
  



