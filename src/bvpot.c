#include "header.h"

/* Censored Likelihood Routines */

void nllbvclog(double *data1, double *data2, int *nn, int *n, double *thid, double *lambda, double *dep, double *scale1, double *shape1, double *scale2, double *shape2, double *dns)
{
  int i;
  double *dvec, *t1, *t2, *v, *v1, *v2, *v12;
  double lambda2[2], zdn;

  dvec = (double *)R_alloc(*nn, sizeof(double));
  t1 = (double *)R_alloc(*nn, sizeof(double));
  t2 = (double *)R_alloc(*nn, sizeof(double));
  v = (double *)R_alloc(*nn, sizeof(double));
  v1 = (double *)R_alloc(*nn, sizeof(double));
  v2 = (double *)R_alloc(*nn, sizeof(double));
  v12 = (double *)R_alloc(*nn, sizeof(double));

  if(*scale1 < 0.01 || *scale2 < 0.01 || *dep < 0.1 || *dep > 1) {
     *dns = 1e6;
     return;
  }

  lambda2[0] = -1/log(1 - lambda[0]);
  lambda2[1] = -1/log(1 - lambda[1]);
  lambda2[0] = R_pow(lambda2[0], -1 / *dep);
  lambda2[1] = R_pow(lambda2[1], -1 / *dep);
  zdn = R_pow(lambda2[0] + lambda2[1], *dep - 1);
  zdn = -zdn * (lambda2[0] + lambda2[1]);

  for(i=0;i<*nn;i++)  {

    data1[i] = data1[i] / *scale1;
    data2[i] = data2[i] / *scale2;
    
    if(*shape1 == 0) 
      t1[i] = exp(-data1[i]);       
    else {
      t1[i] = 1 + *shape1 * data1[i];
      if(t1[i] <= 0) {
        *dns = 1e6;
        return;
      }
      t1[i] = R_pow(t1[i], -1 / *shape1);
    }
    data1[i] = -1/log(1 - lambda[0] * t1[i]);
    
    if(*shape2 == 0) 
      t2[i] = exp(-data2[i]);       
    else {
      t2[i] = 1 + *shape2 * data2[i];
      if(t2[i] <= 0) {
        *dns = 1e6;
        return;
      }
      t2[i] = R_pow(t2[i], -1 / *shape2);
    }
    data2[i] = -1/log(1 - lambda[1] * t2[i]);

    t1[i] = R_pow(data1[i], 2) * R_pow(t1[i], 1 + *shape1) /
      (1 - lambda[0] * t1[i]);
    t1[i] = lambda[0] * t1[i] / *scale1;
    t2[i] = R_pow(data2[i], 2) * R_pow(t2[i], 1 + *shape2) /
      (1 - lambda[1] * t2[i]);
    t2[i] = lambda[1] * t2[i] / *scale2;
    
    v1[i] = R_pow(data1[i], -1 / *dep);
    v2[i] = R_pow(data2[i], -1 / *dep);
    v12[i] = R_pow(v1[i] + v2[i], *dep - 1);
    v[i] = v12[i] * (v1[i] + v2[i]);
    v1[i] = -(v1[i]/data1[i]) * v12[i];
    v2[i] = -(v2[i]/data2[i]) * v12[i];
    v12[i] = (1 - 1 / *dep) * v1[i] * v2[i] / v[i];
    
    if(thid[i] < 1.5) 
      dvec[i] = log(-v1[i]) + log(t1[i]) - v[i];
    if(thid[i] >= 1.5 && thid[i] < 2.5)
      dvec[i] = log(-v2[i]) + log(t2[i]) - v[i];
    if(thid[i] >= 2.5) dvec[i] = log(v1[i] * v2[i] - v12[i]) + 
      log(t1[i]) + log(t2[i]) - v[i];
  }

  
  for(i=0;i<*nn;i++) 
    *dns = *dns - dvec[i];
  *dns = *dns - (*n - *nn) * zdn;
}

void nllbvcbilog(double *data1, double *data2, int *nn, int *n, double *thid, double *lambda, double *alpha, double *beta, double *scale1, double *shape1, double *scale2, double *shape2, double *dns)
{
  int i,j;
  double *dvec, *t1, *t2, *v, *v1, *v2, *v12;
  double *q, *q1, *q2, *q12, *x1, *x2, *qa, *qb;
  double lambda2[2], lambda3[2], lambdaq, zdn;
  double llim,midpt,ilen,lval,midval,uval,delta,eps;

  dvec = (double *)R_alloc(*nn, sizeof(double));
  t1 = (double *)R_alloc(*nn, sizeof(double));
  t2 = (double *)R_alloc(*nn, sizeof(double));
  v = (double *)R_alloc(*nn, sizeof(double));
  v1 = (double *)R_alloc(*nn, sizeof(double));
  v2 = (double *)R_alloc(*nn, sizeof(double));
  v12 = (double *)R_alloc(*nn, sizeof(double));
  q = (double *)R_alloc(*nn, sizeof(double));
  qa = (double *)R_alloc(*nn, sizeof(double));
  qb = (double *)R_alloc(*nn, sizeof(double));
  q1 = (double *)R_alloc(*nn, sizeof(double));
  q2 = (double *)R_alloc(*nn, sizeof(double));
  q12 = (double *)R_alloc(*nn, sizeof(double));
  x1 = (double *)R_alloc(*nn, sizeof(double));
  x2 = (double *)R_alloc(*nn, sizeof(double));

  if(*scale1 < 0.01 || *scale2 < 0.01 || *alpha < 0.1 || *beta < 0.1 ||
     *alpha > 0.999 || *beta > 0.999) {
     *dns = 1e6;
     return;
  }
  delta = eps = R_pow(DOUBLE_EPS, 0.8);

  lambda2[0] = -1/log(1 - lambda[0]);
  lambda2[1] = -1/log(1 - lambda[1]);
  
  llim = 0;
  ilen = 1;
  lval = (1 - *alpha) / lambda2[0];
  uval = (*beta - 1) / lambda2[1];
  if(!(sign(lval) != sign(uval))) 
    error("values at end points are not of opposite sign");
  for(j=0;j<DOUBLE_DIGITS;j++) {
    ilen = ilen/2;
    midpt = llim + ilen;
    midval = (1 - *alpha) / lambda2[0] * R_pow(1 - midpt, *beta) - 
      (1 - *beta) / lambda2[1] * R_pow(midpt, *alpha);
    if(fabs(midval) < eps || fabs(ilen) < delta) 
      break;
    if(sign(lval) != sign(midval)) {
      uval = midval;
    }
    else {
      llim = midpt;
      lval = midval;
    }
  if(j == DOUBLE_DIGITS-1) 
    error("numerical problem in root finding algorithm");
  }
  lambdaq = midpt;

  lambda3[0] = R_pow(lambdaq, *alpha);
  lambda3[1] = R_pow(1 - lambdaq, *beta);
  zdn =  (lambdaq - 1) / (lambda3[1] * lambda2[1]) - lambdaq / 
    (lambda3[0] * lambda2[0]);

  for(i=0;i<*nn;i++)  {

      data1[i] = data1[i] / *scale1;
      data2[i] = data2[i] / *scale2;
    
      if(*shape1 == 0) 
        t1[i] = exp(-data1[i]);       
      else {
        t1[i] = 1 + *shape1 * data1[i];
        if(t1[i] <= 0) {
          *dns = 1e6;
          return;
        }
        t1[i] = R_pow(t1[i], -1 / *shape1);
      }
      data1[i] = -1/log(1 - lambda[0] * t1[i]);
    
      if(*shape2 == 0) 
        t2[i] = exp(-data2[i]);       
      else {
        t2[i] = 1 + *shape2 * data2[i];
        if(t2[i] <= 0) {
          *dns = 1e6;
          return;
        }
        t2[i] = R_pow(t2[i], -1 / *shape2);
      }
      data2[i] = -1/log(1 - lambda[1] * t2[i]);

      t1[i] = R_pow(data1[i], 2) * R_pow(t1[i], 1 + *shape1) /
        (1 - lambda[0] * t1[i]);
      t1[i] = lambda[0] * t1[i] / *scale1;
      t2[i] = R_pow(data2[i], 2) * R_pow(t2[i], 1 + *shape2) /
        (1 - lambda[1] * t2[i]);
      t2[i] = lambda[1] * t2[i] / *scale2;

      llim = 0;
	  ilen = 1;
      lval = (1 - *alpha) / data1[i];
      uval = (*beta - 1) / data2[i];
      if(!(sign(lval) != sign(uval))) 
        error("values at end points are not of opposite sign");
      for(j=0;j<DOUBLE_DIGITS;j++) {
        ilen = ilen/2;
        midpt = llim + ilen;
        midval = (1 - *alpha) / data1[i] * R_pow(1 - midpt, *beta) - 
          (1 - *beta) / data2[i] * R_pow(midpt, *alpha);
        if(fabs(midval) < eps || fabs(ilen) < delta) 
          break;
        if(sign(lval) != sign(midval)) {
          uval = midval;
        }
        else {
          llim = midpt;
          lval = midval;
        }
      if(j == DOUBLE_DIGITS-1) 
        error("numerical problem in root finding algorithm");
      }
      q[i] = midpt;

      qa[i] = R_pow(q[i], *alpha);
      qb[i] = R_pow(1 - q[i], *beta);
      x1[i] = *beta * (1 - *alpha) * qb[i] / ((1 - q[i]) * data1[i]);
      x2[i] = *alpha * (1 - *beta) * qa[i] / (q[i] * data2[i]);
      q1[i] = -(1 - *alpha) * qb[i] / (data1[i] * data1[i] * (x2[i] + x1[i]));
      q2[i] = (1 - *beta) * qa[i] / (data2[i] * data2[i] * (x2[i] + x1[i]));
      q12[i] = x2[i] * (*alpha - 1) * q2[i] / q[i] - x1[i] * (*beta - 1) * 
        q2[i] / (1 - q[i]) - x2[i] / data2[i];
      q12[i] = x1[i] * q2[i] / (data1[i] * (x2[i] + x1[i])) + 
        (1 - *alpha) * qb[i] * q12[i] / (data1[i] * data1[i] * 
	      (x2[i] + x1[i]) * (x2[i] + x1[i]));
      v[i] = q[i] / (qa[i] * data1[i]) + (1 - q[i]) / (qb[i] * data2[i]);
      v1[i] = (1 - *alpha) * q1[i] / (qa[i] * data1[i]) - (1 - *beta) * 
        q1[i] / (qb[i] * data2[i]) - q[i] / (qa[i] * data1[i] * data1[i]);
      v2[i] = (1 - *alpha) * q2[i] / (qa[i] * data1[i]) - (1 - *beta) * 
        q2[i] / (qb[i] * data2[i])-(1 - q[i])/(qb[i] * data2[i] * data2[i]);
      v12[i] = (1 - *alpha) * q12[i] / (qa[i] * data1[i]) - (1 - *alpha) * 
        q2[i] / (qa[i] * data1[i] * data1[i]) - *alpha * (1 - *alpha) * 
        q1[i] * q2[i] / (q[i] * qa[i] * data1[i]) + (1 - *beta) * q1[i] / 
        (qb[i] * data2[i] * data2[i]) - *beta * (1 - *beta) * q1[i] * 
        q2[i] / ((1 - q[i]) * qb[i] * data2[i]) - (1 - *beta) * q12[i] / 
        (qb[i] * data2[i]);

      if(thid[i] < 1.5) 
        dvec[i] = log(-v1[i]) + log(t1[i]) - v[i];
      if(thid[i] >= 1.5 && thid[i] < 2.5)
        dvec[i] = log(-v2[i]) + log(t2[i]) - v[i];
      if(thid[i] >= 2.5) dvec[i] = log(v1[i] * v2[i] - v12[i]) + 
        log(t1[i]) + log(t2[i]) - v[i];
    
  }

  for(i=0;i<*nn;i++) 
    *dns = *dns - dvec[i];
  *dns = *dns - (*n - *nn) * zdn;
}

void nllbvcalog(double *data1, double *data2, int *nn, int *n, double *thid, double *lambda, double *dep, double *asy1, double *asy2, double *scale1, double *shape1, double *scale2, double *shape2, double *dns)
{
  int i;
  double *dvec, *t1, *t2, *v, *v1, *v2, *v12;
  double *x1, *x2, *x12;
  double lambda2[2], lambda3[2], zdn;

  dvec = (double *)R_alloc(*nn, sizeof(double));
  t1 = (double *)R_alloc(*nn, sizeof(double));
  t2 = (double *)R_alloc(*nn, sizeof(double));
  v = (double *)R_alloc(*nn, sizeof(double));
  v1 = (double *)R_alloc(*nn, sizeof(double));
  v2 = (double *)R_alloc(*nn, sizeof(double));
  v12 = (double *)R_alloc(*nn, sizeof(double));
  x1 = (double *)R_alloc(*nn, sizeof(double));
  x2 = (double *)R_alloc(*nn, sizeof(double));
  x12 = (double *)R_alloc(*nn, sizeof(double));

  if(*scale1 < 0.01 || *scale2 < 0.01 || *dep < 0.1 || *dep > 1 ||
    *asy1 < 0.001 || *asy2 < 0.001 || *asy1 > 1 || *asy2 > 1) {
     *dns = 1e6;
     return;
  }

  lambda2[0] = -1/log(1 - lambda[0]);
  lambda2[1] = -1/log(1 - lambda[1]);
  lambda3[0] = R_pow(*asy1 / lambda2[0], 1 / *dep);
  lambda3[1] = R_pow(*asy2 / lambda2[1], 1 / *dep);
  zdn = R_pow(lambda3[0] + lambda3[1], *dep - 1);
  zdn = (*asy1 - 1)/lambda2[0] + (*asy2 - 1)/lambda2[1] - 
    zdn * (lambda3[0] + lambda3[1]);

  for(i=0;i<*nn;i++)  {

      data1[i] = data1[i] / *scale1;
      data2[i] = data2[i] / *scale2;
    
      if(*shape1 == 0) 
        t1[i] = exp(-data1[i]);       
      else {
        t1[i] = 1 + *shape1 * data1[i];
        if(t1[i] <= 0) {
          *dns = 1e6;
          return;
        }
        t1[i] = R_pow(t1[i], -1 / *shape1);
      }
      data1[i] = -1/log(1 - lambda[0] * t1[i]);
    
      if(*shape2 == 0) 
        t2[i] = exp(-data2[i]);       
      else {
        t2[i] = 1 + *shape2 * data2[i];
        if(t2[i] <= 0) {
          *dns = 1e6;
          return;
        }
        t2[i] = R_pow(t2[i], -1 / *shape2);
      }
      data2[i] = -1/log(1 - lambda[1] * t2[i]);

      t1[i] = R_pow(data1[i], 2) * R_pow(t1[i], 1 + *shape1) /
        (1 - lambda[0] * t1[i]);
      t1[i] = lambda[0] * t1[i] / *scale1;
      t2[i] = R_pow(data2[i], 2) * R_pow(t2[i], 1 + *shape2) /
        (1 - lambda[1] * t2[i]);
      t2[i] = lambda[1] * t2[i] / *scale2;
    
      x1[i] = R_pow(*asy1 / data1[i], 1 / *dep);
      x2[i] = R_pow(*asy2 / data2[i], 1 / *dep);
      x12[i] = R_pow(x1[i] + x2[i], *dep - 1);
      v[i] = (1 - *asy1)/data1[i] + (1 - *asy2)/data2[i] + 
        x12[i] * (x1[i] + x2[i]);
      v1[i] = ((*asy1 - 1)/data1[i] - x1[i] * x12[i]) / data1[i];
      v2[i] = ((*asy2 - 1)/data2[i] - x2[i] * x12[i]) / data2[i];
      v12[i] = (1 - 1 / *dep) * x1[i]/data1[i] * x2[i]/data2[i] * 
        x12[i] / (x1[i] + x2[i]);
    
      if(thid[i] < 1.5) 
        dvec[i] = log(-v1[i]) + log(t1[i]) - v[i];
      if(thid[i] >= 1.5 && thid[i] < 2.5)
        dvec[i] = log(-v2[i]) + log(t2[i]) - v[i];
      if(thid[i] >= 2.5) dvec[i] = log(v1[i] * v2[i] - v12[i]) + 
        log(t1[i]) + log(t2[i]) - v[i];
  }

  for(i=0;i<*nn;i++) 
    *dns = *dns - dvec[i];
  *dns = *dns - (*n - *nn) * zdn;
}

void nllbvcneglog(double *data1, double *data2, int *nn, int *n, double *thid, double *lambda, double *dep, double *scale1, double *shape1, double *scale2, double *shape2, double *dns)
{
  int i;
  double *dvec, *t1, *t2, *v, *v1, *v2, *v12;
  double *x1, *x2, *x12;
  double lambda2[2], lambda3[2], zdn;

  dvec = (double *)R_alloc(*nn, sizeof(double));
  t1 = (double *)R_alloc(*nn, sizeof(double));
  t2 = (double *)R_alloc(*nn, sizeof(double));
  v = (double *)R_alloc(*nn, sizeof(double));
  v1 = (double *)R_alloc(*nn, sizeof(double));
  v2 = (double *)R_alloc(*nn, sizeof(double));
  v12 = (double *)R_alloc(*nn, sizeof(double));
  x1 = (double *)R_alloc(*nn, sizeof(double));
  x2 = (double *)R_alloc(*nn, sizeof(double));
  x12 = (double *)R_alloc(*nn, sizeof(double));

  if(*scale1 < 0.01 || *scale2 < 0.01 || *dep < 0.05 || *dep > 5) {
     *dns = 1e6;
     return;
  }

  lambda2[0] = -1/log(1 - lambda[0]);
  lambda2[1] = -1/log(1 - lambda[1]);
  lambda3[0] = R_pow(lambda2[0], *dep);
  lambda3[1] = R_pow(lambda2[1], *dep);
  zdn = R_pow(lambda3[0] + lambda3[1], -1 / *dep - 1);
  zdn = zdn * (lambda3[0] + lambda3[1]) - 1/lambda2[0] - 1/lambda2[1];

  for(i=0;i<*nn;i++)  {

      data1[i] = data1[i] / *scale1;
      data2[i] = data2[i] / *scale2;
    
      if(*shape1 == 0) 
        t1[i] = exp(-data1[i]);       
      else {
        t1[i] = 1 + *shape1 * data1[i];
        if(t1[i] <= 0) {
          *dns = 1e6;
          return;
        }
        t1[i] = R_pow(t1[i], -1 / *shape1);
      }
      data1[i] = -1/log(1 - lambda[0] * t1[i]);
    
      if(*shape2 == 0) 
        t2[i] = exp(-data2[i]);       
      else {
        t2[i] = 1 + *shape2 * data2[i];
        if(t2[i] <= 0) {
          *dns = 1e6;
          return;
        }
        t2[i] = R_pow(t2[i], -1 / *shape2);
      }
      data2[i] = -1/log(1 - lambda[1] * t2[i]);

      t1[i] = R_pow(data1[i], 2) * R_pow(t1[i], 1 + *shape1) /
        (1 - lambda[0] * t1[i]);
      t1[i] = lambda[0] * t1[i] / *scale1;
      t2[i] = R_pow(data2[i], 2) * R_pow(t2[i], 1 + *shape2) /
        (1 - lambda[1] * t2[i]);
      t2[i] = lambda[1] * t2[i] / *scale2;
    
      x1[i] = R_pow(data1[i], *dep);
      x2[i] = R_pow(data2[i], *dep);
      x12[i] = R_pow(x1[i] + x2[i], -1 / *dep - 1);
      v[i] = 1/data1[i] + 1/data2[i] - x12[i] * (x1[i] + x2[i]);
      v1[i] = (x1[i] * x12[i] - 1/data1[i]) / data1[i];
      v2[i] = (x2[i] * x12[i] - 1/data2[i]) / data2[i];
      v12[i] = -(1 + *dep) * (x1[i]/data1[i]) * (x2[i]/data2[i]) *
        x12[i] / (x1[i] + x2[i]);
    
      if(thid[i] < 1.5) 
        dvec[i] = log(-v1[i]) + log(t1[i]) - v[i];
      if(thid[i] >= 1.5 && thid[i] < 2.5)
        dvec[i] = log(-v2[i]) + log(t2[i]) - v[i];
      if(thid[i] >= 2.5) dvec[i] = log(v1[i] * v2[i] - v12[i]) + 
        log(t1[i]) + log(t2[i]) - v[i];
  }

  for(i=0;i<*nn;i++) 
    *dns = *dns - dvec[i];
  *dns = *dns - (*n - *nn) * zdn;
}

void nllbvcnegbilog(double *data1, double *data2, int *nn, int *n, double *thid, double *lambda, double *alpha, double *beta, double *scale1, double *shape1, double *scale2, double *shape2, double *dns)
{
  int i,j;
  double *dvec, *t1, *t2, *v, *v1, *v2, *v12;
  double *q, *q1, *q2, *q12, *x1, *x2, *qa, *qb;
  double lambda2[2], lambda3[2], lambdaq, zdn;
  double llim,midpt,ilen,lval,midval,uval,delta,eps;

  dvec = (double *)R_alloc(*nn, sizeof(double));
  t1 = (double *)R_alloc(*nn, sizeof(double));
  t2 = (double *)R_alloc(*nn, sizeof(double));
  v = (double *)R_alloc(*nn, sizeof(double));
  v1 = (double *)R_alloc(*nn, sizeof(double));
  v2 = (double *)R_alloc(*nn, sizeof(double));
  v12 = (double *)R_alloc(*nn, sizeof(double));
  q = (double *)R_alloc(*nn, sizeof(double));
  qa = (double *)R_alloc(*nn, sizeof(double));
  qb = (double *)R_alloc(*nn, sizeof(double));
  q1 = (double *)R_alloc(*nn, sizeof(double));
  q2 = (double *)R_alloc(*nn, sizeof(double));
  q12 = (double *)R_alloc(*nn, sizeof(double));
  x1 = (double *)R_alloc(*nn, sizeof(double));
  x2 = (double *)R_alloc(*nn, sizeof(double));

  if(*scale1 < 0.01 || *scale2 < 0.01 || *alpha < 0.1 || *beta < 0.1 ||
     *alpha > 20 || *beta > 20) {
     *dns = 1e6;
     return;
  }
  delta = eps = R_pow(DOUBLE_EPS, 0.8);

  lambda2[0] = -1/log(1 - lambda[0]);
  lambda2[1] = -1/log(1 - lambda[1]);
  
  llim = 0;
  ilen = 1;
  uval = (1 + *alpha) / lambda2[0];
  lval = - (1 + *beta) / lambda2[1];
  if(!(sign(lval) != sign(uval))) 
    error("values at end points are not of opposite sign");
  for(j=0;j<DOUBLE_DIGITS;j++) {
    ilen = ilen/2;
    midpt = llim + ilen;
    midval = (1 + *alpha) / lambda2[0] * R_pow(midpt, *alpha) - 
          (1 + *beta) / lambda2[1] * R_pow(1 - midpt, *beta);
    if(fabs(midval) < eps || fabs(ilen) < delta) 
      break;
    if(sign(lval) != sign(midval)) {
      uval = midval;
    }
    else {
      llim = midpt;
      lval = midval;
    }
  if(j == DOUBLE_DIGITS-1) 
    error("numerical problem in root finding algorithm");
  }
  lambdaq = midpt;

  lambda3[0] = R_pow(lambdaq, *alpha);
  lambda3[1] = R_pow(1 - lambdaq, *beta);
  zdn = (lambdaq * lambda3[0] - 1) / lambda2[0] + ((1 - lambdaq) * 
    lambda3[1] - 1) / lambda2[1];

  for(i=0;i<*nn;i++)  {

      data1[i] = data1[i] / *scale1;
      data2[i] = data2[i] / *scale2;
    
      if(*shape1 == 0) 
        t1[i] = exp(-data1[i]);       
      else {
        t1[i] = 1 + *shape1 * data1[i];
        if(t1[i] <= 0) {
          *dns = 1e6;
          return;
        }
        t1[i] = R_pow(t1[i], -1 / *shape1);
      }
      data1[i] = -1/log(1 - lambda[0] * t1[i]);
    
      if(*shape2 == 0) 
        t2[i] = exp(-data2[i]);       
      else {
        t2[i] = 1 + *shape2 * data2[i];
        if(t2[i] <= 0) {
          *dns = 1e6;
          return;
        }
        t2[i] = R_pow(t2[i], -1 / *shape2);
      }
      data2[i] = -1/log(1 - lambda[1] * t2[i]);

      t1[i] = R_pow(data1[i], 2) * R_pow(t1[i], 1 + *shape1) /
        (1 - lambda[0] * t1[i]);
      t1[i] = lambda[0] * t1[i] / *scale1;
      t2[i] = R_pow(data2[i], 2) * R_pow(t2[i], 1 + *shape2) /
        (1 - lambda[1] * t2[i]);
      t2[i] = lambda[1] * t2[i] / *scale2;
      
      llim = 0;
	  ilen = 1;
      uval = (1 + *alpha) / data1[i];
      lval = - (1 + *beta) / data2[i];
      if(!(sign(lval) != sign(uval))) 
        error("values at end points are not of opposite sign");
      for(j=0;j<DOUBLE_DIGITS;j++) {
        ilen = ilen/2;
        midpt = llim + ilen;
        midval = (1 + *alpha) / data1[i] * R_pow(midpt, *alpha) - 
          (1 + *beta) / data2[i] * R_pow(1 - midpt, *beta);
        if(fabs(midval) < eps || fabs(ilen) < delta) 
          break;
        if(sign(lval) != sign(midval)) {
          uval = midval;
        }
        else {
          llim = midpt;
          lval = midval;
        }
      if(j == DOUBLE_DIGITS-1) 
        error("numerical problem in root finding algorithm");
      }
      q[i] = midpt;

      qa[i] = R_pow(q[i], *alpha);
      qb[i] = R_pow(1 - q[i], *beta);
      x1[i] = *alpha * (1 + *alpha) * qa[i] / (q[i] * data1[i]);
      x2[i] = *beta * (1 + *beta) * qb[i] / ((1 - q[i]) * data2[i]);
      q1[i] = (1 + *alpha) * qa[i] / (data1[i] * data1[i] * (x2[i] + x1[i]));
      q2[i] = -(1 + *beta) * qb[i] / (data2[i] * data2[i] * (x2[i] + x1[i]));
      q12[i] = x1[i] * (*alpha - 1) * q2[i] / q[i] - x2[i] * (*beta - 1) * 
        q2[i] / (1 - q[i]) - x2[i] / data2[i];
      q12[i] = x1[i] * q2[i] / (data1[i] * (x2[i] + x1[i])) - 
        (1 + *alpha) * qa[i] * q12[i] / (data1[i] * data1[i] * 
	      (x2[i] + x1[i]) * (x2[i] + x1[i]));
      v[i] = (1 - q[i] * qa[i]) / data1[i] + (1 - (1 - q[i]) * qb[i]) / 
        data2[i];
      v1[i] = (q[i] * qa[i] - 1) / (data1[i] * data1[i]);
      v2[i] = (qb[i] * (1 - q[i]) - 1) / (data2[i] * data2[i]);
      v12[i] = (1 + *alpha) * qa[i] * q2[i] / (data1[i] * data1[i]);

      if(thid[i] < 1.5) 
        dvec[i] = log(-v1[i]) + log(t1[i]) - v[i];
      if(thid[i] >= 1.5 && thid[i] < 2.5)
        dvec[i] = log(-v2[i]) + log(t2[i]) - v[i];
      if(thid[i] >= 2.5) dvec[i] = log(v1[i] * v2[i] - v12[i]) + 
        log(t1[i]) + log(t2[i]) - v[i];
  }

  for(i=0;i<*nn;i++)
    *dns = *dns - dvec[i];
  *dns = *dns - (*n - *nn) * zdn;
}

void nllbvcaneglog(double *data1, double *data2, int *nn, int *n, double *thid, double *lambda, double *dep, double *asy1, double *asy2, double *scale1, double *shape1, double *scale2, double *shape2, double *dns)
{
  int i;
  double *dvec, *t1, *t2, *v, *v1, *v2, *v12;
  double *x1, *x2, *x12;
  double lambda2[2], lambda3[2], zdn;

  dvec = (double *)R_alloc(*nn, sizeof(double));
  t1 = (double *)R_alloc(*nn, sizeof(double));
  t2 = (double *)R_alloc(*nn, sizeof(double));
  v = (double *)R_alloc(*nn, sizeof(double));
  v1 = (double *)R_alloc(*nn, sizeof(double));
  v2 = (double *)R_alloc(*nn, sizeof(double));
  v12 = (double *)R_alloc(*nn, sizeof(double));
  x1 = (double *)R_alloc(*nn, sizeof(double));
  x2 = (double *)R_alloc(*nn, sizeof(double));
  x12 = (double *)R_alloc(*nn, sizeof(double));

  if(*scale1 < 0.01 || *scale2 < 0.01 || *dep < 0.05 || *dep > 5 ||
     *asy1 < 0.001 || *asy2 < 0.001 || *asy1 > 1 || *asy2 > 1) {
     *dns = 1e6;
     return;
  }

  lambda2[0] = -1/log(1 - lambda[0]);
  lambda2[1] = -1/log(1 - lambda[1]);
  lambda3[0] = R_pow(lambda2[0] / *asy1, *dep);
  lambda3[1] = R_pow(lambda2[1] / *asy2, *dep);
  zdn = R_pow(lambda3[0] + lambda3[1], -1 / *dep - 1);
  zdn = zdn * (lambda3[0] + lambda3[1]) - 1/lambda2[0] - 1/lambda2[1];

  for(i=0;i<*nn;i++)  {

      data1[i] = data1[i] / *scale1;
      data2[i] = data2[i] / *scale2;
    
      if(*shape1 == 0) 
        t1[i] = exp(-data1[i]);       
      else {
        t1[i] = 1 + *shape1 * data1[i];
        if(t1[i] <= 0) {
          *dns = 1e6;
          return;
        }
        t1[i] = R_pow(t1[i], -1 / *shape1);
      }
      data1[i] = -1/log(1 - lambda[0] * t1[i]);
    
      if(*shape2 == 0) 
        t2[i] = exp(-data2[i]);       
      else {
        t2[i] = 1 + *shape2 * data2[i];
        if(t2[i] <= 0) {
          *dns = 1e6;
          return;
        }
        t2[i] = R_pow(t2[i], -1 / *shape2);
      }
      data2[i] = -1/log(1 - lambda[1] * t2[i]);

      t1[i] = R_pow(data1[i], 2) * R_pow(t1[i], 1 + *shape1) /
        (1 - lambda[0] * t1[i]);
      t1[i] = lambda[0] * t1[i] / *scale1;
      t2[i] = R_pow(data2[i], 2) * R_pow(t2[i], 1 + *shape2) /
        (1 - lambda[1] * t2[i]);
      t2[i] = lambda[1] * t2[i] / *scale2;
    
      x1[i] = R_pow(data1[i] / *asy1, *dep);
      x2[i] = R_pow(data2[i] / *asy2, *dep);
      x12[i] = R_pow(x1[i] + x2[i], -1 / *dep - 1);
      v[i] = 1/data1[i] + 1/data2[i] - x12[i] * (x1[i] + x2[i]);
      v1[i] = (x1[i] * x12[i] - 1/data1[i]) / data1[i];
      v2[i] = (x2[i] * x12[i] - 1/data2[i]) / data2[i];
      v12[i] = -(1 + *dep) * x1[i]/data1[i] * x2[i]/data2[i] * 
        x12[i] / (x1[i] + x2[i]);
    
      if(thid[i] < 1.5) 
        dvec[i] = log(-v1[i]) + log(t1[i]) - v[i];
      if(thid[i] >= 1.5 && thid[i] < 2.5)
        dvec[i] = log(-v2[i]) + log(t2[i]) - v[i];
      if(thid[i] >= 2.5) dvec[i] = log(v1[i] * v2[i] - v12[i]) + 
        log(t1[i]) + log(t2[i]) - v[i];
  }

  for(i=0;i<*nn;i++) 
    *dns = *dns - dvec[i];
  *dns = *dns - (*n - *nn) * zdn;
}

void nllbvcct(double *data1, double *data2, int *nn, int *n, double *thid, double *lambda, double *alpha, double *beta, double *scale1, double *shape1, double *scale2, double *shape2, double *dns)
{
  int i;
  double *dvec, *t1, *t2, *v, *v1, *v2, *v12;
  double *x;
  double lambda2[2], zdn;

  dvec = (double *)R_alloc(*nn, sizeof(double));
  t1 = (double *)R_alloc(*nn, sizeof(double));
  t2 = (double *)R_alloc(*nn, sizeof(double));
  v = (double *)R_alloc(*nn, sizeof(double));
  v1 = (double *)R_alloc(*nn, sizeof(double));
  v2 = (double *)R_alloc(*nn, sizeof(double));
  v12 = (double *)R_alloc(*nn, sizeof(double));
  x = (double *)R_alloc(*nn, sizeof(double));

  if(*scale1 < 0.01 || *scale2 < 0.01 || *alpha < 0.001 || *beta < 0.001
    || *alpha > 30 || *beta > 30) {
     *dns = 1e6;
     return;
  }

  lambda2[0] = -1/log(1 - lambda[0]);
  lambda2[1] = -1/log(1 - lambda[1]);
  zdn = *alpha * lambda2[0] / (*alpha * lambda2[0] + *beta * lambda2[1]);
  zdn = -pbeta(zdn, *alpha + 1, *beta, 0, 0) / lambda2[0] -
    pbeta(zdn, *alpha, *beta + 1, 1, 0) / lambda2[1]; 

  for(i=0;i<*nn;i++)  {

      data1[i] = data1[i] / *scale1;
      data2[i] = data2[i] / *scale2;
    
      if(*shape1 == 0) 
        t1[i] = exp(-data1[i]);       
      else {
        t1[i] = 1 + *shape1 * data1[i];
        if(t1[i] <= 0) {
          *dns = 1e6;
          return;
        }
        t1[i] = R_pow(t1[i], -1 / *shape1);
      }
      data1[i] = -1/log(1 - lambda[0] * t1[i]);
    
      if(*shape2 == 0) 
        t2[i] = exp(-data2[i]);       
      else {
        t2[i] = 1 + *shape2 * data2[i];
        if(t2[i] <= 0) {
          *dns = 1e6;
          return;
        }
        t2[i] = R_pow(t2[i], -1 / *shape2);
      }
      data2[i] = -1/log(1 - lambda[1] * t2[i]);

      t1[i] = R_pow(data1[i], 2) * R_pow(t1[i], 1 + *shape1) /
        (1 - lambda[0] * t1[i]);
      t1[i] = lambda[0] * t1[i] / *scale1;
      t2[i] = R_pow(data2[i], 2) * R_pow(t2[i], 1 + *shape2) /
        (1 - lambda[1] * t2[i]);
      t2[i] = lambda[1] * t2[i] / *scale2;
    
      x[i] = *alpha * data1[i] / (*alpha * data1[i] + *beta * data2[i]);
      v[i] = pbeta(x[i], *alpha + 1, *beta, 0, 0) / data1[i] +
        pbeta(x[i], *alpha, *beta + 1, 1, 0) / data2[i];  
      v1[i] = -pbeta(x[i], *alpha + 1, *beta, 0, 0) / R_pow(data1[i], 2);
      v2[i] = -pbeta(x[i], *alpha, *beta + 1, 1, 0) / R_pow(data2[i], 2);
      v12[i] = -(*alpha * *beta) * dbeta(x[i], *alpha + 1, *beta, 0) / 
        (data1[i] * R_pow(*alpha * data1[i] + *beta * data2[i], 2));
    
      if(thid[i] < 1.5) 
        dvec[i] = log(-v1[i]) + log(t1[i]) - v[i];
      if(thid[i] >= 1.5 && thid[i] < 2.5)
        dvec[i] = log(-v2[i]) + log(t2[i]) - v[i];
      if(thid[i] >= 2.5) dvec[i] = log(v1[i] * v2[i] - v12[i]) + 
        log(t1[i]) + log(t2[i]) - v[i];
  }

  for(i=0;i<*nn;i++) 
    *dns = *dns - dvec[i];
  *dns = *dns - (*n - *nn) * zdn;
}

void nllbvchr(double *data1, double *data2, int *nn, int *n, double *thid, double *lambda, double *dep, double *scale1, double *shape1, double *scale2, double *shape2, double *dns)
{
  int i;
  double *dvec, *t1, *t2, *v, *v1, *v2, *v12;
  double lambda2[2], zdn, idep;

  dvec = (double *)R_alloc(*nn, sizeof(double));
  t1 = (double *)R_alloc(*nn, sizeof(double));
  t2 = (double *)R_alloc(*nn, sizeof(double));
  v = (double *)R_alloc(*nn, sizeof(double));
  v1 = (double *)R_alloc(*nn, sizeof(double));
  v2 = (double *)R_alloc(*nn, sizeof(double));
  v12 = (double *)R_alloc(*nn, sizeof(double));

  if(*scale1 < 0.01 || *scale2 < 0.01 || *dep < 0.2 || *dep > 10) {
     *dns = 1e6;
     return;
  }

  idep = 1/ *dep;
  lambda2[0] = -1/log(1 - lambda[0]);
  lambda2[1] = -1/log(1 - lambda[1]);
  zdn = -1/lambda2[0] * pnorm(idep + *dep * (log(lambda2[1]) - 
    log(lambda2[0]))/2, 0, 1, 1, 0) - 1/lambda2[1] * pnorm(idep + 
    *dep * (log(lambda2[0]) - log(lambda2[1]))/2, 0, 1, 1, 0);

  for(i=0;i<*nn;i++)  {

    data1[i] = data1[i] / *scale1;
    data2[i] = data2[i] / *scale2;
    
    if(*shape1 == 0) 
      t1[i] = exp(-data1[i]);       
    else {
      t1[i] = 1 + *shape1 * data1[i];
      if(t1[i] <= 0) {
        *dns = 1e6;
        return;
      }
      t1[i] = R_pow(t1[i], -1 / *shape1);
    }
    data1[i] = -1/log(1 - lambda[0] * t1[i]);
    
    if(*shape2 == 0) 
      t2[i] = exp(-data2[i]);       
    else {
      t2[i] = 1 + *shape2 * data2[i];
      if(t2[i] <= 0) {
        *dns = 1e6;
        return;
      }
      t2[i] = R_pow(t2[i], -1 / *shape2);
    }
    data2[i] = -1/log(1 - lambda[1] * t2[i]);

    t1[i] = R_pow(data1[i], 2) * R_pow(t1[i], 1 + *shape1) /
      (1 - lambda[0] * t1[i]);
    t1[i] = lambda[0] * t1[i] / *scale1;
    t2[i] = R_pow(data2[i], 2) * R_pow(t2[i], 1 + *shape2) /
      (1 - lambda[1] * t2[i]);
    t2[i] = lambda[1] * t2[i] / *scale2;
    
    idep = 1/ *dep;
    v[i] = 1/data1[i] * pnorm(idep + *dep * (log(data2[i]) - 
      log(data1[i]))/2, 0, 1, 1, 0) + 1/data2[i] * pnorm(idep + 
      *dep * (log(data1[i]) - log(data2[i]))/2, 0, 1, 1, 0);
    v1[i] = -1/R_pow(data1[i], 2) *
      pnorm(idep + *dep * (log(data2[i]) - log(data1[i]))/2, 0, 1, 1, 0);
    v2[i] = -1/R_pow(data2[i], 2) *
      pnorm(idep + *dep * (log(data1[i]) - log(data2[i]))/2, 0, 1, 1, 0);
    v12[i] = - *dep / (2 * data1[i] * data2[i]) / data1[i] *
      dnorm(idep + *dep * (log(data2[i]) - log(data1[i]))/2, 0, 1, 0);
    
    if(thid[i] < 1.5) 
      dvec[i] = log(-v1[i]) + log(t1[i]) - v[i];
    if(thid[i] >= 1.5 && thid[i] < 2.5)
      dvec[i] = log(-v2[i]) + log(t2[i]) - v[i];
    if(thid[i] >= 2.5) dvec[i] = log(v1[i] * v2[i] - v12[i]) + 
      log(t1[i]) + log(t2[i]) - v[i];
  }

  for(i=0;i<*nn;i++) 
    *dns = *dns - dvec[i];
  *dns = *dns - (*n - *nn) * zdn;
}

void nllbvcamix(double *data1, double *data2, int *nn, int *n, double *thid, double *lambda, double *alpha, double *beta, double *scale1, double *shape1, double *scale2, double *shape2, double *dns)
{
  int i;
  double *dvec, *t1, *t2, *v, *v1, *v2, *v12;
  double *x;
  double lambda2[2], zdn;

  dvec = (double *)R_alloc(*nn, sizeof(double));
  t1 = (double *)R_alloc(*nn, sizeof(double));
  t2 = (double *)R_alloc(*nn, sizeof(double));
  v = (double *)R_alloc(*nn, sizeof(double));
  v1 = (double *)R_alloc(*nn, sizeof(double));
  v2 = (double *)R_alloc(*nn, sizeof(double));
  v12 = (double *)R_alloc(*nn, sizeof(double));
  x = (double *)R_alloc(*nn, sizeof(double));

  if(*scale1 < 0.01 || *scale2 < 0.01 || 
     *alpha < 0 || *alpha + 3 * *beta < 0 || 
     *alpha + *beta > 1 || *alpha + 2 * *beta > 1) {
     *dns = 1e6;
     return;
  }

  lambda2[0] = -1/log(1 - lambda[0]);
  lambda2[1] = -1/log(1 - lambda[1]);
  lambda2[0] = 1/lambda2[0];
  lambda2[1] = 1/lambda2[1];
  zdn = lambda2[0]/(lambda2[0] + lambda2[1]);
  zdn = -lambda2[0] - lambda2[1] + (*alpha + *beta) * lambda2[0] -
    *alpha * lambda2[0] * zdn - *beta * lambda2[0] * zdn * zdn; 

  for(i=0;i<*nn;i++)  {

      data1[i] = data1[i] / *scale1;
      data2[i] = data2[i] / *scale2;
    
      if(*shape1 == 0) 
        t1[i] = exp(-data1[i]);       
      else {
        t1[i] = 1 + *shape1 * data1[i];
        if(t1[i] <= 0) {
          *dns = 1e6;
          return;
        }
        t1[i] = R_pow(t1[i], -1 / *shape1);
      }
      data1[i] = -1/log(1 - lambda[0] * t1[i]);
    
      if(*shape2 == 0) 
        t2[i] = exp(-data2[i]);       
      else {
        t2[i] = 1 + *shape2 * data2[i];
        if(t2[i] <= 0) {
          *dns = 1e6;
          return;
        }
        t2[i] = R_pow(t2[i], -1 / *shape2);
      }
      data2[i] = -1/log(1 - lambda[1] * t2[i]);

      t1[i] = R_pow(data1[i], 2) * R_pow(t1[i], 1 + *shape1) /
        (1 - lambda[0] * t1[i]);
      t1[i] = lambda[0] * t1[i] / *scale1;
      t2[i] = R_pow(data2[i], 2) * R_pow(t2[i], 1 + *shape2) /
        (1 - lambda[1] * t2[i]);
      t2[i] = lambda[1] * t2[i] / *scale2;
      
      x[i] = 1 / (data1[i] + data2[i]);
      v[i] = 1/data1[i] + 1/data2[i] - (*alpha + *beta) / data1[i] +
        *alpha * data2[i] * x[i] / data1[i] + 
        *beta * data2[i] * data2[i] * x[i] * x[i] / data1[i];  
      v1[i] = -1 / (data1[i] * data1[i]) + *alpha * x[i] * x[i] +
        *beta * x[i] * x[i] * x[i] * (data1[i] + 3 * data2[i]);
      v2[i] = -1 / (data2[i] * data2[i]) + *alpha * x[i] * x[i] +
        2 * *beta * x[i] * x[i] * x[i] * data2[i];
      v12[i] = -2 * *alpha * x[i] * x[i] * x[i] - 
        6 * *beta * x[i] * x[i] * x[i] * x[i] * data2[i];
      
      if(thid[i] < 1.5) 
        dvec[i] = log(-v1[i]) + log(t1[i]) - v[i];
      if(thid[i] >= 1.5 && thid[i] < 2.5)
        dvec[i] = log(-v2[i]) + log(t2[i]) - v[i];
      if(thid[i] >= 2.5) dvec[i] = log(v1[i] * v2[i] - v12[i]) + 
        log(t1[i]) + log(t2[i]) - v[i];
  }

  for(i=0;i<*nn;i++) 
    *dns = *dns - dvec[i];
  *dns = *dns - (*n - *nn) * zdn;
}

/* Point Process Likelihood Routines */

void nllbvplog(double *data1, double *data2, int *nn, double *thid, double *r1, double *r2, double *p, double *dep, double *scale1, double *shape1, double *scale2, double *shape2, double *dns)
{
  int i;
  double *dvec, *r, *w, *jac, *h;
  double idep, v, utt[2];

  dvec = (double *)R_alloc(*nn, sizeof(double));
  r = (double *)R_alloc(*nn, sizeof(double));
  w = (double *)R_alloc(*nn, sizeof(double));
  jac = (double *)R_alloc(*nn, sizeof(double));
  h = (double *)R_alloc(*nn, sizeof(double));

  if(*scale1 < 0.01 || *scale2 < 0.01 || *dep < 0.1 || *dep > 1) {
     *dns = 1e6;
     return;
  }

  for(i=0;i<*nn;i++)  {

    data1[i] = data1[i] / *scale1;
    data2[i] = data2[i] / *scale2;
    
    if(*shape1 == 0) 
      data1[i] = exp(-data1[i]);       
    else {
      data1[i] = 1 + *shape1 * data1[i];
      if(data1[i] <= 0) {
        *dns = 1e6;
        return;
      }
      data1[i] = R_pow(data1[i], -1 / *shape1);
    }
    data1[i] = -1/log(1 - r1[i] * data1[i]);
    
    if(*shape2 == 0) 
      data2[i] = exp(-data2[i]);       
    else {
      data2[i] = 1 + *shape2 * data2[i];
      if(data2[i] <= 0) {
        *dns = 1e6;
        return;
      }
      data2[i] = R_pow(data2[i], -1 / *shape2);
    }
    data2[i] = -1/log(1 - r2[i] * data2[i]);
		
		r[i] = log(data1[i] + data2[i]);
    w[i] = data1[i] / exp(r[i]);

    if(thid[i] < 1.5) 
      jac[i] = 2 * log(data1[i]) + 1 / data1[i] + (1 + *shape1) * log(1 - 
        exp(-1 / data1[i])) - log(*scale1) - *shape1 * log(p[0]);
    if(thid[i] >= 1.5 && thid[i] < 2.5)
      jac[i] = 2 * log(data2[i]) + 1 / data2[i] + (1 + *shape2) * log(1 - 
        exp(-1 / data2[i])) - log(*scale2) - *shape2 * log(p[1]);
    if(thid[i] >= 2.5)
      jac[i] = 2 * log(data1[i]) + 1 / data1[i] + (1 + *shape1) * log(1 - 
        exp(-1 / data1[i])) - log(*scale1) - *shape1 * log(p[0]) +
        2 * log(data2[i]) + 1 / data2[i] + (1 + *shape2) * log(1 - 
        exp(-1 / data2[i])) - log(*scale2) - *shape2 * log(p[1]);

    idep = 1 / *dep;
    h[i] = log(idep - 1) - (1+idep) * log(w[i] * (1-w[i])) + 
      (*dep - 2) * log(R_pow(w[i],-idep) + R_pow(1-w[i],-idep));

    dvec[i] = jac[i] + h[i] - 3 * r[i];
  }
  
  for(i=0;i<*nn;i++)
    *dns = *dns - dvec[i];  

  utt[0] = -1 / log(1 - p[0]);
  utt[1] = -1 / log(1 - p[1]);
  v = R_pow(R_pow(utt[0],-1 / *dep) + R_pow(utt[1],-1 / *dep), *dep);
  
	*dns = *dns + v;
}

void nllbvpneglog(double *data1, double *data2, int *nn, double *thid, double *r1, double *r2, double *p, double *dep, double *scale1, double *shape1, double *scale2, double *shape2, double *dns)
{
  int i;
  double *dvec, *r, *w, *jac, *h;
  double v, utt[2];

  dvec = (double *)R_alloc(*nn, sizeof(double));
  r = (double *)R_alloc(*nn, sizeof(double));
  w = (double *)R_alloc(*nn, sizeof(double));
  jac = (double *)R_alloc(*nn, sizeof(double));
  h = (double *)R_alloc(*nn, sizeof(double));

  if(*scale1 < 0.01 || *scale2 < 0.01 || *dep < 0.05 || *dep > 5) {
     *dns = 1e6;
     return;
  }

  for(i=0;i<*nn;i++)  {

    data1[i] = data1[i] / *scale1;
    data2[i] = data2[i] / *scale2;
    
    if(*shape1 == 0) 
      data1[i] = exp(-data1[i]);       
    else {
      data1[i] = 1 + *shape1 * data1[i];
      if(data1[i] <= 0) {
        *dns = 1e6;
        return;
      }
      data1[i] = R_pow(data1[i], -1 / *shape1);
    }
    data1[i] = -1/log(1 - r1[i] * data1[i]);
    
    if(*shape2 == 0) 
      data2[i] = exp(-data2[i]);       
    else {
      data2[i] = 1 + *shape2 * data2[i];
      if(data2[i] <= 0) {
        *dns = 1e6;
        return;
      }
      data2[i] = R_pow(data2[i], -1 / *shape2);
    }
    data2[i] = -1/log(1 - r2[i] * data2[i]);

		r[i] = log(data1[i] + data2[i]);
    w[i] = data1[i] / exp(r[i]);

    if(thid[i] < 1.5) 
      jac[i] = 2 * log(data1[i]) + 1 / data1[i] + (1 + *shape1) * log(1 - 
        exp(-1 / data1[i])) - log(*scale1) - *shape1 * log(p[0]);
    if(thid[i] >= 1.5 && thid[i] < 2.5)
      jac[i] = 2 * log(data2[i]) + 1 / data2[i] + (1 + *shape2) * log(1 - 
        exp(-1 / data2[i])) - log(*scale2) - *shape2 * log(p[1]);
    if(thid[i] >= 2.5)
      jac[i] = 2 * log(data1[i]) + 1 / data1[i] + (1 + *shape1) * log(1 - 
        exp(-1 / data1[i])) - log(*scale1) - *shape1 * log(p[0]) +
        2 * log(data2[i]) + 1 / data2[i] + (1 + *shape2) * log(1 - 
        exp(-1 / data2[i])) - log(*scale2) - *shape2 * log(p[1]);

    h[i] = 1 / (1 + R_pow(1 / w[i] - 1, *dep));
    h[i] = log(*dep + 1) + log(1 - h[i]) + (1 + 1 / *dep) * log(h[i]) - 
      log(1 - w[i]) - 2 * log(w[i]);

    dvec[i] = jac[i] + h[i] - 3 * r[i];
  }
  
  for(i=0;i<*nn;i++)
    *dns = *dns - dvec[i];  

  utt[0] = -1 / log(1 - p[0]);
  utt[1] = -1 / log(1 - p[1]);
  v = 1 / utt[0] + 1 / utt[1] - R_pow(R_pow(utt[0], *dep) + 
    R_pow(utt[1], *dep), -1 / *dep);
  
  *dns = *dns + v;
}

void nllbvpct(double *data1, double *data2, int *nn, double *thid, double *r1, double *r2, double *p, double *alpha, double *beta, double *scale1, double *shape1, double *scale2, double *shape2, double *dns)
{
  int i;
  double *dvec, *r, *w, *jac, *h;
  double v, utt[2];

  dvec = (double *)R_alloc(*nn, sizeof(double));
  r = (double *)R_alloc(*nn, sizeof(double));
  w = (double *)R_alloc(*nn, sizeof(double));
  jac = (double *)R_alloc(*nn, sizeof(double));
  h = (double *)R_alloc(*nn, sizeof(double));

  if(*scale1 < 0.01 || *scale2 < 0.01 || *alpha < 0.001 || *beta < 0.001
    || *alpha > 30 || *beta > 30) {
     *dns = 1e6;
     return;
  }

  for(i=0;i<*nn;i++)  {

    data1[i] = data1[i] / *scale1;
    data2[i] = data2[i] / *scale2;
    
    if(*shape1 == 0) 
      data1[i] = exp(-data1[i]);       
    else {
      data1[i] = 1 + *shape1 * data1[i];
      if(data1[i] <= 0) {
        *dns = 1e6;
        return;
      }
      data1[i] = R_pow(data1[i], -1 / *shape1);
    }
    data1[i] = -1/log(1 - r1[i] * data1[i]);
    
    if(*shape2 == 0) 
      data2[i] = exp(-data2[i]);       
    else {
      data2[i] = 1 + *shape2 * data2[i];
      if(data2[i] <= 0) {
        *dns = 1e6;
        return;
      }
      data2[i] = R_pow(data2[i], -1 / *shape2);
    }
    data2[i] = -1/log(1 - r2[i] * data2[i]);

		r[i] = log(data1[i] + data2[i]);
    w[i] = data1[i] / exp(r[i]);

    if(thid[i] < 1.5) 
      jac[i] = 2 * log(data1[i]) + 1 / data1[i] + (1 + *shape1) * log(1 - 
        exp(-1 / data1[i])) - log(*scale1) - *shape1 * log(p[0]);
    if(thid[i] >= 1.5 && thid[i] < 2.5)
      jac[i] = 2 * log(data2[i]) + 1 / data2[i] + (1 + *shape2) * log(1 - 
        exp(-1 / data2[i])) - log(*scale2) - *shape2 * log(p[1]);
    if(thid[i] >= 2.5)
      jac[i] = 2 * log(data1[i]) + 1 / data1[i] + (1 + *shape1) * log(1 - 
        exp(-1 / data1[i])) - log(*scale1) - *shape1 * log(p[0]) +
        2 * log(data2[i]) + 1 / data2[i] + (1 + *shape2) * log(1 - 
        exp(-1 / data2[i])) - log(*scale2) - *shape2 * log(p[1]);

    h[i] = (*alpha + *beta + 1) * log(*alpha * w[i] + *beta * (1-w[i]))
      + lgammafn(*alpha) + lgammafn(*beta);
    h[i] = lgammafn(*alpha + *beta + 1) + *alpha * log(*alpha) + *beta * 
      log(*beta) + (*alpha - 1) * log(w[i]) + (*beta - 1) * log(1-w[i])
      - h[i];

    dvec[i] = jac[i] + h[i] - 3 * r[i];
  }
  
  for(i=0;i<*nn;i++)
    *dns = *dns - dvec[i];  

  utt[0] = -1 / log(1 - p[0]);
  utt[1] = -1 / log(1 - p[1]);
  v = *alpha * utt[0] /(*alpha * utt[0] + *beta * utt[1]);
  v = pbeta(v, *alpha + 1, *beta, 0, 0) / utt[0] + pbeta(v, *alpha, 
    *beta + 1, 1, 0) / utt[1];
  
  *dns = *dns + v;
}

void nllbvpbilog(double *data1, double *data2, int *nn, double *thid, double *r1, double *r2, double *p, double *alpha, double *beta, double *scale1, double *shape1, double *scale2, double *shape2, double *dns)
{
  int i,j;
  double *dvec, *r, *w, *jac, *h;
  double v, utt[2];
  double llim,midpt,ilen,lval,midval,uval,delta,eps;

  dvec = (double *)R_alloc(*nn, sizeof(double));
  r = (double *)R_alloc(*nn, sizeof(double));
  w = (double *)R_alloc(*nn, sizeof(double));
  jac = (double *)R_alloc(*nn, sizeof(double));
  h = (double *)R_alloc(*nn, sizeof(double));

  if(*scale1 < 0.01 || *scale2 < 0.01 || *alpha < 0.1 || *beta < 0.1
    || *alpha > 0.999 || *beta > 0.999) {
     *dns = 1e6;
     return;
  }
  delta = eps = R_pow(DOUBLE_EPS, 0.8);

  for(i=0;i<*nn;i++)  {

    data1[i] = data1[i] / *scale1;
    data2[i] = data2[i] / *scale2;
    
    if(*shape1 == 0) 
      data1[i] = exp(-data1[i]);       
    else {
      data1[i] = 1 + *shape1 * data1[i];
      if(data1[i] <= 0) {
        *dns = 1e6;
        return;
      }
      data1[i] = R_pow(data1[i], -1 / *shape1);
    }
    data1[i] = -1/log(1 - r1[i] * data1[i]);
    
    if(*shape2 == 0) 
      data2[i] = exp(-data2[i]);       
    else {
      data2[i] = 1 + *shape2 * data2[i];
      if(data2[i] <= 0) {
        *dns = 1e6;
        return;
      }
      data2[i] = R_pow(data2[i], -1 / *shape2);
    }
    data2[i] = -1/log(1 - r2[i] * data2[i]);

		r[i] = log(data1[i] + data2[i]);
    w[i] = data1[i] / exp(r[i]);

    if(thid[i] < 1.5) 
      jac[i] = 2 * log(data1[i]) + 1 / data1[i] + (1 + *shape1) * log(1 - 
        exp(-1 / data1[i])) - log(*scale1) - *shape1 * log(p[0]);
    if(thid[i] >= 1.5 && thid[i] < 2.5)
      jac[i] = 2 * log(data2[i]) + 1 / data2[i] + (1 + *shape2) * log(1 - 
        exp(-1 / data2[i])) - log(*scale2) - *shape2 * log(p[1]);
    if(thid[i] >= 2.5)
      jac[i] = 2 * log(data1[i]) + 1 / data1[i] + (1 + *shape1) * log(1 - 
        exp(-1 / data1[i])) - log(*scale1) - *shape1 * log(p[0]) +
        2 * log(data2[i]) + 1 / data2[i] + (1 + *shape2) * log(1 - 
        exp(-1 / data2[i])) - log(*scale2) - *shape2 * log(p[1]);
     
    llim = 0;
	  ilen = 1;
    lval = (1 - *alpha) * (1 - w[i]);
    uval = (*beta - 1) * w[i];
    if(!(sign(lval) != sign(uval))) 
      error("values at end points are not of opposite sign");
    for(j=0;j<DOUBLE_DIGITS;j++) {
      ilen = ilen/2;
      midpt = llim + ilen;
      midval = (1 - *alpha) * (1-w[i]) * R_pow(1 - midpt, *beta) - 
        (1 - *beta) * w[i] * R_pow(midpt, *alpha);
      if(fabs(midval) < eps || fabs(ilen) < delta) 
        break;
      if(sign(lval) != sign(midval)) {
        uval = midval;
      }
      else {
        llim = midpt;
        lval = midval;
      }
    if(j == DOUBLE_DIGITS-1) 
      error("numerical problem in root finding algorithm");
    }

    h[i] = log(1 - *alpha) + log(1 - midpt) + (1 - *alpha) * log(midpt) - 
      2 * log(w[i]) - log(1-w[i]) - log(*alpha * (1-midpt) + *beta * midpt);

    dvec[i] = jac[i] + h[i] - 3 * r[i];
  }
  
  for(i=0;i<*nn;i++)
    *dns = *dns - dvec[i];  

  utt[0] = -1 / log(1 - p[0]);
  utt[1] = -1 / log(1 - p[1]);
  llim = 0;
  ilen = 1;
  lval = (1 - *alpha) / utt[0];
  uval = (*beta - 1) / utt[1];
  if(!(sign(lval) != sign(uval))) 
    error("values at end points are not of opposite sign");
  for(j=0;j<DOUBLE_DIGITS;j++) {
    ilen = ilen/2;
    midpt = llim + ilen;
    midval = (1 - *alpha) / utt[0] * R_pow(1 - midpt, *beta) - 
      (1 - *beta) / utt[1] * R_pow(midpt, *alpha);
    if(fabs(midval) < eps || fabs(ilen) < delta) 
      break;
    if(sign(lval) != sign(midval)) {
      uval = midval;
    }
    else {
      llim = midpt;
      lval = midval;
    }
  if(j == DOUBLE_DIGITS-1) 
    error("numerical problem in root finding algorithm");
  }

  v = R_pow(midpt, 1 - *alpha) / utt[0] + R_pow(1-midpt, 1 - *beta) / 
    utt[1];
  
  *dns = *dns + v;
}

void nllbvpnegbilog(double *data1, double *data2, int *nn, double *thid, double *r1, double *r2, double *p, double *alpha, double *beta, double *scale1, double *shape1, double *scale2, double *shape2, double *dns)
{
  int i,j;
  double *dvec, *r, *w, *jac, *h;
  double v, utt[2];
  double llim,midpt,ilen,lval,midval,uval,delta,eps;

  dvec = (double *)R_alloc(*nn, sizeof(double));
  r = (double *)R_alloc(*nn, sizeof(double));
  w = (double *)R_alloc(*nn, sizeof(double));
  jac = (double *)R_alloc(*nn, sizeof(double));
  h = (double *)R_alloc(*nn, sizeof(double));

  if(*scale1 < 0.01 || *scale2 < 0.01 || *alpha < 0.1 || *beta < 0.1
    || *alpha > 20 || *beta > 20) {
     *dns = 1e6;
     return;
  }
  delta = eps = R_pow(DOUBLE_EPS, 0.8);

  for(i=0;i<*nn;i++)  {

    data1[i] = data1[i] / *scale1;
    data2[i] = data2[i] / *scale2;
    
    if(*shape1 == 0) 
      data1[i] = exp(-data1[i]);       
    else {
      data1[i] = 1 + *shape1 * data1[i];
      if(data1[i] <= 0) {
        *dns = 1e6;
        return;
      }
      data1[i] = R_pow(data1[i], -1 / *shape1);
    }
    data1[i] = -1/log(1 - r1[i] * data1[i]);
    
    if(*shape2 == 0) 
      data2[i] = exp(-data2[i]);       
    else {
      data2[i] = 1 + *shape2 * data2[i];
      if(data2[i] <= 0) {
        *dns = 1e6;
        return;
      }
      data2[i] = R_pow(data2[i], -1 / *shape2);
    }
    data2[i] = -1/log(1 - r2[i] * data2[i]);

		r[i] = log(data1[i] + data2[i]);
    w[i] = data1[i] / exp(r[i]);

    if(thid[i] < 1.5) 
      jac[i] = 2 * log(data1[i]) + 1 / data1[i] + (1 + *shape1) * log(1 - 
        exp(-1 / data1[i])) - log(*scale1) - *shape1 * log(p[0]);
    if(thid[i] >= 1.5 && thid[i] < 2.5)
      jac[i] = 2 * log(data2[i]) + 1 / data2[i] + (1 + *shape2) * log(1 - 
        exp(-1 / data2[i])) - log(*scale2) - *shape2 * log(p[1]);
    if(thid[i] >= 2.5)
      jac[i] = 2 * log(data1[i]) + 1 / data1[i] + (1 + *shape1) * log(1 - 
        exp(-1 / data1[i])) - log(*scale1) - *shape1 * log(p[0]) +
        2 * log(data2[i]) + 1 / data2[i] + (1 + *shape2) * log(1 - 
        exp(-1 / data2[i])) - log(*scale2) - *shape2 * log(p[1]);
     
    llim = 0;
	  ilen = 1;
    uval = (1 + *alpha) * (1 - w[i]);
    lval = - (1 + *beta) * w[i];
    if(!(sign(lval) != sign(uval))) 
      error("values at end points are not of opposite sign");
    for(j=0;j<DOUBLE_DIGITS;j++) {
      ilen = ilen/2;
      midpt = llim + ilen;
      midval = (1 + *alpha) * (1-w[i]) * R_pow(midpt, *alpha) - 
        (1 + *beta) * w[i] * R_pow(1-midpt, *beta);
      if(fabs(midval) < eps || fabs(ilen) < delta) 
        break;
      if(sign(lval) != sign(midval)) {
        uval = midval;
      }
      else {
        llim = midpt;
        lval = midval;
      }
    if(j == DOUBLE_DIGITS-1) 
      error("numerical problem in root finding algorithm");
    }

    h[i] = log(1 + *alpha) + log(1 - midpt) + (1 + *alpha) * log(midpt) - 
      2 * log(w[i]) - log(1-w[i]) - log(*alpha * (1-midpt) + *beta * midpt);

    dvec[i] = jac[i] + h[i] - 3 * r[i];
  }
  
  for(i=0;i<*nn;i++)
    *dns = *dns - dvec[i];  

  utt[0] = -1 / log(1 - p[0]);
  utt[1] = -1 / log(1 - p[1]);
  llim = 0;
  ilen = 1;
  uval = (1 + *alpha) / utt[0];
  lval = - (1 + *beta) / utt[1];
  if(!(sign(lval) != sign(uval))) 
    error("values at end points are not of opposite sign");
  for(j=0;j<DOUBLE_DIGITS;j++) {
    ilen = ilen/2;
    midpt = llim + ilen;
    midval = (1 + *alpha) / utt[0] * R_pow(midpt, *alpha) - 
        (1 + *beta) / utt[1] * R_pow(1-midpt, *beta);
    if(fabs(midval) < eps || fabs(ilen) < delta) 
      break;
    if(sign(lval) != sign(midval)) {
      uval = midval;
    }
    else {
      llim = midpt;
      lval = midval;
    }
  if(j == DOUBLE_DIGITS-1) 
    error("numerical problem in root finding algorithm");
  }

  v = (1 - R_pow(midpt, 1 + *alpha)) / utt[0] + (1 - R_pow(1-midpt, 
    1 + *beta)) / utt[1];
  
  *dns = *dns + v;
}

void nllbvphr(double *data1, double *data2, int *nn, double *thid, double *r1, double *r2, double *p, double *dep, double *scale1, double *shape1, double *scale2, double *shape2, double *dns)
{
  int i;
  double *dvec, *r, *w, *jac, *h;
  double idep, v, utt[2];

  dvec = (double *)R_alloc(*nn, sizeof(double));
  r = (double *)R_alloc(*nn, sizeof(double));
  w = (double *)R_alloc(*nn, sizeof(double));
  jac = (double *)R_alloc(*nn, sizeof(double));
  h = (double *)R_alloc(*nn, sizeof(double));

  if(*scale1 < 0.01 || *scale2 < 0.01 || *dep < 0.2 || *dep > 10) {
     *dns = 1e6;
     return;
  } 

  for(i=0;i<*nn;i++)  {

    data1[i] = data1[i] / *scale1;
    data2[i] = data2[i] / *scale2;
    
    if(*shape1 == 0) 
      data1[i] = exp(-data1[i]);       
    else {
      data1[i] = 1 + *shape1 * data1[i];
      if(data1[i] <= 0) {
        *dns = 1e6;
        return;
      }
      data1[i] = R_pow(data1[i], -1 / *shape1);
    }
    data1[i] = -1/log(1 - r1[i] * data1[i]);
    
    if(*shape2 == 0) 
      data2[i] = exp(-data2[i]);       
    else {
      data2[i] = 1 + *shape2 * data2[i];
      if(data2[i] <= 0) {
        *dns = 1e6;
        return;
      }
      data2[i] = R_pow(data2[i], -1 / *shape2);
    }
    data2[i] = -1/log(1 - r2[i] * data2[i]);

		r[i] = log(data1[i] + data2[i]);
    w[i] = data1[i] / exp(r[i]);

    if(thid[i] < 1.5) 
      jac[i] = 2 * log(data1[i]) + 1 / data1[i] + (1 + *shape1) * log(1 - 
        exp(-1 / data1[i])) - log(*scale1) - *shape1 * log(p[0]);
    if(thid[i] >= 1.5 && thid[i] < 2.5)
      jac[i] = 2 * log(data2[i]) + 1 / data2[i] + (1 + *shape2) * log(1 - 
        exp(-1 / data2[i])) - log(*scale2) - *shape2 * log(p[1]);
    if(thid[i] >= 2.5)
      jac[i] = 2 * log(data1[i]) + 1 / data1[i] + (1 + *shape1) * log(1 - 
        exp(-1 / data1[i])) - log(*scale1) - *shape1 * log(p[0]) +
        2 * log(data2[i]) + 1 / data2[i] + (1 + *shape2) * log(1 - 
        exp(-1 / data2[i])) - log(*scale2) - *shape2 * log(p[1]);

    idep = 1 / *dep;
    h[i] = log(*dep / 2) - 2 * log(w[i]) - log(1-w[i]) +
      dnorm(idep + *dep * (log(1-w[i]) - log(w[i]))/2, 0, 1, 1);

    dvec[i] = jac[i] + h[i] - 3 * r[i];
  }
  
  for(i=0;i<*nn;i++)
    *dns = *dns - dvec[i];  

  utt[0] = -1 / log(1 - p[0]);
  utt[1] = -1 / log(1 - p[1]);
  v = pnorm(1 / *dep + *dep * log(utt[1]/utt[0]) / 2, 0, 1, 1, 0) / utt[0] + 
    pnorm(1 / *dep + *dep * log(utt[0]/utt[1]) / 2, 0, 1, 1, 0) / utt[1];
  
  *dns = *dns + v;
}


