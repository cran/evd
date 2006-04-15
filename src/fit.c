#include "header.h" 

void nlgev(double *data, int *n, double *loc, double *scale, double *shape, 
           double *dns)
{
  int i;
  double *dvec;

  dvec = (double *)R_alloc(*n, sizeof(double));

  for(i=0;i<*n;i++)  {
    data[i] = (data[i] - loc[i]) / *scale;
    if(*shape == 0) 
      dvec[i] = log(1 / *scale) - data[i] - exp(-data[i]);
    else {
      data[i] = 1 + *shape * data[i];
      if(data[i] <= 0) {
        *dns = 1e6;
        return;
      }
      dvec[i] = log(1 / *scale) - R_pow(data[i], -1 / *shape) -
                (1 / *shape + 1) * log(data[i]);
    }
  }

  for(i=0;i<*n;i++) 
    *dns = *dns - dvec[i];
}

void nlbvlog(double *datam1, double *datam2, int *n, int *si, double *dep, 
             double *loc1, double *scale1, double *shape1, double *loc2, 
             double *scale2, double *shape2, int *split, double *dns)
{
  int i;
  double idep, *dvec, *z;

  dvec = (double *)R_alloc(*n, sizeof(double));
  z = (double *)R_alloc(*n, sizeof(double));
  
  for(i=0;i<*n;i++) {
    datam1[i] = (datam1[i] - loc1[i]) / *scale1;
    datam2[i] = (datam2[i] - loc2[i]) / *scale2;
    if(*shape1 == 0) 
        datam1[i] = -datam1[i];
    else {
      datam1[i] = 1 + *shape1 * datam1[i];
      if(datam1[i] <= 0) {
        *dns = 1e6;
        return;
      }
      else datam1[i] = -1 / *shape1 * log(datam1[i]);
    }
    if(*shape2 == 0) 
      datam2[i] = -datam2[i];
    else {
      datam2[i] = 1 + *shape2 * datam2[i];
      if(datam2[i] <= 0) {
        *dns = 1e6;
        return;
      }
      else datam2[i] = -1 / *shape2 * log(datam2[i]);
    }
  }

  idep = 1 / *dep;
  for(i=0;i<*n;i++) {
    z[i] = R_pow(exp(idep * datam1[i]) + exp(idep * datam2[i]), *dep);
    dvec[i] = (idep + *shape1) * datam1[i] + (idep + *shape2) *
      datam2[i] - log(*scale1 * *scale2);
    dvec[i] = dvec[i] + (1-2*idep)*log(z[i]) - z[i];
    if(si[i] == 0) dvec[i] = dvec[i] + log(z[i]);
    else if(si[i] == 1) dvec[i] = dvec[i] + log(idep-1);
    else dvec[i] = dvec[i] + log(idep-1+z[i]);
  }
  
  if(*split > 0.5) {
    for(i=0;i<*n;i++) dns[i] = - dvec[i];
  }
  else {
    for(i=0;i<*n;i++) *dns = *dns - dvec[i];
  }
}

void nlbvalog(double *datam1, double *datam2, int *n, int *si, double *dep,
	      double *asy1, double *asy2, double *loc1, double *scale1, 
              double *shape1, double *loc2, double *scale2, double *shape2, 
              int *split, double *dns)
{
  int i;
  double idep,c1,c2,c3,c4;
  double *e1,*e2,*e3,*e4,*z,*v,*jc,*dvec;

  e1 = (double *)R_alloc(*n, sizeof(double));
  e2 = (double *)R_alloc(*n, sizeof(double));
  e3 = (double *)R_alloc(*n, sizeof(double));
  e4 = (double *)R_alloc(*n, sizeof(double));
  z = (double *)R_alloc(*n, sizeof(double));
  v = (double *)R_alloc(*n, sizeof(double));
  jc = (double *)R_alloc(*n, sizeof(double));
  dvec = (double *)R_alloc(*n, sizeof(double));
  
  for(i=0;i<*n;i++) {
    datam1[i] = (datam1[i] - loc1[i]) / *scale1;
    datam2[i] = (datam2[i] - loc2[i]) / *scale2;
    if(*shape1 == 0) 
        datam1[i] = -datam1[i];
    else {
      datam1[i] = 1 + *shape1 * datam1[i];
      if(datam1[i] <= 0) {
        *dns = 1e6;
        return;
      }
      else datam1[i] = -1 / *shape1 * log(datam1[i]);
    }
    if(*shape2 == 0) 
      datam2[i] = -datam2[i];
    else {
      datam2[i] = 1 + *shape2 * datam2[i];
      if(datam2[i] <= 0) {
        *dns = 1e6;
        return;
      }
      else datam2[i] = -1 / *shape2 * log(datam2[i]);
    }
  }

  idep = 1 / *dep;
  c1 = log(1 - *asy1) + log(1 - *asy2);
  c2 = idep * log(*asy1) + idep * log(*asy2);
  c3 = log(1 - *asy1) + idep * log(*asy2);
  c4 = log(1 - *asy2) + idep * log(*asy1);

  for(i=0;i<*n;i++) {
    z[i] = R_pow(exp(idep * (log(*asy1) + datam1[i])) + 
      exp(idep * (log(*asy2) + datam2[i])), *dep);
    v[i] = (1 - *asy1) * exp(datam1[i]) + (1 - *asy2) * exp(datam2[i]) + z[i];
    jc[i] = (1 + *shape1) * datam1[i] + (1 + *shape2) * datam2[i] -
      log(*scale1 * *scale2);
    e1[i] = c3 + (idep - 1) * datam2[i];
    e2[i] = c4 + (idep - 1) * datam1[i];
    e3[i] = (1 - idep) * log(z[i]) + log(exp(e1[i]) + exp(e2[i]));
    e4[i] = c2 + (idep - 1) * datam1[i] + (idep - 1) * datam2[i] +
      (1 - 2*idep) * log(z[i]);
    
    dvec[i] =  jc[i] - v[i];
    if(si[i] == 0) {
      e4[i] = e4[i] + log(z[i]);
      dvec[i] = dvec[i] + log(exp(c1) + exp(e3[i]) + exp(e4[i]));
    }
    else if(si[i] == 1) {
      e4[i] = e4[i] + log(idep-1);
      dvec[i] = dvec[i] + e4[i];
    }
    else {
      e4[i] = e4[i] + log(idep-1+z[i]);
      dvec[i] = dvec[i] + log(exp(c1) + exp(e3[i]) + exp(e4[i]));
    }
  }
  
  if(*split > 0.5) {
    for(i=0;i<*n;i++) dns[i] = - dvec[i];
  }
  else {
    for(i=0;i<*n;i++) *dns = *dns - dvec[i];
  }
}

void nlbvhr(double *datam1, double *datam2, int *n, int *si, double *dep, 
            double *loc1, double *scale1, double *shape1, double *loc2, 
            double *scale2, double *shape2, int *split, double *dns)
{
  int i;
  double idep;
  double *e1,*e2,*e3,*v,*jc,*dvec;

  e1 = (double *)R_alloc(*n, sizeof(double));
  e2 = (double *)R_alloc(*n, sizeof(double));
  e3 = (double *)R_alloc(*n, sizeof(double));
  v = (double *)R_alloc(*n, sizeof(double));
  jc = (double *)R_alloc(*n, sizeof(double));
  dvec = (double *)R_alloc(*n, sizeof(double));
  
  for(i=0;i<*n;i++) {
    datam1[i] = (datam1[i] - loc1[i]) / *scale1;
    datam2[i] = (datam2[i] - loc2[i]) / *scale2;
    if(*shape1 == 0) 
        datam1[i] = -datam1[i];
    else {
      datam1[i] = 1 + *shape1 * datam1[i];
      if(datam1[i] <= 0) {
        *dns = 1e6;
        return;
      }
      else datam1[i] = -1 / *shape1 * log(datam1[i]);
    }
    if(*shape2 == 0) 
      datam2[i] = -datam2[i];
    else {
      datam2[i] = 1 + *shape2 * datam2[i];
      if(datam2[i] <= 0) {
        *dns = 1e6;
        return;
      }
      else datam2[i] = -1 / *shape2 * log(datam2[i]);
    }
  }

  idep = 1 / *dep;
  for(i=0;i<*n;i++) {
    e1[i] = exp(datam1[i]) * 
      pnorm(idep + *dep * (datam1[i] - datam2[i]) / 2, 0, 1, 1, 0);
    e2[i] = exp(datam2[i]) * 
      pnorm(idep + *dep * (datam2[i] - datam1[i]) / 2, 0, 1, 1, 0);
    e3[i] = exp(datam1[i]) * 
      dnorm(idep + *dep * (datam1[i] - datam2[i]) / 2, 0, 1, 0);
    v[i] = e1[i] + e2[i];
    if(si[i] == 0) dvec[i] = e1[i] * e2[i];
    else if(si[i] == 1) dvec[i] = *dep * e3[i] / 2;
    else dvec[i] = e1[i] * e2[i] + *dep * e3[i] / 2;
    jc[i] = *shape1 * datam1[i] + *shape2 * datam2[i] -
      log(*scale1 * *scale2);
    dvec[i] = log(dvec[i]) + jc[i] - v[i];
  }

  if(*split > 0.5) {
    for(i=0;i<*n;i++) dns[i] = - dvec[i];
  }
  else {
    for(i=0;i<*n;i++) *dns = *dns - dvec[i];
  }
}

void nlbvneglog(double *datam1, double *datam2, int *n, int *si, double *dep, 
                double *loc1, double *scale1, double *shape1, double *loc2, 
                double *scale2, double *shape2, int *split, double *dns)
{
  int i;
  double idep;
  double *e1,*e2,*z,*v,*jc,*dvec;

  e1 = (double *)R_alloc(*n, sizeof(double));
  e2 = (double *)R_alloc(*n, sizeof(double));
  z = (double *)R_alloc(*n, sizeof(double));
  v = (double *)R_alloc(*n, sizeof(double));
  jc = (double *)R_alloc(*n, sizeof(double));
  dvec = (double *)R_alloc(*n, sizeof(double));
  
  for(i=0;i<*n;i++) {
    datam1[i] = (datam1[i] - loc1[i]) / *scale1;
    datam2[i] = (datam2[i] - loc2[i]) / *scale2;
    if(*shape1 == 0) 
        datam1[i] = -datam1[i];
    else {
      datam1[i] = 1 + *shape1 * datam1[i];
      if(datam1[i] <= 0) {
        *dns = 1e6;
        return;
      }
      else datam1[i] = -1 / *shape1 * log(datam1[i]);
    }
    if(*shape2 == 0) 
      datam2[i] = -datam2[i];
    else {
      datam2[i] = 1 + *shape2 * datam2[i];
      if(datam2[i] <= 0) {
        *dns = 1e6;
        return;
      }
      else datam2[i] = -1 / *shape2 * log(datam2[i]);
    }
  }

  idep = 1 / *dep;
  for(i=0;i<*n;i++) {
    z[i] = R_pow(exp(-*dep * datam1[i]) + exp(-*dep * datam2[i]), -idep);
    v[i] = exp(datam1[i]) + exp(datam2[i]) - z[i];
    jc[i] = (1 + *shape1)*datam1[i] + (1 + *shape2)*datam2[i] -
      log(*scale1 * *scale2);
    e1[i] = (1 + *dep) * log(z[i]) +
      log(exp((-*dep-1) * datam1[i]) + exp((-*dep-1) * datam2[i]));
    e2[i] = (-*dep-1) * datam1[i] + (-*dep-1) * datam2[i] +
      (1 + 2 * *dep) * log(z[i]);

    dvec[i] =  jc[i] - v[i];
    if(si[i] == 0) {
      e2[i] = e2[i] + log(z[i]);
      dvec[i] = dvec[i] + log(1 - exp(e1[i]) + exp(e2[i]));
    }
    else if(si[i] == 1) {
      e2[i] = e2[i] + log(1 + *dep);
      dvec[i] = dvec[i] + e2[i];
    }
    else {
      e2[i] = e2[i] + log(1 + *dep + z[i]);
      dvec[i] = dvec[i] + log(1 - exp(e1[i]) + exp(e2[i]));
    }
  }
 
  if(*split > 0.5) {
    for(i=0;i<*n;i++) dns[i] = - dvec[i];
  }
  else {
    for(i=0;i<*n;i++) *dns = *dns - dvec[i];
  }
}

void nlbvaneglog(double *datam1, double *datam2, int *n, int *si, double *dep,
	         double *asy1, double *asy2, double *loc1, double *scale1, 
                 double *shape1, double *loc2, double *scale2, double *shape2, 
                 int *split, double *dns)
{
  int i;
  double idep;
  double *e1,*e2,*e3,*e4,*z,*v,*jc,*dvec;

  e1 = (double *)R_alloc(*n, sizeof(double));
  e2 = (double *)R_alloc(*n, sizeof(double));
  e3 = (double *)R_alloc(*n, sizeof(double));
  e4 = (double *)R_alloc(*n, sizeof(double));
  z = (double *)R_alloc(*n, sizeof(double));
  v = (double *)R_alloc(*n, sizeof(double));
  jc = (double *)R_alloc(*n, sizeof(double));
  dvec = (double *)R_alloc(*n, sizeof(double));
  
  for(i=0;i<*n;i++) {
    datam1[i] = (datam1[i] - loc1[i]) / *scale1;
    datam2[i] = (datam2[i] - loc2[i]) / *scale2;
    if(*shape1 == 0) 
        datam1[i] = -datam1[i];
    else {
      datam1[i] = 1 + *shape1 * datam1[i];
      if(datam1[i] <= 0) {
        *dns = 1e6;
        return;
      }
      else datam1[i] = -1 / *shape1 * log(datam1[i]);
    }
    if(*shape2 == 0) 
      datam2[i] = -datam2[i];
    else {
      datam2[i] = 1 + *shape2 * datam2[i];
      if(datam2[i] <= 0) {
        *dns = 1e6;
        return;
      }
      else datam2[i] = -1 / *shape2 * log(datam2[i]);
    }
  }
  idep = 1 / *dep;
  for(i=0;i<*n;i++) {
    z[i] = R_pow(exp(-*dep * (log(*asy1) + datam1[i])) + 
      exp(-*dep * (log(*asy2) + datam2[i])), -idep);
    v[i] = exp(datam1[i]) + exp(datam2[i]) - z[i];
    jc[i] = (1 + *shape1)*datam1[i] + (1 + *shape2)*datam2[i] -
      log(*scale1 * *scale2);
    e1[i] = -*dep * log(*asy1) + (-*dep - 1) * datam1[i];
    e2[i] = -*dep * log(*asy2) + (-*dep - 1) * datam2[i];
    e3[i] = (1 + *dep) * log(z[i]) + log(exp(e1[i]) + exp(e2[i]));
    e4[i] = e1[i] + e2[i] + (1 + 2 * *dep) * log(z[i]);

    dvec[i] =  jc[i] - v[i];
    if(si[i] == 0) {
      e4[i] = e4[i] + log(z[i]);
      dvec[i] = dvec[i] + log(1 - exp(e3[i]) + exp(e4[i]));
    }
    else if(si[i] == 1) {
      e4[i] = e4[i] + log(1 + *dep);
      dvec[i] = dvec[i] + e4[i];
    }
    else {
      e4[i] = e4[i] + log(1 + *dep + z[i]);
      dvec[i] = dvec[i] + log(1 - exp(e3[i]) + exp(e4[i]));
    }
  }

  if(*split > 0.5) {
    for(i=0;i<*n;i++) dns[i] = - dvec[i];
  }
  else {
    for(i=0;i<*n;i++) *dns = *dns - dvec[i];
  }
}

void nlbvbilog(double *datam1, double *datam2, int *n, int *si, double *alpha,
	      double *beta, double *loc1, double *scale1, double *shape1, 
              double *loc2, double *scale2, double *shape2, int *split, 
              double *dns)
{
  int i,j;
  double *e1,*e2,*v,*jc,*dvec,*gma;
  double llim,midpt,ulim,ilen,lval,midval,uval,delta,eps;

  gma = (double *)R_alloc(*n, sizeof(double));
  e1 = (double *)R_alloc(*n, sizeof(double));
  e2 = (double *)R_alloc(*n, sizeof(double));
  v = (double *)R_alloc(*n, sizeof(double));
  jc = (double *)R_alloc(*n, sizeof(double));
  dvec = (double *)R_alloc(*n, sizeof(double));
  
  for(i=0;i<*n;i++) {
    datam1[i] = (datam1[i] - loc1[i]) / *scale1;
    datam2[i] = (datam2[i] - loc2[i]) / *scale2;
    if(*shape1 == 0) 
        datam1[i] = -datam1[i];
    else {
      datam1[i] = 1 + *shape1 * datam1[i];
      if(datam1[i] <= 0) {
        *dns = 1e6;
        return;
      }
      else datam1[i] = -1 / *shape1 * log(datam1[i]);
    }
    if(*shape2 == 0) 
      datam2[i] = -datam2[i];
    else {
      datam2[i] = 1 + *shape2 * datam2[i];
      if(datam2[i] <= 0) {
        *dns = 1e6;
        return;
      }
      else datam2[i] = -1 / *shape2 * log(datam2[i]);
    }
  }

  delta = eps = R_pow(DOUBLE_EPS, 0.5);
  for(i=0;i<*n;i++) {
    llim = 0;
    ulim = ilen = 1;
    lval = (1 - *alpha) * exp(datam1[i]);
    uval = (*beta - 1) * exp(datam2[i]);
    if(!(sign(lval) != sign(uval))) 
      error("values at end points are not of opposite sign");
    for(j=0;j<DOUBLE_DIGITS;j++) {
      ilen = ilen/2;
      midpt = llim + ilen;
      midval = (1 - *alpha) * exp(datam1[i]) * R_pow(1 - midpt, *beta) - 
               (1 - *beta) * exp(datam2[i]) * R_pow(midpt, *alpha);
      if(fabs(midval) < eps || fabs(ilen) < delta) 
        break;
      if(sign(lval) != sign(midval)) {
        ulim = midpt;
        uval = midval;
      }
      else {
        llim = midpt;
        lval = midval;
      }
    if(j == DOUBLE_DIGITS-1) 
      error("numerical problem in root finding algorithm");
    }
    gma[i] = midpt;
  }

  for(i=0;i<*n;i++) {
    v[i] = exp((1 - *alpha) * log(gma[i]) + datam1[i]) + 
      exp((1 - *beta) * log(1-gma[i]) + datam2[i]);
    jc[i] = (1 + *shape1) * datam1[i] + (1 + *shape2) * datam2[i] -
      log(*scale1 * *scale2);
    e1[i] = exp((1 - *alpha) * log(gma[i]) + (1 - *beta) * log(1 - gma[i]));
    e2[i] = exp(log(1 - *alpha) + log(*beta) + (*beta - 1) * log(1-gma[i]) +
            datam1[i]) + 
            exp(log(1 - *beta) + log(*alpha) + (*alpha - 1) * log(gma[i]) +
            datam2[i]);
    if(si[i] == 0) dvec[i] = log(e1[i]) - v[i] + jc[i];
    else if(si[i] == 1) 
      dvec[i] = log((1 - *alpha) * (1 - *beta)/e2[i]) - v[i] + jc[i];
    else dvec[i] = log(e1[i] + (1 - *alpha) * (1 - *beta) / e2[i]) - v[i] + jc[i];
  }

  if(*split > 0.5) {
    for(i=0;i<*n;i++) dns[i] = - dvec[i];
  }
  else {
    for(i=0;i<*n;i++) *dns = *dns - dvec[i];
  }
}

void nlbvnegbilog(double *datam1, double *datam2, int *n, int *si, double *alpha,
	         double *beta, double *loc1, double *scale1, double *shape1, 
                 double *loc2, double *scale2, double *shape2, int *split, 
                 double *dns)
{
  int i,j;
  double *e1,*e2,*e3,*v,*jc,*dvec,*gma;
  double llim,midpt,ulim,ilen,lval,midval,uval,delta,eps;

  gma = (double *)R_alloc(*n, sizeof(double));
  e1 = (double *)R_alloc(*n, sizeof(double));
  e2 = (double *)R_alloc(*n, sizeof(double));
  e3 = (double *)R_alloc(*n, sizeof(double));
  v = (double *)R_alloc(*n, sizeof(double));
  jc = (double *)R_alloc(*n, sizeof(double));
  dvec = (double *)R_alloc(*n, sizeof(double));
  
  for(i=0;i<*n;i++) {
    datam1[i] = (datam1[i] - loc1[i]) / *scale1;
    datam2[i] = (datam2[i] - loc2[i]) / *scale2;
    if(*shape1 == 0) 
        datam1[i] = -datam1[i];
    else {
      datam1[i] = 1 + *shape1 * datam1[i];
      if(datam1[i] <= 0) {
        *dns = 1e6;
        return;
      }
      else datam1[i] = -1 / *shape1 * log(datam1[i]);
    }
    if(*shape2 == 0) 
      datam2[i] = -datam2[i];
    else {
      datam2[i] = 1 + *shape2 * datam2[i];
      if(datam2[i] <= 0) {
        *dns = 1e6;
        return;
      }
      else datam2[i] = -1 / *shape2 * log(datam2[i]);
    }
  }

  delta = eps = R_pow(DOUBLE_EPS, 0.5);
  for(i=0;i<*n;i++) {
    llim = 0;
    ulim = ilen = 1;
    uval = (1 + *alpha) * exp(datam1[i]);
    lval = - (1 + *beta) * exp(datam2[i]);
    if(!(sign(lval) != sign(uval))) 
      error("values at end points are not of opposite sign");
    for(j=0;j<DOUBLE_DIGITS;j++) {
      ilen = ilen/2;
      midpt = llim + ilen;
      midval = (1 + *alpha) * exp(datam1[i]) * R_pow(midpt, *alpha) - 
               (1 + *beta) * exp(datam2[i]) * R_pow(1 - midpt, *beta);
      if(fabs(midval) < eps || fabs(ilen) < delta) 
        break;
      if(sign(lval) != sign(midval)) {
        ulim = midpt;
        uval = midval;
      }
      else {
        llim = midpt;
        lval = midval;
      }
    if(j == DOUBLE_DIGITS-1) 
      error("numerical problem in root finding algorithm");
    }
    gma[i] = midpt;
  }

  for(i=0;i<*n;i++) {
    v[i] = exp(datam1[i]) + exp(datam2[i]) - 
     exp((1 + *alpha) * log(gma[i]) + datam1[i]) - 
     exp((1 + *beta) * log(1-gma[i]) + datam2[i]);
    jc[i] = (1 + *shape1) * datam1[i] + (1 + *shape2) * datam2[i] -
      log(*scale1 * *scale2);
    e1[i] = (1 - R_pow(gma[i], 1 + *alpha)) * 
      (1 - R_pow(1 - gma[i], 1 + *beta));
    e2[i] = exp(log(1 + *alpha) + log(1 + *beta) + 
      *alpha * log(gma[i]) + *beta * log(1 - gma[i]));
    e3[i] = exp(log(1 + *alpha) + log(*alpha) + (*alpha - 1) * log(gma[i]) +
            datam1[i]) + 
            exp(log(1 + *beta) + log(*beta) + (*beta - 1) * log(1-gma[i]) +
            datam2[i]);
    if(si[i] == 0) dvec[i] = log(e1[i]) - v[i] + jc[i];
    else if(si[i] == 1) 
      dvec[i] = dvec[i] = log(e2[i] / e3[i]) - v[i] + jc[i];
    else dvec[i] = log(e1[i] + e2[i] / e3[i]) - v[i] + jc[i];
  }
  
  if(*split > 0.5) {
    for(i=0;i<*n;i++) dns[i] = - dvec[i];
  }
  else {
    for(i=0;i<*n;i++) *dns = *dns - dvec[i];
  }
}

void nlbvct(double *datam1, double *datam2, int *n, int *si, double *alpha,
	    double *beta, double *loc1, double *scale1, double *shape1, 
            double *loc2, double *scale2, double *shape2, int *split, 
            double *dns)
{
  int i;
  double *e1,*e2,*u,*v,*jc,*dvec;
  double c;

  e1 = (double *)R_alloc(*n, sizeof(double));
  e2 = (double *)R_alloc(*n, sizeof(double));
  u = (double *)R_alloc(*n, sizeof(double));
  v = (double *)R_alloc(*n, sizeof(double));
  jc = (double *)R_alloc(*n, sizeof(double));
  dvec = (double *)R_alloc(*n, sizeof(double));
  
  for(i=0;i<*n;i++) {
    datam1[i] = (datam1[i] - loc1[i]) / *scale1;
    datam2[i] = (datam2[i] - loc2[i]) / *scale2;
    if(*shape1 == 0) 
        datam1[i] = -datam1[i];
    else {
      datam1[i] = 1 + *shape1 * datam1[i];
      if(datam1[i] <= 0) {
        *dns = 1e6;
        return;
      }
      else datam1[i] = -1 / *shape1 * log(datam1[i]);
    }
    if(*shape2 == 0) 
      datam2[i] = -datam2[i];
    else {
      datam2[i] = 1 + *shape2 * datam2[i];
      if(datam2[i] <= 0) {
        *dns = 1e6;
        return;
      }
      else datam2[i] = -1 / *shape2 * log(datam2[i]);
    }
  }

  c = *alpha * *beta / (*alpha + *beta + 1);
  for(i=0;i<*n;i++) {
    u[i] = (*alpha * exp(datam2[i])) / 
      (*alpha * exp(datam2[i]) + *beta * exp(datam1[i]));
    v[i] = exp(datam2[i]) * pbeta(u[i], *alpha, *beta + 1, 1, 0) + 
      exp(datam1[i]) * pbeta(u[i], *alpha + 1, *beta, 0, 0);
    jc[i] = (1 + *shape1) * datam1[i] + (1 + *shape2) * datam2[i] -
      log(*scale1 * *scale2);
    e1[i] = pbeta(u[i], *alpha, *beta + 1, 1, 0) *
          pbeta(u[i], *alpha + 1, *beta, 0, 0);
    e2[i] = dbeta(u[i], *alpha + 1, *beta + 1, 0) /
          (*alpha * exp(datam2[i]) + *beta * exp(datam1[i]));
    if(si[i] == 0) dvec[i] = log(e1[i]) - v[i] + jc[i];
    else if(si[i] == 1) dvec[i] = log(c * e2[i]) - v[i] + jc[i];
    else dvec[i] = log(e1[i] + c * e2[i]) - v[i] + jc[i];
  }
  
  if(*split > 0.5) {
    for(i=0;i<*n;i++) dns[i] = - dvec[i];
  }
  else {
    for(i=0;i<*n;i++) *dns = *dns - dvec[i];
  }
}

void nlbvamix(double *datam1, double *datam2, int *n, int *si, double *alpha,
	     double *beta, double *loc1, double *scale1, double *shape1, 
             double *loc2, double *scale2, double *shape2, int *split, 
             double *dns)
{
  int i;
  double *v,*v1,*v2,*v12,*u,*u1,*u2,*jc,*dvec;
  double apb;

  v1 = (double *)R_alloc(*n, sizeof(double));
  v2 = (double *)R_alloc(*n, sizeof(double));
  v12 = (double *)R_alloc(*n, sizeof(double));
  u = (double *)R_alloc(*n, sizeof(double));
  u1 = (double *)R_alloc(*n, sizeof(double));
  u2 = (double *)R_alloc(*n, sizeof(double));
  v = (double *)R_alloc(*n, sizeof(double));
  jc = (double *)R_alloc(*n, sizeof(double));
  dvec = (double *)R_alloc(*n, sizeof(double));
  
  for(i=0;i<*n;i++) {
    datam1[i] = (datam1[i] - loc1[i]) / *scale1;
    datam2[i] = (datam2[i] - loc2[i]) / *scale2;
    if(*shape1 == 0) 
        datam1[i] = -datam1[i];
    else {
      datam1[i] = 1 + *shape1 * datam1[i];
      if(datam1[i] <= 0) {
        *dns = 1e6;
        return;
      }
      else datam1[i] = -1 / *shape1 * log(datam1[i]);
    }
    if(*shape2 == 0) 
      datam2[i] = -datam2[i];
    else {
      datam2[i] = 1 + *shape2 * datam2[i];
      if(datam2[i] <= 0) {
        *dns = 1e6;
        return;
      }
      else datam2[i] = -1 / *shape2 * log(datam2[i]);
    }
  }

  apb = *alpha + *beta;
  for(i=0;i<*n;i++) {
    jc[i] = (1 + *shape1) * datam1[i] + (1 + *shape2) * datam2[i] -
      log(*scale1 * *scale2);
    u[i] = exp(datam1[i]) + exp(datam2[i]);
    u1[i] = exp(datam1[i])/u[i];
    u2[i] = exp(datam2[i])/u[i];
    v[i] = u[i] - exp(datam1[i]) * (apb - *alpha * u1[i] - 
      *beta * u1[i] * u1[i]);
    v1[i] = 1 - *alpha * u2[i] * u2[i] - *beta * (3 * u2[i]*u2[i] -
      2 * u2[i]*u2[i]*u2[i]);
    v2[i] = 1 - *alpha * u1[i]*u1[i] - 2 * *beta * u1[i]*u1[i]*u1[i];
    v12[i] = (-2 * *alpha * u1[i] * u2[i] - 6 * *beta * u1[i]*u1[i] * 
      u2[i]) / u[i];
    if(si[i] == 0) dvec[i] = log(v1[i] * v2[i]) - v[i] + jc[i];
    else if(si[i] == 1) dvec[i] = log(- v12[i]) - v[i] + jc[i];
    else dvec[i] = log(v1[i] * v2[i] - v12[i]) - v[i] + jc[i];
  }
  
  if(*split > 0.5) {
    for(i=0;i<*n;i++) dns[i] = - dvec[i];
  }
  else {
    for(i=0;i<*n;i++) *dns = *dns - dvec[i];
  }
}
 

void nslmvalog(double *data, int *n, int *d, double *deps, double *thetas, 
    double *mpar, double *psrvs, int *q, int *nslocid, double *nsloc, 
    int *depindx, int *thetaindx, double *dns)
{
  int i,j,k,l,dd,nn,qq,niinb,niinbm,ndepp,nthetap,nmp;
  double iterm1, iterm2, term1, term2, eps;
  double thetasum, psrv, repdens;
  double dep, theta, loc;
  double *tdata, *dvec;
  int tmp1, tmp2;
  
  dd = *d; nn = *n; qq = *q;
  eps = R_pow(DOUBLE_EPS, 0.3);
  ndepp = R_pow(2, dd) - 1 - dd; 
  niinb = R_pow(2, dd - 1);
  nthetap = dd * (niinb - 1);
  niinbm = niinb-1;
  if(*nslocid) nmp = 4;
  else nmp = 3;
  *dns = 0;
  
  tdata = (double *)Calloc(nn * dd * sizeof(double), double);
  dvec = (double *)Calloc(nn * sizeof(double), double);

  for(i=0;i<nn;i++) dvec[i] = 0;

  for(i=0;i<nn;i++) {

    for(l=0;l<qq;l++) {

      repdens = 0;
      for(j=0;j<dd;j++) {

        if(!ISNA(data[i*dd+j])) {
          if(*nslocid) loc = mpar[4*j] + mpar[4*j+3] * nsloc[i];
          else loc = mpar[3*j]; 
          tdata[i*dd+j] = (data[i*dd+j] - loc) / mpar[nmp*j+1];
          if(fabs(mpar[nmp*j+2]) <= eps) tdata[i*dd+j] = exp(tdata[i*dd+j]);
          else {
            tdata[i*dd+j] = 1 + mpar[nmp*j+2] * tdata[i*dd+j];
            if(tdata[i*dd+j] <= 0) {
              *dns = 1e6;
              return;
	    }
            tdata[i*dd+j] = R_pow(tdata[i*dd+j], 1 / mpar[nmp*j+2]); 
          }

          thetasum = iterm1 = iterm2 = 0; 
          for(k=0;k<niinbm;k++) {
            tmp1 = depindx[j*niinbm + k];
            tmp2 = thetaindx[j*niinbm + k];
            dep = deps[tmp1];
          
            if(dep < 0.2) {
              *dns = 1e6;
              return;
	    }	  
            theta = thetas[tmp2];
            psrv = psrvs[tmp1 + ndepp * (l + i * qq)];
            term1 = psrv * R_pow(theta/tdata[i*dd+j], 1/dep);
            term2 = term1/dep;
            thetasum = thetasum + theta;
            iterm1 = iterm1 + term1;
            iterm2 = iterm2 + term2;    
          }
          if(thetasum > 1) {
            *dns = 1e6;
            return;
	  }
          else {
            iterm1 = iterm1 + (1-thetasum)/tdata[i*dd+j];
            iterm2 = iterm2 + (1-thetasum)/tdata[i*dd+j]; 
          }
          repdens = repdens + log(iterm2) - iterm1 - log(mpar[nmp*j+1]) -
            mpar[nmp*j+2] * log(tdata[i*dd+j]);
        }
        else tdata[i*dd+j] = NA_REAL;
      }
      dvec[i] = dvec[i] + exp(repdens);
    }
    dvec[i] = log(dvec[i]) - log(qq);
  }
 
  for(i=0;i<nn;i++) {
    *dns = *dns - dvec[i];
  }

  /*Rprintf("%f\n",*dns);
    error("stop");*/
  Free(dvec); Free(tdata);
  if(!R_FINITE(*dns) || R_IsNaN(*dns)) error("density is NaN/Inf");
}
