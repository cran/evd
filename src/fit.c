#include <R.h>
#include <Rmath.h>

void nlgev(double *data, int *n, double *loc, double *scale, double *shape, 
           double *dns);
void nlbvalog(double *datam1, double *datam2, int *n, double *dep,
	      double *asy1, double *asy2, double *loc1, double *scale1, 
              double *shape1, double *loc2, double *scale2, double *shape2, 
              double *dns);
void nlbvlog(double *datam1, double *datam2, int *n, double *dep, 
	     double *loc1, double *scale1, double *shape1, double *loc2, 
             double *scale2, double *shape2,  double *dns);
void nlbvhr(double *datam1, double *datam2, int *n, double *dep, 
            double *loc1, double *scale1, double *shape1, double *loc2, 
            double *scale2, double *shape2,  double *dns);
void nlbvneglog(double *datam1, double *datam2, int *n, double *dep, 
                double *loc1, double *scale1, double *shape1, double *loc2, 
                double *scale2, double *shape2,  double *dns);
void nlbvaneglog(double *datam1, double *datam2, int *n, double *dep,
	         double *asy1, double *asy2, double *loc1, double *scale1, 
                 double *shape1, double *loc2, double *scale2, double *shape2, 
                 double *dns);
void nlbvbilog(double *datam1, double *datam2, int *n, double *alpha,
	      double *beta, double *loc1, double *scale1, double *shape1, 
              double *loc2, double *scale2, double *shape2, double *dns);
void nlbvnegbilog(double *datam1, double *datam2, int *n, double *alpha,
	          double *beta, double *loc1, double *scale1, double *shape1, 
		  double *loc2, double *scale2, double *shape2, double *dns);
void nlbvct(double *datam1, double *datam2, int *n, double *alpha,
	    double *beta, double *loc1, double *scale1, double *shape1, 
            double *loc2, double *scale2, double *shape2, double *dns);

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

void nlbvlog(double *datam1, double *datam2, int *n, double *dep, 
             double *loc1, double *scale1, double *shape1, double *loc2, 
             double *scale2, double *shape2,  double *dns)
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
    dvec[i] = dvec[i] + (1-2*idep)*log(z[i]) + log(idep-1+z[i]) - z[i];
  }

  for(i=0;i<*n;i++) 
    *dns = *dns - dvec[i];
}

void nlbvalog(double *datam1, double *datam2, int *n, double *dep,
	      double *asy1, double *asy2, double *loc1, double *scale1, 
              double *shape1, double *loc2, double *scale2, double *shape2, 
              double *dns)
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
      (1 - 2*idep) * log(z[i]) + log(idep-1+z[i]);
    dvec[i] = log(exp(c1) + exp(e3[i]) + exp(e4[i])) - v[i] + jc[i];
  }
  for(i=0;i<*n;i++) 
    *dns = *dns - dvec[i];
}

void nlbvhr(double *datam1, double *datam2, int *n, double *dep, 
            double *loc1, double *scale1, double *shape1, double *loc2, 
            double *scale2, double *shape2,  double *dns)
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
    dvec[i] = e1[i] * e2[i] + *dep * e3[i] / 2;
    jc[i] = *shape1 * datam1[i] + *shape2 * datam2[i] -
      log(*scale1 * *scale2);
    dvec[i] = log(dvec[i]) + jc[i] - v[i];
  }

  for(i=0;i<*n;i++) 
    *dns = *dns - dvec[i];
}

void nlbvneglog(double *datam1, double *datam2, int *n, double *dep, 
                double *loc1, double *scale1, double *shape1, double *loc2, 
                double *scale2, double *shape2,  double *dns)
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
      (1 + 2 * *dep) * log(z[i]) + log(1 + *dep + z[i]);
    /*if(*dep > 0.95 && i == 9) {
      Rprintf("%f ",datam1[9]);
      Rprintf("%f ",datam2[9]);
      Rprintf("%f ",z[9]);
      Rprintf("%f ",e1[9]);
      Rprintf("%f ",e2[9]);
      Rprintf("\n");
      }*/
    dvec[i] = log(1 - exp(e1[i]) + exp(e2[i])) - v[i] + jc[i];
  }

  for(i=0;i<*n;i++) 
    *dns = *dns - dvec[i];
}

void nlbvaneglog(double *datam1, double *datam2, int *n, double *dep,
	         double *asy1, double *asy2, double *loc1, double *scale1, 
                 double *shape1, double *loc2, double *scale2, double *shape2, 
                 double *dns)
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
    e4[i] = e1[i] + e2[i] + (1 + 2 * *dep) * log(z[i]) + log(1 + *dep + z[i]);
    dvec[i] = log(1 - exp(e3[i]) + exp(e4[i])) - v[i] + jc[i];
  }

  for(i=0;i<*n;i++) 
    *dns = *dns - dvec[i];
}

void nlbvbilog(double *datam1, double *datam2, int *n, double *alpha,
	      double *beta, double *loc1, double *scale1, double *shape1, 
              double *loc2, double *scale2, double *shape2, double *dns)
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
    dvec[i] = log(e1[i] + (1 - *alpha) * (1 - *beta) / e2[i]) - v[i] + jc[i];
  }
  for(i=0;i<*n;i++) 
    *dns = *dns - dvec[i];
}

void nlbvnegbilog(double *datam1, double *datam2, int *n, double *alpha,
	         double *beta, double *loc1, double *scale1, double *shape1, 
                 double *loc2, double *scale2, double *shape2, double *dns)
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
    dvec[i] = log(e1[i] + e2[i] / e3[i]) - v[i] + jc[i];
  }
  for(i=0;i<*n;i++) 
    *dns = *dns - dvec[i];
}

void nlbvct(double *datam1, double *datam2, int *n, double *alpha,
	    double *beta, double *loc1, double *scale1, double *shape1, 
            double *loc2, double *scale2, double *shape2, double *dns)
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
    dvec[i] = log(e1[i] + c * e2[i]) - v[i] + jc[i];
  }
  for(i=0;i<*n;i++) 
    *dns = *dns - dvec[i];
}
 
