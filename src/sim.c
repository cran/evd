#include "header.h"


/* produces standard Frechet margins */
void rbvlog_shi(int *n, double *alpha, double *sim)
{
  double u,z;
  int i;
  
  RANDIN;
  for(i=0;i<*n;i++) 
  { 
    u = UNIF;
    if(UNIF < *alpha) z = EXP+EXP;
    else z = EXP;
    sim[2*i] = 1/(z * R_pow(u,*alpha));
    sim[2*i+1] = 1/(z * R_pow(1-u,*alpha));
  }
  RANDOUT;
}

/* produces standard Frechet margins */
void rbvalog_shi(int *n, double *alpha, double *asy, double *sim)
{
  double v1_1,v2_2,v1_12,v2_12,u,z;
  int i;
    
  RANDIN;

  if(*alpha == 1)
    for(i=0;i<2*(*n);i++) sim[i] = 1/EXP;
  else {
    for(i=0;i<*n;i++) 
    {
      v1_1 = (1-asy[0]) / EXP;
      v2_2 = (1-asy[1]) / EXP;
      u = UNIF;
      if(UNIF < *alpha) z = EXP+EXP;
      else z = EXP;
      v1_12 = asy[0] / (z * R_pow(u,*alpha));
      v2_12 = asy[1] / (z * R_pow(1-u,*alpha));
      sim[2*i] = fmax2(v1_1,v1_12); 
      sim[2*i+1] = fmax2(v2_2,v2_12);
    }
  }
  RANDOUT;
}

/* produces standard Frechet margins */
void rmvlog_tawn(int *n, int *d, double *alpha, double *sim)
{
  double s;
  int i,j;
    
  RANDIN;
  for(i=0;i<*n;i++) 
  {
    s = rpstable(*alpha);
    for(j=0;j<*d;j++) 
      sim[i*(*d)+j] = R_pow(s/EXP,*alpha);
  }
  RANDOUT;
}

/* produces standard Frechet margins */
void rmvalog_tawn(int *n, int *d, int *nb, double *alpha, double *asy, 
                  double *sim)
{
  double s;
  double *gevsim;
  double *maxsim;
  int i,j,k;

  gevsim = (double *)R_alloc((*nb)*(*d), sizeof(double));
  maxsim = (double *)R_alloc(*nb, sizeof(double));

  for(i=0;i<(*nb)*(*d);i++)  
    gevsim[i] = 0;
    
  RANDIN;
  for(i=0;i<*n;i++) 
  {
    for(j=0;j<*nb;j++) {
      if(alpha[j] != 1) 
        s = rpstable(alpha[j]);
      else s = 1;
      for(k=0;k<*d;k++) {
	if(asy[j*(*d) + k] != 0) 
	    gevsim[j*(*d) + k] = asy[j*(*d) + k] * R_pow((s/EXP),alpha[j]);
      }
    }
    for(j=0;j<*d;j++) {
      for(k=0;k<*nb;k++) 
        maxsim[k] = gevsim[k*(*d) + j];
      sim[i*(*d)+j] = maximum_n(*nb,maxsim);
    }
  }
  RANDOUT;
}

double rpstable(double cexp)
{
  double tcexp,u,w,a;

  if(cexp == 1) return 1;
  tcexp = 1-cexp;
  u = M_PI*UNIF;
  w = EXP; 
  a = sin(tcexp*u) * R_pow(sin(cexp*u),(cexp/tcexp)) / 
      R_pow(sin(u),(1/tcexp));
  return R_pow((a/w),(tcexp/cexp));
}

double maximum_n(int n, double *x)
{
  int i;

  for(i=1;i<n;i++)
    if(x[0] < x[i]) x[0] = x[i];
  return x[0];
}

/* produces uniform margins; needed for evmc */
void rbvlog(int *n, double *dep, double *sim)
{
  double delta,eps,llim,midpt,ulim,ilen,lval,midval,uval;
  int i,j;

  for(i=0;i<*n;i++) 
  {
    delta = eps = llim = R_pow(DOUBLE_EPS, 0.5);
    ulim = 1 - llim;
    ilen = 1;
    midpt = 0.5;
    lval = ccbvlog(llim, sim[2*i+1], sim[2*i+0], *dep);
    uval = ccbvlog(ulim, sim[2*i+1], sim[2*i+0], *dep);
    if(!(sign(lval) != sign(uval))) 
      error("values at end points are not of opposite sign");
    for(j=0;j<DOUBLE_DIGITS;j++) {
      ilen = ilen/2;
      midpt = llim + ilen;
      midval = ccbvlog(midpt, sim[2*i+1], sim[2*i+0], *dep);
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
    sim[2*i+0] = midpt;
  }
}

/* produces uniform margins; needed for evmc */
void rbvalog(int *n, double *dep, double *asy, double *sim)
{
  double delta,eps,llim,midpt,ulim,ilen,lval,midval,uval;
  int i,j;

  for(i=0;i<*n;i++) 
  {
    delta = eps = llim = R_pow(DOUBLE_EPS, 0.5);
    ulim = 1 - llim;
    ilen = 1;
    midpt = 0.5;
    lval = ccbvalog(llim, sim[2*i+1], sim[2*i+0], *dep,asy[0],asy[1]);
    uval = ccbvalog(ulim, sim[2*i+1], sim[2*i+0], *dep,asy[0],asy[1]);
    if(!(sign(lval) != sign(uval))) 
      error("values at end points are not of opposite sign");
    for(j=0;j<DOUBLE_DIGITS;j++) {
      ilen = ilen/2;
      midpt = llim + ilen;
      midval = ccbvalog(midpt, sim[2*i+1], sim[2*i+0], *dep,asy[0],asy[1]);
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
    sim[2*i+0] = midpt;
  }
}

/* produces uniform margins */
void rbvhr(int *n, double *dep, double *sim)
{
  double delta,eps,llim,midpt,ulim,ilen,lval,midval,uval;
  int i,j;

  for(i=0;i<*n;i++) 
  {
    delta = eps = llim = R_pow(DOUBLE_EPS, 0.5);
    ulim = 1 - llim;
    ilen = 1;
    midpt = 0.5;
    lval = ccbvhr(llim, sim[2*i+1], sim[2*i+0], *dep);
    uval = ccbvhr(ulim, sim[2*i+1], sim[2*i+0], *dep);
    if(!(sign(lval) != sign(uval))) 
      error("values at end points are not of opposite sign");
    for(j=0;j<DOUBLE_DIGITS;j++) {
      ilen = ilen/2;
      midpt = llim + ilen;
      midval = ccbvhr(midpt, sim[2*i+1], sim[2*i+0], *dep);
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
    sim[2*i+0] = midpt;
  }
}

/* produces uniform margins */
void rbvneglog(int *n, double *dep, double *sim)
{
  double delta,eps,llim,midpt,ulim,ilen,lval,midval,uval;
  int i,j;

  for(i=0;i<*n;i++) 
  {
    delta = eps = llim = R_pow(DOUBLE_EPS, 0.5);
    ulim = 1 - llim;
    ilen = 1;
    midpt = 0.5;
    lval = ccbvneglog(llim, sim[2*i+1], sim[2*i+0], *dep);
    uval = ccbvneglog(ulim, sim[2*i+1], sim[2*i+0], *dep);
    if(!(sign(lval) != sign(uval))) 
      error("values at end points are not of opposite sign");
    for(j=0;j<DOUBLE_DIGITS;j++) {
      ilen = ilen/2;
      midpt = llim + ilen;
      midval = ccbvneglog(midpt, sim[2*i+1], sim[2*i+0], *dep);
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
    sim[2*i+0] = midpt;
  }
}

/* produces uniform margins */
void rbvaneglog(int *n, double *dep, double *asy, double *sim)
{
  double delta,eps,llim,midpt,ulim,ilen,lval,midval,uval;
  int i,j;

  for(i=0;i<*n;i++) 
  {
    delta = eps = llim = R_pow(DOUBLE_EPS, 0.5);
    ulim = 1 - llim;
    ilen = 1;
    midpt = 0.5;
    lval = ccbvaneglog(llim, sim[2*i+1], sim[2*i+0], *dep,asy[0],asy[1]);
    uval = ccbvaneglog(ulim, sim[2*i+1], sim[2*i+0], *dep,asy[0],asy[1]);
    if(!(sign(lval) != sign(uval))) 
      error("values at end points are not of opposite sign");
    for(j=0;j<DOUBLE_DIGITS;j++) {
      ilen = ilen/2;
      midpt = llim + ilen;
      midval = ccbvaneglog(midpt, sim[2*i+1], sim[2*i+0], *dep,asy[0],asy[1]);
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
    sim[2*i+0] = midpt;
  }
}

/* produces uniform margins */
void rbvbilog(int *n, double *alpha, double *beta, double *sim)
{
  double delta,eps,llim,midpt,ulim,ilen,lval,midval,uval;
  int i,j;

  for(i=0;i<*n;i++) 
  {
    delta = eps = llim = R_pow(DOUBLE_EPS, 0.5);
    ulim = 1 - llim;
    ilen = 1;
    midpt = 0.5;
    lval = ccbvbilog(llim, sim[2*i+1], sim[2*i+0], *alpha, *beta);
    uval = ccbvbilog(ulim, sim[2*i+1], sim[2*i+0], *alpha, *beta);
    if(!(sign(lval) != sign(uval))) 
      error("values at end points are not of opposite sign");
    for(j=0;j<DOUBLE_DIGITS;j++) {
      ilen = ilen/2;
      midpt = llim + ilen;
      midval = ccbvbilog(midpt, sim[2*i+1], sim[2*i+0], *alpha, *beta);
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
    sim[2*i+0] = midpt;
  }
}

/* produces uniform margins */
void rbvnegbilog(int *n, double *alpha, double *beta, double *sim)
{
  double delta,eps,llim,midpt,ulim,ilen,lval,midval,uval;
  int i,j;

  for(i=0;i<*n;i++) 
  {
    delta = eps = llim = R_pow(DOUBLE_EPS, 0.5);
    ulim = 1 - llim;
    ilen = 1;
    midpt = 0.5;
    lval = ccbvnegbilog(llim, sim[2*i+1], sim[2*i+0], *alpha, *beta);
    uval = ccbvnegbilog(ulim, sim[2*i+1], sim[2*i+0], *alpha, *beta);
    if(!(sign(lval) != sign(uval))) 
      error("values at end points are not of opposite sign2");
    for(j=0;j<DOUBLE_DIGITS;j++) {
      ilen = ilen/2;
      midpt = llim + ilen;
      midval = ccbvnegbilog(midpt, sim[2*i+1], sim[2*i+0], *alpha, *beta);
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
    sim[2*i+0] = midpt;
  }
}

/* produces uniform margins */
void rbvct(int *n, double *alpha, double *beta, double *sim)
{
  double delta,eps,llim,midpt,ulim,ilen,lval,midval,uval;
  int i,j;

  for(i=0;i<*n;i++) 
  {
    delta = eps = llim = R_pow(DOUBLE_EPS, 0.5);
    ulim = 1 - llim;
    ilen = 1;
    midpt = 0.5;
    lval = ccbvct(llim, sim[2*i+1], sim[2*i+0], *alpha, *beta);
    uval = ccbvct(ulim, sim[2*i+1], sim[2*i+0], *alpha, *beta);
    if(!(sign(lval) != sign(uval))) 
      error("values at end points are not of opposite sign2");
    for(j=0;j<DOUBLE_DIGITS;j++) {
      ilen = ilen/2;
      midpt = llim + ilen;
      midval = ccbvct(midpt, sim[2*i+1], sim[2*i+0], *alpha, *beta);
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
    sim[2*i+0] = midpt;
  }
}













