#include <R.h>
#include <Rmath.h>

#define RANDIN GetRNGstate()
#define RANDOUT PutRNGstate()
#define UNIF unif_rand()
#define EXP exp_rand()

void rbvlog_shi(int *n, double *alpha, double *sim);
void rbvalog_shi(int *n, double *alpha, double *asy, double *sim);
void rmvlog_tawn(int *n, int *d, double *alpha, double *sim);
void rmvalog_tawn(int *n, int *d, int *nb, double *alpha, double *asy, 
                  double *sim);
double rpstable(double cexp);
double maximum_n(int n, double *x);
void rbvhr(double *n, double *dep, double *sim);
double ccbvhr(double m1, double m2, double oldm1, double dep);
void rbvneglog(double *n, double *dep, double *sim);
double ccbvneglog(double m1, double m2, double oldm1, double dep);
void rbvaneglog(double *n, double *dep, double *asy, double *sim);
double ccbvaneglog(double m1, double m2, double oldm1, double dep, 
                   double asy1, double asy2);
void rbvbilog(double *n, double *alpha, double *beta, double *sim);
double ccbvbilog(double m1, double m2, double oldm1, double alpha, 
                double beta);
void rbvnegbilog(double *n, double *alpha, double *beta, double *sim);
double ccbvnegbilog(double m1, double m2, double oldm1, double alpha, 
                   double beta);
void rbvct(double *n, double *alpha, double *beta, double *sim);
double ccbvct(double m1, double m2, double oldm1, double alpha, double beta);


/* 
   All simulation functions produce standard Frechet margins 
   unless explicitly stated otherwise.
*/

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

/* produces uniform margins */
void rbvhr(double *n, double *dep, double *sim)
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

double ccbvhr(double m1, double m2, double oldm1, double dep)
{
  double tm1,tm2,v,idep,q,fval;

  tm1 = -log(m1);
  tm2 = -log(m2);
  idep = 1 / dep;
  v = tm2 * pnorm(idep + (log(tm2) - log(tm1)) * dep/2, 0, 1, 1, 0) +
    tm1 * pnorm(idep + (log(tm1) - log(tm2)) * dep/2, 0, 1, 1, 0);
  q = pnorm(idep + (log(tm2) - log(tm1)) * dep/2, 0, 1, 1, 0) *
    exp(-v) / m2 - oldm1;
  fval = pnorm(idep + (log(tm2) - log(tm1)) * dep/2, 0, 1, 1, 0) *
    exp(-v) / m2 - oldm1;
  return fval;
}

/* produces uniform margins */
void rbvneglog(double *n, double *dep, double *sim)
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

/* produces uniform margins */
void rbvaneglog(double *n, double *dep, double *asy, double *sim)
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

/* produces uniform margins */
void rbvbilog(double *n, double *alpha, double *beta, double *sim)
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

double ccbvbilog(double m1, double m2, double oldm1, double alpha, 
                double beta)
{
  int i;
  double tm1,tm2,v,fval;
  double delta,eps,llim,midpt,ulim,ilen,lval,midval,uval;

  tm1 = -log(m1);
  tm2 = -log(m2);

  delta = eps = R_pow(DOUBLE_EPS, 0.75);
  llim = 0;
  ulim = ilen = 1;
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
      ulim = midpt;
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

/* produces uniform margins */
void rbvnegbilog(double *n, double *alpha, double *beta, double *sim)
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

double ccbvnegbilog(double m1, double m2, double oldm1, double alpha, 
                   double beta)
{
  int i;
  double tm1,tm2,v,fval;
  double delta,eps,llim,midpt,ulim,ilen,lval,midval,uval;

  tm1 = -log(m1);
  tm2 = -log(m2);

  delta = eps = R_pow(DOUBLE_EPS, 0.75);
  llim = 0;
  ulim = ilen = 1;
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
      ulim = midpt;
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

/* produces uniform margins */
void rbvct(double *n, double *alpha, double *beta, double *sim)
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



