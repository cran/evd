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

/* All simulation functions produce standard Frechet margins */

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

  gevsim = (double *)malloc(sizeof(double)*(*nb)*(*d));
  maxsim = (double *)malloc(sizeof(double)*(*nb));

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



















