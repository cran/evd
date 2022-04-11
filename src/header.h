#include <R.h>
#include <Rmath.h>

#define RANDIN GetRNGstate()
#define RANDOUT PutRNGstate()
#define UNIF unif_rand()
#define EXP exp_rand()

/* from pot.c */

void nlgpd(double *data, int *n, double *loc, double *scale, 
	   double *shape, double *dns);
void nlpp(double *exceed, int *nhigh, double *loc, double *scale, 
          double *shape, double *thresh, double *nop, double *dns);
void clusters(double *high, double *high2, int *n, int *r, 
	      int *rlow, double *clstrs);

/* from ccop.c */

double ccbvlog(double m1, double m2, double oldm1, double dep);
double ccbvalog(double m1, double m2, double oldm1, double dep, double asy1, 
                double asy2);
double ccbvhr(double m1, double m2, double oldm1, double dep);
double ccbvneglog(double m1, double m2, double oldm1, double dep);
double ccbvaneglog(double m1, double m2, double oldm1, double dep, 
                   double asy1, double asy2);
double ccbvbilog(double m1, double m2, double oldm1, double alpha, 
                double beta);
double ccbvnegbilog(double m1, double m2, double oldm1, double alpha, 
                   double beta);
double ccbvct(double m1, double m2, double oldm1, double alpha, double beta);
double ccbvamix(double m1, double m2, double oldm1, double alpha, double beta);
void ccop(double *m1, double *m2, int *cnd, double *dep, double *asy1, 
          double *asy2, double *alpha, double *beta, int *n, int *model, 
          double *ccop);

/* from sim.c */

void rbvlog_shi(int *n, double *alpha, double *sim);
void rbvalog_shi(int *n, double *alpha, double *asy, double *sim);
void rmvlog_tawn(int *n, int *d, double *alpha, double *sim);
void rmvalog_tawn(int *n, int *d, int *nb, double *alpha, double *asy, 
                  double *sim);
double rpstable(double cexp);
double maximum_n(int n, double *x);
void rbvlog(int *n, double *dep, double *sim);
void rbvalog(int *n, double *dep, double *asy, double *sim);
void rbvhr(int *n, double *dep, double *sim);
void rbvneglog(int *n, double *dep, double *sim);
void rbvaneglog(int *n, double *dep, double *asy, double *sim);
void rbvbilog(int *n, double *alpha, double *beta, double *sim);
void rbvnegbilog(int *n, double *alpha, double *beta, double *sim);
void rbvct(int *n, double *alpha, double *beta, double *sim);
void rbvamix(int *n, double *alpha, double *beta, double *sim);

/* from fit.c */

void nlgev(double *data, int *n, double *loc, double *scale, double *shape, 
           double *dns);		   
void nlgumbelx(double *data, int *n, double *loc1, double *scale1, double *loc2, double *scale2, 
           double *dns);
void nlbvalog(double *datam1, double *datam2, int *n, int *si, double *dep,
	      double *asy1, double *asy2, double *loc1, double *scale1, 
              double *shape1, double *loc2, double *scale2, double *shape2, 
              int *split, double *dns);
void nlbvlog(double *datam1, double *datam2, int *n, int *si, double *dep, 
	     double *loc1, double *scale1, double *shape1, double *loc2, 
             double *scale2, double *shape2, int *split, double *dns);
void nlbvhr(double *datam1, double *datam2, int *n, int *si, double *dep, 
            double *loc1, double *scale1, double *shape1, double *loc2, 
            double *scale2, double *shape2, int *split, double *dns);
void nlbvneglog(double *datam1, double *datam2, int *n, int *si, double *dep, 
                double *loc1, double *scale1, double *shape1, double *loc2, 
                double *scale2, double *shape2, int *split, double *dns);
void nlbvaneglog(double *datam1, double *datam2, int *n, int *si, double *dep,
	         double *asy1, double *asy2, double *loc1, double *scale1, 
                 double *shape1, double *loc2, double *scale2, double *shape2, 
                 int *split, double *dns);
void nlbvbilog(double *datam1, double *datam2, int *n, int *si, double *alpha,
	      double *beta, double *loc1, double *scale1, double *shape1, 
              double *loc2, double *scale2, double *shape2, int *split, 
              double *dns);
void nlbvnegbilog(double *datam1, double *datam2, int *n, int *si, double *alpha,
	          double *beta, double *loc1, double *scale1, double *shape1, 
		  double *loc2, double *scale2, double *shape2, int *split, 
                  double *dns);
void nlbvct(double *datam1, double *datam2, int *n, int *si, double *alpha,
	    double *beta, double *loc1, double *scale1, double *shape1, 
            double *loc2, double *scale2, double *shape2, int *split, 
            double *dns);
void nlbvamix(double *datam1, double *datam2, int *n, int *si, double *alpha,
	      double *beta, double *loc1, double *scale1, double *shape1, 
              double *loc2, double *scale2, double *shape2, int *split, 
	      double *dns);

void nslmvalog(double *data, int *n, int *d, double *deps, double *thetas, 
               double *mpar, double *psrvs, int *q, int *nslocid, double *nsloc, 
	       int *depindx, int *thetaindx, double *dns);

/* from bvpot.c (censored) */

void nllbvclog(double *data1, double *data2, int *nn, int *n, double *thid, 
              double *lambda, double *dep, double *scale1, double *shape1, 
              double *scale2, double *shape2, double *dns);
void nllbvcbilog(double *data1, double *data2, int *nn, int *n, double *thid, 
                 double *lambda, double *alpha, double *beta, 
                 double *scale1, double *shape1, double *scale2, 
                 double *shape2, double *dns);
void nllbvcalog(double *data1, double *data2, int *nn, int *n, double *thid, 
                double *lambda, double *dep, double *asy1, double *asy2, 
                double *scale1, double *shape1, double *scale2, 
                double *shape2, double *dns);
void nllbvcneglog(double *data1, double *data2, int *nn, int *n, double *thid, 
                  double *lambda, double *dep, double *scale1, 
                  double *shape1, double *scale2, double *shape2, 
                  double *dns);
void nllbvcnegbilog(double *data1, double *data2, int *nn, int *n, 
                    double *thid, double *lambda, double *alpha, double *beta, 
                    double *scale1, double *shape1, double *scale2, 
                    double *shape2, double *dns);
void nllbvcaneglog(double *data1, double *data2, int *nn, int *n, 
                   double *thid, double *lambda, double *dep, double *asy1, 
                   double *asy2, double *scale1, double *shape1, 
                   double *scale2, double *shape2, double *dns);
void nllbvcct(double *data1, double *data2, int *nn, int *n, double *thid, 
              double *lambda, double *alpha, double *beta, double *scale1, 
              double *shape1, double *scale2, double *shape2, double *dns);
void nllbvchr(double *data1, double *data2, int *nn, int *n, double *thid, 
              double *lambda, double *dep, double *scale1, double *shape1, 
              double *scale2, double *shape2, double *dns);
void nllbvcamix(double *data1, double *data2, int *nn, int *n, double *thid, 
                double *lambda, double *alpha, double *beta, double *scale1, 
                double *shape1, double *scale2, double *shape2, double *dns);


/* from bvpot.c (poisson) */

void nllbvplog(double *data1, double *data2, int *nn, 
               double *thid, double *r1, double *r2, double *p, 
               double *dep, double *scale1, double *shape1, 
               double *scale2, double *shape2, double *dns);
void nllbvpneglog(double *data1, double *data2, int *nn, 
                  double *thid, double *r1, double *r2, double *p, 
                  double *dep, double *scale1, double *shape1, 
                  double *scale2, double *shape2, double *dns);
void nllbvpct(double *data1, double *data2, int *nn, double *thid, 
              double *r1, double *r2, double *p, double *alpha, 
              double *beta, double *scale1, double *shape1, double *scale2, 
              double *shape2, double *dns);
void nllbvpbilog(double *data1, double *data2, int *nn, 
                 double *thid, double *r1, double *r2, double *p, 
                 double *alpha, double *beta, double *scale1, 
                 double *shape1, double *scale2, double *shape2, 
                 double *dns);
void nllbvpnegbilog(double *data1, double *data2, int *nn, 
                    double *thid, double *r1, double *r2, double *p, 
                    double *alpha, double *beta, double *scale1, 
                    double *shape1, double *scale2, double *shape2, 
                    double *dns);
void nllbvphr(double *data1, double *data2, int *nn, double *thid, 
              double *r1, double *r2, double *p, double *dep, 
              double *scale1, double *shape1, double *scale2, 
              double *shape2, double *dns);
