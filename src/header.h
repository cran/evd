#include <R.h>
#include <Rmath.h>

#define RANDIN GetRNGstate()
#define RANDOUT PutRNGstate()
#define UNIF unif_rand()
#define EXP exp_rand()

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

/* from fit.c */

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







