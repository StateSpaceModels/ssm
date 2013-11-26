//Special thanks to Aaron King for this template!

#include <pomp.h>
#include <R_ext/Rdynload.h>

#define r0	(p[parindex[0]])
#define v	(p[parindex[1]])
#define rep	(p[parindex[2]])
#define S_0	(p[parindex[3]])
#define I_0	(p[parindex[4]])
#define R_0	(p[parindex[5]])

#define S	(x[stateindex[0]])
#define I	(x[stateindex[1]])
#define R	(x[stateindex[2]])
#define Inc	(x[stateindex[3]])

#define cases	(y[obsindex[0]])

#define DS	(f[stateindex[0]])
#define DI	(f[stateindex[1]])
#define DR	(f[stateindex[2]])
#define DInc	(f[stateindex[3]])


#define Tr0	(pt[parindex[0]])
#define Tv	(pt[parindex[1]])
#define Trep	(pt[parindex[2]])
#define TS_0	(pt[parindex[3]])
#define TI_0	(pt[parindex[4]])
#define TR_0	(pt[parindex[5]])


#define lik	(like[0])


double logit_ab(double x, double a, double b)
{
    if (a == b)
        return x; // nothing will happen in the transformed space for x, so no need to transform it
    else{
        double ratio = (x-a)/(b-x);
        if(ratio < 1e-17){
            ratio = 1e-17;
        } else if(ratio > (1.0/1e-17)) {
            ratio = 1.0/1e-17;
        }
        return log(ratio);
    }
}


double inv_logit_ab(double x, double a, double b)
{
    if (a == b) {
        return x ;
    } else {
        if (x < 0) {
            return (b*exp(x)+a)/(1.0+exp(x));
        } else {
            return (b+a*exp(-x))/(1.0+exp(-x));
        }
    }
}


void SIR_par_trans (double *pt, double *p, int *parindex)
{

  Tr0 = inv_logit_ab(r0, 15, 35);
  Tv = exp(v); 
}

void SIR_par_untrans (double *pt, double *p, int *parindex)
{
  Tr0 = logit_ab(r0, 15, 35);
  Tv = log(v);  
}


void SIR_rmeasure (double *y, double *x, double *p, 
		   int *obsindex, int *stateindex, int *parindex, int *covindex, 
		   int ncovars, double *covars, double t)
{
    double mu = rep*Inc;
    double sd = sqrt((1-rep)*rep*Inc);

    double yobs = rnorm(mu, sd);
    if(yobs>0){
	cases = yobs;
    } else {
	cases = 0.0;
    }        

}


void SIR_dmeasure (double *like, double *y, double *x, double *p, int give_log, 
		   int *obsindex, int *stateindex, int *parindex, int *covindex, 
		   int ncovars, double *covars, double t)
{

    double mu = rep*Inc;
    double sd = sqrt((1-rep)*rep*Inc);

    if (cases > 0.0) {
        lik = pnorm(cases + 0.5, mu, sd, 1 , give_log) - pnorm(cases - 0.5, mu, sd, 1 , give_log);
    } else {
        lik = pnorm(cases + 0.5, mu, sd, 1 , give_log);
    }

}


void SIR_stepfn (double *x, const double *p, 
		 const int *stateindex, const int *parindex, const int *covindex, 
		 int covdim, const double *covars, double t, double dt)
{
  void (*reulermultinom)(int,double,double*,double,double*);
  reulermultinom = (void (*)(int,double,double*,double,double*)) R_GetCCallable("pomp","reulermultinom");

  double rate[2];		// transition rates
  double trans[2];		// transition numbers

  double myv = 1.0/(v/7.0);

  // compute the transition rates
  rate[0] = r0*myv*I/(S+I+R); 
  rate[1] = myv;			

  // compute the transition numbers
  reulermultinom(1, S, &rate[0], dt, &trans[0]);
  reulermultinom(1, I, &rate[1], dt, &trans[1]);

  // balance the equations
  S += -trans[0];
  I += trans[0]-trans[1];
  R += trans[1];
  Inc += trans[0];
 
}


void SIR_skelfn (double *f, double *x, double *p, 
		 int *stateindex, int *parindex, int *covindex, 
		 int ncovars, double *covars, double t)
{

  double rate[2];
  double term[2];

  double myv = 1.0/(v/7.0);

  // compute the transition rates
  rate[0] = r0*myv*I/(S+I+R); 
  rate[1] = myv;			
  
  // compute the several terms
  term[0] = rate[0]*S;
  term[1] = rate[1]*I;

  // assemble the differential equations
  DS = -term[0];
  DI = term[0]-term[1];
  DR = term[1];
  DInc = term[0]; 
}

#undef beta
#undef v
#undef rep
#undef S_0
#undef I_0
#undef R_0

#undef S
#undef I
#undef R
#undef Inc

#undef cases

#undef DS
#undef DI
#undef DR
#undef DInc


#undef Tbeta
#undef Tv
#undef Trep
#undef TS_0
#undef TI_0
#undef TR_0

#undef lik
