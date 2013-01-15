#ifndef CASPER
#define CASPER

#include "dataframe.h"
#include "cstat.h"
#include <vector>
#include <math.h>

using namespace std;

class Casper
{
 public:
  Casper(Model* model, DataFrame* frame);

  // current model (set of variants) that tries to explain the data
  Model* model;
  // data
  DataFrame* frame;

  // EM to find posterior mode under the current model
  double* calculateMode();  //initialize to equal expression for all variants
  void calculateMode(double* pi);  //use pi as initial value, and return it updated with solution

  // Asymptotic standard error for mixture proportions (uses delta method)
  void asymptoticSE(double *se, double *mode, int n);

  // Normal approximation to posterior on logit re-parameterization
  void normapprox(double **S, double *mode, int n, int Sidx_ini=0);  //note: S is Hessian at the mode, i.e. inverse of covariance matrix
  void normapprox(double **S, double** G, double*** H, double* mode, double* thmode, int n, int Sidx_ini=0); //note: S is Hessian at the mode, i.e. inverse of covariance matrix

  // Independent proposal Metropolis-Hastings
  void IPMH(double *pi, double *paccept, double *integralIS, int niter, int burnin);  //stores sample in pi, proportion of accepted proposals in paccept
  void IPMH(double *pi, double *paccept, double *integralIS, int niter, int burnin, double *mode);  //same but uses pre-computed mode
  void IPMH(double *pi, double *paccept, double *integralIS, int niter, int burnin, double *mode, double **Sinv); //same using pre-computed mode & hessian

  // gives the integral given the current model and data
  double calculateIntegral(int method=1);  //uses mode=calculateMode() and n= model->count()
  double calculateIntegral(double* mode, int n, int method=1);  //do integral with pre-computed mode
  double LaplaceApprox(double *mode, int n);

  // evaluate the likelihood & prior
  double priorLikelihoodLn(double* pi);
  double priorLn(double* pi);
  double likelihoodLn(double* pi);

  // indicates whether the model can explain all fragments
  bool isValid();

  static double priorq;
  static int em_maxruns;
  static double em_tol;
  static int is_runs;

 private:

  map<Fragment*, map<Variant*, double> > mempprobs;
  map<Variant*, map<Fragment*, double> > memvprobs;

  static const double mh_gammah;

  map<Fragment*, double> fragdist(double* pi);

  void mlogit(double *theta, double* pi, int n);
  void milogit(double *pi, double* theta, int n);
  void vtGradG(double **G, double* th, int n);
  double vtGradLogdet(double** G, int n);
  void vtHess(double ***H, double* th, int n);
  double det(double** a, int n, bool *posdef);
};

#endif
